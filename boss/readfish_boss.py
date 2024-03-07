"""
wrapper around readfish functionality

- based on the `targets.py` entry-point of readfish
- instead of launching it as an entry-point, we run _cli_base.main()
- the only mod in that function is to return parser and args instead of sending it
- from targets.py: subclass Analysis to modify run(): this contains most changes
- main run() function of entry-point is only modified to use the AnalysisMod object
- also overwrite the make_decision function -> BossBits.make_decision_boss()
- that function originally lives in _config.py

for future changes check:
_cli_base.main()
targets.py (Analysis.run(), run())
_config.py.make_decision()


"""

# Core imports
from __future__ import annotations
import argparse
import logging
import time
from timeit import default_timer as timer
from pathlib import Path
from typing import Any

# Third party imports
from readfish.read_until.read_cache import AccumulatingCache
from readfish.read_until import ReadUntilClient
from minknow_api import protocol_service

# Library
from readfish._cli_args import DEVICE_BASE_ARGS, DEFAULT_UNBLOCK
from readfish._read_until_client import RUClient
from readfish._config import Action, Conf, make_decision, _Condition
from readfish._statistics import ReadfishStatistics
from readfish._utils import (
    get_device,
    send_message,
    ChunkTracker,
    Severity,
)
from readfish.plugins.abc import AlignerABC, CallerABC
from readfish.plugins.utils import Decision, PreviouslySentActionTracker, Result


# TODO custom import
from readfish.entry_points import targets
from _cli_base import main as main_args
import numpy as np
import sys
import os
import mappy


class AnalysisMod(targets.Analysis):

    def run(self):
        """Run the read until loop, in one continuous while loop."""

        # TODO: Swap this for a CSV record later
        self.conf.write_channels_toml(self.client.mk_run_dir)

        # TODO: This could still be passed through to the basecaller to prevent
        #       rebasecalling data that is already being unblocked or sequenced
        loop_counter = 0

        last_live_toml_mtime = 0
        self.logger.info("Starting main loop")
        mapper_description = self.mapper.describe(self.conf.regions, self.conf.barcodes)
        self.logger.info(mapper_description)
        send_message(self.client.connection, mapper_description, Severity.INFO)

        # boss - initialise parts needed for dynamic decisions
        boss = BossBits(conf=self.conf, logger=self.logger, mapper=self.mapper)


        while self.client.is_sequencing:
            t0 = timer()
            # Check if we have started readfish before PHASE_SEQUENCING,
            if self.wait_for_sequencing:
                time.sleep(self.throttle)
                continue
            # Set back to true for when we re-enter a non sequencing phase
            self.log_once_in_loop = True

            if self.readfish_started_during_sequencing and loop_counter == 0:
                self.logger.info(
                    "readfish started in PHASE_SEQUENCING. Fully sequencing first read from each channel."
                )
            if not self.mapper.initialised:
                self.logger.warning(
                    "readfish main loop started but mapper is not initialised. Please check your aligners plugin documentation."
                    "If you are using mappy or mappy-rs this is definitely an error. Please open an issue here - "
                    "https://github.com/LooseLab/readfish/issues"
                )
                time.sleep(self.throttle)
                continue

            # boss load updated decision masks and contigs
            reload_mapper, mapper = boss.reload(conf=self.conf)
            if reload_mapper:
                self.mapper = mapper  # replace the mapper with new one

            last_live_toml_mtime = self.reload_toml(last_live_toml_mtime)
            ########### Main Loop ###########
            loop_counter += 1
            number_reads = 0
            unblock_batch_action_list = []
            stop_receiving_action_list = []

            chunks = self.client.get_read_chunks(self.client.channel_count, last=True)
            calls = self.caller.basecall(
                chunks, self.client.signal_dtype, self.client.calibration_values
            )
            aligns = self.mapper.map_reads(calls)

            #######################################################################
            for result in aligns:
                number_reads += 1
                control, condition = self.conf.get_conditions(
                    result.channel, result.barcode
                )
                # result.decision = make_decision(self.conf, result)
                # boss decisions
                result.decision = boss.make_decision_boss(self.conf, result)
                action = condition.get_action(result.decision)
                seen_count = self.chunk_tracker.seen(result.channel, result.read_number)
                #  Check if there any conditions that override the action chose, exceed_max_chunks etc...
                (
                    previous_action,
                    action_overridden,
                    overridden_action_name,
                ) = self.check_override_action(
                    control,
                    action,
                    result,
                    seen_count,
                    condition,
                    stop_receiving_action_list,
                    unblock_batch_action_list,
                )
                self.loop_statistics.log_read(
                    client_iteration=loop_counter,
                    read_in_loop=number_reads,
                    read_id=result.read_id,
                    channel=result.channel,
                    read_number=result.read_number,
                    seq_len=len(result.seq),
                    counter=seen_count,
                    mode=result.decision.name,
                    decision=action.name,
                    condition=condition.name,
                    barcode=result.barcode,
                    previous_action=(
                        previous_action.name
                        if previous_action is not None
                        else previous_action
                    ),
                    action_overridden=action_overridden,
                    timestamp=time.time(),
                    # Anything below here is not included in the Debug log
                    region_name=(
                        _region.name
                        if (_region := self.conf.get_region(result.channel))
                        else "flowcell"
                    ),
                    overridden_action_name=overridden_action_name,
                )

            #######################################################################
            # Compile actions to be sent
            self.client.unblock_read_batch(
                unblock_batch_action_list, duration=self.unblock_duration
            )
            self.client.stop_receiving_batch(stop_receiving_action_list)

            t1 = timer()
            if number_reads > 0:
                self.loop_statistics.add_batch_performance(
                    number_of_reads=number_reads, batch_time=t1 - t0
                )
                self.logger.info(self.loop_statistics.get_batch_performance())

            # limit the rate at which we make requests
            if t0 + self.throttle > t1:
                time.sleep(self.throttle + t0 - t1)
        else:
            send_message(
                self.client.connection,
                "Readfish client stopped.",
                Severity.WARN,
            )
            self.caller.disconnect()
            self.mapper.disconnect()
            self.logger.info("Finished analysis of reads as client stopped.")






# TODO the only change in this function is the renamed AnalysisMod object
# TODO to inject the modified run() function above
def run(
    parser: argparse.ArgumentParser, args: argparse.Namespace, extras: list[Any]
) -> int:
    """Run function for targets.py

    Imported in `_cli_base.py`.
    Sets up the read until client and starts the analysis thread above.

    :param parser: Argparse onject - unused but must be taken due as may be needed
    :param args: The arguments passed to ArgParse
    :param extras: Extra stuff, I guess

    :returns: An exit code integer, 0 for success
    """
    # Setup logger used in this entry point, this one should be passed through
    logger = logging.getLogger(f"readfish.{args.command}")

    # Fetch sequencing device
    position = get_device(args.device, host=args.host, port=args.port)

    # Create a read until client
    read_until_client = RUClient(
        mk_host=position.host,
        mk_port=position.description.rpc_ports.secure,
        filter_strands=True,
        cache_type=AccumulatingCache,
        timeout=args.wait_for_ready,
    )

    # Load TOML configuration
    conf = Conf.from_file(args.toml, read_until_client.channel_count, logger=logger)

    # Set the padding if it is specified.
    if padding := getattr(args, "padding", None):
        for region in conf.regions:
            region.targets.padding = padding
        for barcode in conf.barcodes:
            conf.barcodes[barcode].targets.padding = padding
    logger.info(conf.describe_experiment())

    send_message(
        read_until_client.connection,
        f"'readfish {args.command}' connected to this device.",
        Severity.WARN,
    )

    # start the client running
    read_until_client.run(
        # TODO: Set correct channel range
        # first_channel=186,
        # last_channel=187,
        first_channel=1,
        last_channel=read_until_client.channel_count,
        max_unblock_read_length_seconds=args.max_unblock_read_length_seconds,
    )

    worker = AnalysisMod(
        read_until_client,
        conf=conf,
        logger=logger,
        debug_log=args.debug_log,
        unblock_duration=args.unblock_duration,
        throttle=args.throttle,
        dry_run=args.dry_run,
        toml=args.toml,
    )

    # begin readfish function
    try:
        worker.run()
    except KeyboardInterrupt:
        logger.info("Keyboard interrupt received, stopping readfish.")
        pass
    finally:
        read_until_client.reset()

    send_message(
        read_until_client.connection,
        "Readfish disconnected from this device. Sequencing will proceed normally.",
        Severity.WARN,
    )
    return 0






class BossBits:

    def __init__(self, conf, logger, mapper):
        # initialise the bits that are necessary for dynamic updates
        self.scale_factor = 100  # make param in future maybe
        self.masks = {}
        self.logger = logger
        self.mapper = mapper
        self.last_mask_mtime = 0
        self.last_contig_mtime = 0

        # Overwrite the default strand_converter
        # used to index into the strategies,
        # which are shaped 2xN for forward & reverse in 1st dim
        # readfish uses 1 for forward, -1 for reverse
        # For BOSS we use 0/False for forward, 1/True for reverse
        self.strand_converter = {1: False, -1: True}

        # grab path to masks and contigs, requires presence of 1 BOSS region in toml
        for region in conf.regions:
            # get the regions that is not the control
            # this currently limits using BOSS in a single region
            if region.control is False:
                self.mask_path = Path(f"out_{region.name}/masks")
                self.cont_path = Path(f"out_{region.name}/contigs")
                if not self.mask_path.is_dir():
                    os.makedirs(self.mask_path)
                    os.makedirs(self.cont_path)

        self.logger.info(f"Looking for BOSS related files at: \n"
                         f" strategies: {self.mask_path} \n"
                         f" contigs: {self.cont_path}")
        # if there is no strategy, generate dummy to start with
        # is_empty = not any(self.mask_path.iterdir())
        # if is_empty:
        self.logger.info("Creating dummy strategy")
        contig_strats = {'init': np.ones(1)}
        np.savez(self.mask_path / "boss", **contig_strats)


    @staticmethod
    def _reload_npy(mask_files):
        return {path.stem: np.load(path) for path in mask_files}


    @staticmethod
    def _reload_npz(mask_files):
        mask_container = np.load(mask_files[0])
        return {name: mask_container[name] for name in mask_container}


    @staticmethod
    def _idx(fa):
        idx_name = fa.stem + ".mmi"
        mappy.Aligner(fn_idx_in=str(fa), fn_idx_out=str(idx_name), preset="map-ont")
        return idx_name


    def _reload_masks(self):
        """
        Reload updated decision masks.
        Loads all present .npy/.npz files
        """
        # can use multiple .npy or a single .npz file
        # BOSS* uses a single .npz in all modes now
        new_masks = list(self.mask_path.glob("*.npy"))
        if new_masks:
            reload_func = self._reload_npy
        elif not new_masks:
            new_masks = list(self.mask_path.glob("*.npz"))
            reload_func = self._reload_npz
        else:
            raise FileNotFoundError("No mask files present")

        # Do we actually update this time?
        if not new_masks[0].stat().st_mtime > self.last_mask_mtime:
            return 0

        try:
            mask_dict = reload_func(mask_files=new_masks)
            self.masks = mask_dict
            self.logger.info(f"Reloaded strategies for {len(set(self.masks.keys()))} sequences")
        except Exception as e:
            self.logger.error(f"Error reading strategy array ->>> {repr(e)}")
            self.masks = {"exception": True}
        # update last mtime
        self.last_mask_mtime = new_masks[0].stat().st_mtime
        return 1


    def _reload_mapper(self, conf):
        """
        Reload mapping index for decision-making. Only load if marker file exists.
        """
        # Do we actually update this time?
        # check if contigs are newer than before
        cs = self.cont_path / "aeons.fa"

        if not cs.stat().st_mtime > self.last_contig_mtime:
            return

        try:
            self.logger.info(f"Regenerating mapping index")
            idx_name = self._idx(fa=cs)
            self.logger.info(f"Reloading mapping index")
            conf.mapper_settings.parameters['fn_idx_in'] = idx_name
            mapper: AlignerABC = conf.mapper_settings.load_object("Aligner")
            self.mapper = mapper
            self.logger.info("Aligner initialised")
            while not self.mapper.initialised:
                time.sleep(1)
            # new last mtime
            self.last_contig_mtime = cs.stat().st_mtime

        except Exception as e:
            self.logger.error(f"Error loading mapping index ->>> {repr(e)}")


    def _check_names(self):
        """
        Check that the names of masks and sequences in the mapper are the same.
        """
        if not self.mapper.initialised:
            return

        mask_names = sorted(list(self.masks.keys()))
        contig_names = sorted(list(self.mapper.aligner.seq_names))
        same_names = mask_names == contig_names
        if not same_names:
            self.logger.error(f"Error loading masks and contigs: discrepancy in names {len(mask_names)} {len(contig_names)} \n\n\n\n\n")


    def reload(self, conf):
        # masks are reloaded from mask files
        reloaded_masks = self._reload_masks()
        if not reloaded_masks:
            return 0, self.mapper
        # return if running BOSS-RUNS
        if not (self.cont_path / "aeons.fa").exists():
            return 0, self.mapper
        self._reload_mapper(conf)
        # check that the names of the masks and contigs are the same
        self._check_names()
        return 1, self.mapper


    def _check_coord(self, contig, start_pos, reverse):
        """
        Query numpy mask array and return decision to keep sequencing
        Parameters
        ----------
        start_pos: int
            Mapping start coordinate
        reverse: bool
            Whether the coordinate is on the forward or reverse strand
        contig: str
            The name of the contig the mask is for, used to lookup correct mask in the directory
        Returns
        -------
        bool
            Decision to keep sequencing
        """
        # this happens when there was an error reading the mask file, then accept everything
        if self.masks.get("exception", False):
            return 1
        if contig not in self.masks:
            self.logger.warning(f"{contig} is not in mask dict")
            return 1
        arr = self.masks.get(contig)
        try:
            d = arr[:, int(reverse)][start_pos // self.scale_factor]
            return d
        except Exception as e:
            return 1


    def make_decision_boss(self, conf: Conf, result: Result) -> Decision:
        # TODO this is modified from readfishes _config.py
        """
        The main decision making function for readfish.
        Chooses the decision that is looked up in the TOML based
        on the mapped coordinates of the read, checked against the targets.
        Decision is one of single_on, multi_on, single_off, multi_off, no_map, no_seq.

        :param conf: The Conf object for the experiment
        :param result: The result of the alignment and base calling for the read
        :raises ValueError: If readfish fails to make a decision based on
            the passed Result object length
        :return: The decision that readfish has made for this read
        """
        if result.alignment_data is None:
            result.alignment_data = []
        # targets = conf.get_targets(result.channel, result.barcode)
        results = result.alignment_data
        matches = []
        for al in results:
            contig = al.ctg
            strand = al.strand
            coord = al.r_st if al.strand == 1 else al.r_en  # TODO this is different from readfish now
            # matches.append(targets.check_coord(contig, strand, coord))
            strand_conv = self.strand_converter[al.strand]
            matches.append(self._check_coord(contig=contig, start_pos=coord, reverse=strand_conv))
        coord_match = any(matches)

        if not results:
            if len(result.seq) > 0:
                return Decision.no_map
            else:
                return Decision.no_seq
        elif len(results) == 1:
            return Decision.single_on if coord_match else Decision.single_off
        elif len(results) > 1:
            return Decision.multi_on if coord_match else Decision.multi_off
        raise ValueError()


    @staticmethod
    def gen_dummy_idx() -> None:
        """
        Generate a dummy index for readfish to start out with

        :return:
        """
        # generate dummy index
        with open("dummy.fa", "w") as dfa:
            dfa.write(">init\nAAAAAAAAAAAAAAAAAAAAAAAAA")
        mappy.Aligner(fn_idx_in="dummy.fa", fn_idx_out="readfish.mmi", preset="map-ont")
        Path("dummy.fa").unlink()




def get_args(arg_list: list = None) -> tuple[argparse.ArgumentParser, argparse.Namespace]:
    """
    Get the default arguments and extend with toml file using _cli_base.main()
    :param arg_list: List to pass in arguments for debugging
    :return: parser and arguments to pass to the entry-point
    """
    if not arg_list:
        toml_readfish = sys.argv[1]
        device = sys.argv[2]
        name = sys.argv[3]
        # host = sys.argv[3]
        # port = sys.argv[4]
    else:
        toml_readfish = arg_list[0]
        device = arg_list[1]
        name = arg_list[2]

    # check arguments to add in future with: readfish targets -h
    argv = [
        'targets',
        '--toml', toml_readfish,
        '--device', device,
        '--experiment-name', name,
        # '--log-file', 'readfish.log'
        # '--host', host,
        # '--port', port,
    ]
    if arg_list:
        argv.append("--log-file")
        argv.append("readfish.log")
    parser, args = main_args(argv=argv)
    return parser, args



def main(arg_list: list = None) -> None:
    """
    This is equivalent to launching entry-point in readfish from _cli_base.py
    in readfish main() of _cli_base.py would launch the entry-point,
    whereas here we grab the arguments and then launch a modified run function
    :param arg_list: List to pass in arguments for debugging
    :return:
    """
    BossBits.gen_dummy_idx()
    # get arguments
    if not arg_list:
        arg_list = sys.argv[1:]
    parser, args = get_args(arg_list)
    # launch entry-point
    run(parser=parser, args=args, extras=[])


if __name__ == "__main__":
    main()
