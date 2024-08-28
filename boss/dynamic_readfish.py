import sys
import os
from pathlib import Path
import time
import argparse

import numpy as np
import mappy

from boss._cli_base import main as main_args

from readfish._config import Conf  # noqa
from readfish.plugins.abc import AlignerABC
from readfish.plugins.utils import Decision, Result



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
        if not Path(self.mask_path / 'boss.npz').is_file():
            self.logger.info("Creating dummy strategy")
            contig_strats = {'init': np.ones(1)}
            np.savez(self.mask_path / "boss", **contig_strats)
        else:
            self.logger.info("Loading starting strategy")


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
        """
        # BOSS* uses a single .npz in all modes now
        new_masks = list(self.mask_path.glob("boss.npz"))
        reload_func = self._reload_npz
        if not new_masks:
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
        # reject refs have array size 1
        if arr.shape[0] == 1:
            return 0
        # otherwise query the strategy array
        try:
            d = arr[:, int(reverse)][start_pos // self.scale_factor]
            return d
        except Exception as e:  # noqa
            return 1


    def make_decision_boss(self, conf: Conf, result: Result) -> Decision:  # noqa
        """
        The main decision-making function for readfish.
        Modified from _config.py in readfish
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
            strand = al.strand  # noqa
            coord = al.r_st if al.strand == 1 else al.r_en
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
        argv.append("--debug-log")
        argv.append("readfish_chunks.tsv")
    parser, args = main_args(argv=argv)
    return parser, args


