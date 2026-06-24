import time
import logging
from pathlib import Path
from typing import Callable

from boss.live import LiveRun, Sequencer
from boss.readlengthdist import ReadlengthDist
from boss.batch import FastqBatch
from boss.config import BossConfig



class Boss:

    def __init__(self, args: BossConfig):
        """
        Main initialisation of experiment

        :param args: Configuration of the experiment
        """
        self.args = args
        self.name = args.general.name
        self.processed_files = set()
        self.n_fastq = 0
        self.batch = 0
        # initial strategy is to accept
        self.strat = 1
        # initialise output directory structure
        self._init_file_struct()
        # object to record read lengths
        self.rl_dist = ReadlengthDist()  # NOTE: readlength distribution will not differ by barcodes for now



    def _init_file_struct(self) -> None:
        """
        Set up the required directory structure

        :return:
        """
        # make sure the run name does not have any spaces
        assert ' ' not in self.name

        self.out_dir = f'./out_{self.name}'
        out_path = Path(self.out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        (out_path / "masks").mkdir(parents=True, exist_ok=True)
        # NOTE: Could process barcodes into separate outputfiles, will not affect the directory structure, but something to keep in mind
        (out_path / "fq").mkdir(parents=True, exist_ok=True)
        (out_path / "logs").mkdir(parents=True, exist_ok=True)
        (out_path / "contigs").mkdir(parents=True, exist_ok=True)
        (out_path / "contigs" / "prev").mkdir(parents=True, exist_ok=True)
        (out_path / "contigs" / "init").mkdir(parents=True, exist_ok=True)
        (out_path / "metrics").mkdir(parents=True, exist_ok=True)
        (out_path / "tmp").mkdir(parents=True, exist_ok=True)


    def _init_live(self) -> None:
        """
        Wrapper to initialise a live experiment, incl. device related things
        find output path where the fastq files are placed
        find channels toml for runs with multiple regions

        :return:
        """
        # connect to sequencing machine
        # early exit for pytests
        if not self.args.live.device:
            sequencer = Sequencer()

        else:
            sequencer = LiveRun.connect_sequencer(device=self.args.live.device,
                                                  host=self.args.live.host,
                                                  port=self.args.live.port)
            sequencer.grab_channels(run_name=self.name)

        # get the relevant infos from the Sequencer
        self.fq = f'{sequencer.out_path}/fastq_pass' 
        # TODO: Find out if Minknow always places barcoded samples in subfolders 
        # or whether there are config options. Lukas believes we should be able to handle default settings
        assert Path(self.fq).is_dir()
        # readfish needs to have placed the channels.toml at this point
        # channels can be an empty set if there is only one condition
        # then data from all channels will be used (i.e. later regex skipped)
        self.channels = sequencer.channels



    def _get_new_data(self) -> tuple[dict, dict]:
        """
        Get new read data and parse it into a batch object

        :return: Dictionaries for read sequences and their qualities
        """
        # find new fastq files
        new_fastq = LiveRun.scan_dir(
            fastq_pass=self.fq, processed_files=self.processed_files)
        if not new_fastq:
            logging.info("no new files, deferring update ")
            return {}, {}
        # add the new files to the set of processed files
        self.processed_files.update(new_fastq)
        self.n_fastq += len(new_fastq)
        fq_batch = FastqBatch(fq_files=new_fastq, channels=self.channels)
        # read length dist is updated straight away, not in func of subclass
        self.rl_dist.update(read_lengths=fq_batch.read_lengths)
        return fq_batch.read_sequences, fq_batch.read_qualities


    def _time_to_next_update(self, tic: float) -> int:
        """
        Calculate how long to wait until the next update given duration of previous

        :param tic: Start time of previous update
        :return: Seconds to wait until next update
        """
        toc = time.time()
        passed = toc - tic
        next_update = int(self.args.general.wait - passed)
        logging.info(f"batch took: {passed}")
        logging.info(f"finished update, waiting for {next_update}s ... \n")
        return next_update


    def launch_live_components(self):
        # launch readfish and initialise connection to sequencer if live experiment
        assert self.args.general.toml_readfish is not None
        if self.args.live.device:
            LiveRun.launch_readfish(
                toml=self.args.general.toml_readfish,
                device=self.args.live.device,
                name=self.name
            )
        self._init_live()


    def process_batch(self, main_processing_func: Callable) -> int:
        """
        Main processing function for all live BOSS experiments.
        Variants can pass processing functions in

        :param main_processing_func:
        :return: Seconds to wait until next update
        """
        logging.info("")
        logging.info(f"Next batch ---------------------------- # {self.batch}")
        tic = time.time()

        new_reads, new_quals = self._get_new_data()
        if not new_reads:
            return self.args.general.wait

        main_processing_func(new_reads=new_reads, new_quals=new_quals)

        next_update = self._time_to_next_update(tic)
        self.batch += 1
        return next_update


    def process_batch_sim(self, main_processing_func: Callable) -> int:
        """
        Main entry point for simulated experiments. Variants can
        pass processing functions in. A sampler is also required.

        :param main_processing_func:
        :return: Seconds to wait until next update
        """
        logging.info("")
        logging.info(f"Next batch ---------------------------- # {self.batch}")
        tic = time.time()

        main_processing_func()

        next_update = self._time_to_next_update(tic)
        self.batch += 1
        return next_update








