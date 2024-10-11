import argparse
from types import SimpleNamespace
from pathlib import Path
from datetime import datetime
import sys

import rtoml

from boss.utils import init_logger


"""
Configuration:
1. load a template toml that contains defaults
2. parse a toml given on command line
3. overwrite defaults with CL toml
4. parse a toml for readfish separately
5. exchange args between the two
"""



class Config:
    def __init__(self, parse: bool = False, arg_list: list = None):
        """
        Initialise configuration by loading defaults
        When debugging, paths to the TOMLs can be given as a list or arguments to parse
        For simulations, providing a readfish TOML is not necessary

        :param parse: Parse from command line/toml_paths
        :param arg_list: Pass-through list of paths to toml files
        """
        self.template_toml = """
        [general]
        name = "boss"                   # experiment name
        wait = 60                       # waiting time between updates in live version
        ref = ""                        # reference file (fasta or None). Not specifying a file switches to AEONS
        mmi = ""                        # index of reference (will be built if ref is given but not mmi)
        
        [live]
        device = "TEST"                 # position on sequencer
        host = "localhost"              # host of sequencing device
        port = 9502                     # port of sequencing device
        data_wait = 100                 # wait for X Mb of data before first update

        [optional]
        reject_refs = ""                # comma-separated list of headers in reference from which to always reject
        ploidy = 1                      # 1 or 2
        lowcov = 10                     # target coverage for assemblies
        temperature = 60                # max updates during which to consider fragments (higher number might decrease update speeds)
        min_seq_len = 2500              # min sequence length used during contig reconstruction
        min_contig_len = 10_000         # min length to consider contigs
        min_s1 = 200                    # min alignment scores to consider overlaps
        min_map_len = 2000              # min alignment length to consider overlaps
        tetra = true                    # perform tetranucleotide frequency tests
        filter_repeats = false          # perform repeat filtering
        bucket_threshold = 5            # at which coverage to switch on the strategy in a bucket (debug)

        ########################################################################################

        [simulation]                    # simulation arguments
        fq = ""
        batchsize = 4000                    
        maxb = 400
        binit = 5
        dumptime = 200000000
        paf_full = ""                   # giving pafs triggers BOSS-RUNS
        paf_trunc = ""
        """
        # load the default from above
        self.args = rtoml.loads(self.template_toml)
        # convert toml dict to SimpleNamespace
        self.args = self._convert_to_namespace()
        # do we parse toml paths to overwrite defaults?
        if parse:
            self.arg_list = arg_list
            # load TOMLs from command line if not given as list
            if not self.arg_list:
                self.arg_list = sys.argv[1:]

            toml_paths = self._parse_toml_args()
            args_file, args_readfish = self._load_tomls(toml_paths)
            self._overwrite_defaults(args_file)
            # check if we are simulating or running real experiment
            self._check_run_type()
            # initialise a log file in the output folder
            stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            self.logfile = f'{stamp}_boss.log'
            init_logger(logfile=self.logfile, args=self.args)
            # config settings for readfish
            if self.args.live_run:
                # add path to readfish toml as arg
                self.args.toml_readfish = toml_paths.toml_readfish
                # check that experiment exists as readfish region
                self._verify_region_names(self.args, args_readfish)
                # validate readfish args
                self._validate_readfish_conf(args_readfish)



    def _check_run_type(self) -> None:
        """
        Check if we are running a simulation or a real experiment
        If "fq" under the simulation header is given, we simulate

        :return:
        """
        if self.args.fq:
            self.args.sim_run = True
            self.args.live_run = False
        elif self.args.device:
            self.args.live_run = True
            self.args.sim_run = False
        else:
            raise ValueError("Need either fastq for simulation or device for live run")



    def _convert_to_namespace(self) -> SimpleNamespace:
        """
        Convert arguments from a parsed Dict to a Namespace object
        For method-style access to attributes

        :return: Arguments as SimpleNamespace object
        """
        args = SimpleNamespace()
        for category, subdict in self.args.items():
            if not type(subdict) is dict:
                setattr(args, category, subdict)
            else:
                for k, v in subdict.items():
                    setattr(args, k, v)
        return args



    def _parse_toml_args(self) -> argparse.Namespace:
        """
        Parse TOML paths given as arguments on the command line

        :return: argparse.Namespace with paths to the two TOMLs
        """
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--toml',
            type=str,
            default=None,
            required=True,
            help='TOML configuration file')
        parser.add_argument(
            '--toml_readfish',
            type=str,
            default=None,
            help='TOML configuration file for readfish')
        toml_paths = parser.parse_args(self.arg_list)
        return toml_paths



    @staticmethod
    def _load_tomls(toml_paths: argparse.Namespace) -> tuple[dict, dict]:
        """
        Load the TOML files into dictionaries using rTOML

        :param toml_paths: Paths to TOMLs as parsed arguments
        :return: Tuple with parsed TOML dictionaries
        """
        args_file = rtoml.loads(Path(toml_paths.toml).read_text(encoding="utf-8"))
        if toml_paths.toml_readfish:
            args_readfish_file = rtoml.loads(Path(toml_paths.toml_readfish).read_text(encoding="utf-8"))
        else:
            args_readfish_file = dict()
        return args_file, args_readfish_file



    def _overwrite_defaults(self, args_from_file: dict) -> None:
        """
        Use TOML given on CL to overwrite defaults

        :param args_from_file: parsed contents of given TOML
        :return:
        """
        for category in args_from_file:
            for k, v in args_from_file[category].items():
                setattr(self.args, k, v)



    @staticmethod
    def _verify_region_names(args: SimpleNamespace, args_readfish: dict) -> None:
        """
        Verify that the experiment name of BOSS exists as region in readfish

        :param args: Config dictionary for BOSS
        :param args_readfish: Config dictionary for readfish
        :return:
        """
        if type(args_readfish['regions']) is not list:
            raise ValueError("Readfish regions must be specified as array")

        # make sure the names of BOSS and regions on flowcell are the same
        region_names = {r['name'] for r in args_readfish['regions']}
        if args.name not in region_names:
            raise ValueError("One of the regions in readfish needs the same name as the experiment in BOSS")



    @staticmethod
    def _validate_readfish_conf(args_rf: dict) -> int:
        """
        Minimalist version of readfish validate entry-point
        to validate parsed toml config for readfish

        :param args_rf: Dict of arguments parsed for readfish
        :return: indicator
        """
        from readfish._config import Conf
        channels = 512
        try:
            _ = Conf.from_dict(args_rf, channels)
        except:
            raise ValueError("Could not load TOML config for readfish")
        return 0





