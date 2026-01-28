import argparse
from pathlib import Path
from datetime import datetime
import sys
import tomllib

from pydantic import BaseModel, ValidationError, Field

from boss.utils import init_logger


"""
Configuration:
1. load a template toml that contains defaults
2. parse a toml given on command line
3. overwrite defaults with CL toml
4. parse a toml for readfish separately
5. exchange args between the two
"""




class GeneralConfig(BaseModel):
    name : str = Field(default='boss', description="Experiment name. Used as output prefix and to match readfish region name")
    ref: str | None = Field(default=None, description="Reference file (fasta or None). Not specifying a file switches operation to AEONS")
    mmi: str | None = Field(default=None, description="Index of reference (will be built if not provided)")
    toml_readfish: str | None = Field(default=None, description="TOML config file for readfish. Not required for simulations.")
    wait: int = Field(default=60, description="Waiting time between updates in live version")


class LiveConfig(BaseModel):
    device: str | None = Field(default=None, description="Position on sequencing device")
    host: str = Field(default='localhost', description="Host of sequencing device")
    port: int = Field(default=9502, description="Port of sequencing device")
    data_wait: int = Field(default=100, description="Wait for X Mb of data before first strategy update")


class OptionalConfig(BaseModel):
    reject_refs: str | None = Field(default=None, description="Comma-separated list of headers in reference from which to always reject")
    ploidy: int = Field(default=1, description="Ploidy level")
    lowcov: int = Field(default=10, description="[debug] Minimum coverage")
    temperature: int = Field(default=60, description="[debug] Temperature")
    min_seq_len: int = Field(default=2500, description="[debug] Minimum sequence length")
    min_contig_len: int = Field(default=10_000, description="[debug] Minimum contig length")
    min_s1: int = Field(default=200, description="[debug] Minimum S1")
    min_map_len: int = Field(default=2000, description="[debug] Minimum mapping length")
    tetra: bool = Field(default=True, description="[debug] Switch tetranucleotide frequency tests")
    filter_repeats: bool = Field(default=False, description="[debug] Switch repeat filtering")
    bucket_threshold: int = Field(default=5, description="[debug] At which coverage to switch on the strategy in a bucket")


class SimulationConfig(BaseModel):
    fq: str | None = Field(default=None, description="Input fastq file")
    batchsize: int = Field(default=4000, description="Number of reads per update")
    maxb: int = Field(default=400, description="Maximum number of batches")
    binit: int = Field(default=5, description="Initial batch size")
    dumptime: int = Field(default=200000000, description="Time (in units of psudo-sequencing time) between writing output fastq files")
    paf_full: str | None = Field(default=None, description="Mappings (PAF) of full-length reads for fast sampling")
    paf_trunc: str | None = Field(default=None, description="Mappings (PAF) of truncated reads for fast sampling")
    accept_unmapped: bool = Field(default=False, description="Accept unmapped reads")


class BossConfig(BaseModel):
    general : GeneralConfig = GeneralConfig()          
    live : LiveConfig = LiveConfig()                   
    optional : OptionalConfig = OptionalConfig()       
    simulation : SimulationConfig = SimulationConfig() 





class Config:
    def __init__(self, parse: bool = False):
        """
        Initialise configuration by loading defaults
        When debugging, paths to the TOMLs can be given as a list or arguments to parse
        For simulations, providing a readfish TOML is not necessary

        :param parse: Whether to parse command line arguments for TOML config
        """
        # init defaults
        self.args = BossConfig()

        if parse:
            toml_path = self._parse_toml_arg()

            # load and validate main config
            try:
                with Path(toml_path).open("rb") as f:
                    conf = tomllib.load(f)
                self.args = BossConfig.model_validate(conf)
            except ValidationError as e:
                print("Invalid configuration:")
                print(e)
                sys.exit(1)

        # load readfish config if given, just to validate it here
        if self.args.general.toml_readfish:
            args_readfish = tomllib.loads(Path(self.args.general.toml_readfish).read_text(encoding="utf-8"))
        else:
            args_readfish = dict()

        # initialise a log file in the output folder
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        Path("./logs").mkdir(parents=True, exist_ok=True)
        self.logfile = f'./logs/{stamp}_boss.log'
        init_logger(logfile=self.logfile, args=self.args)
        # config settings for readfish
        if self.args.live.device:
            # check that experiment exists as readfish region
            self._verify_region_names(self.args, args_readfish)
            # validate readfish args
            self._validate_readfish_conf(args_readfish)



    @staticmethod
    def write_template(path: Path = Path("config_template.toml")) -> None:
        """
        Write a template config file to the specified path.

        :param path: Path to the template config file
        """        
        VALUE_COL = 30
        config = ""
        for section_name, section in BossConfig.model_fields.items():
            config += f"\n[{section_name}]"
            section_type = section.annotation
            for key, field in section_type.model_fields.items():   # type: ignore
                key_val = f"{key} = {field.default!r}"
                config += f"\n{key_val:<{VALUE_COL}}  # {field.description}"
            config += "\n"

        print(config)
        with path.open("w") as f:
            f.write(config)




    def _parse_toml_arg(self) -> str:
        """
        Parse TOML path given as argument on the command line

        :return: path to TOML config file
        """
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--toml',
            type=str,
            default=None,
            required=True,
            help='TOML configuration file')
        toml_path = parser.parse_args()
        return toml_path.toml




    @staticmethod
    def _verify_region_names(args, args_readfish: dict) -> None:
        """
        Verify that the experiment name of BOSS exists as region in readfish
        TODO this is not relevant if barcodes are used, then there might not be regions

        :param args: Config dictionary for BOSS
        :param args_readfish: Config dictionary for readfish
        :return:
        """
        if type(args_readfish['regions']) is not list:
            raise ValueError("Readfish regions must be specified as array")

        # make sure the names of BOSS and regions on flowcell are the same
        region_names = {r['name'] for r in args_readfish['regions']}
        if args.general.name not in region_names:
            raise ValueError("One of the regions in readfish needs the same name as the experiment in BOSS")

        # TODO: Add check that all boss-runs barcodes are also in readfish



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
        except:  # noqa: E722
            raise ValueError("Could not load TOML config for readfish")
        return 0


if __name__ == "__main__":
    Config.write_template()

