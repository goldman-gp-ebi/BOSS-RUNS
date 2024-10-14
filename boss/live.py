import logging
import subprocess
import os
import time
import glob
import inspect
from pathlib import Path
from datetime import datetime

import rtoml
from minknow_api.manager import Manager, FlowCellPosition
from minknow_api import __version__ as minknow_api_version




class Sequencer:

    def __init__(self, position: FlowCellPosition = None):
        """
        class that represents a connected sequencing device

        :param position: FlowCellPosition instance from minknow_api
        """
        self.position = position

        if self.position:
            self._grab_output_dir()
            self._grab_device_type()
        else:   # testing pass-through
            self.channels = set()
            self.out_path = '../data'
            self.device_type = 'min'
            if not os.path.exists(f'{self.out_path}/fastq_pass'):
                os.mkdir(f'{self.out_path}/fastq_pass')




    def _grab_output_dir(self) -> None:
        """
        Capture the output directory of MinKNOW,
        i.e. where fastq files are deposited during sequencing
        :return:
        """
        # connect to the device and navigate api to get output path
        device_connection = self.position.connect()
        current_run = device_connection.protocol.get_current_protocol_run()
        run_id = current_run.run_id
        logging.info(f"connected to run_id: {run_id}")
        self.out_path = str(current_run.output_path)
        logging.info(f"grabbing Minknow's output path: \n{self.out_path}\n")



    def _grab_device_type(self) -> None:
        """
        Get info on the type of device that the sequencing is running on.
        Can be either minion, gridion, promethion, p2_solo
        Currently not used for anything

        :return:
        """
        pro_device_types = {'P2_SOLO', 'PROMETHION'}
        min_device_types = {'MINION'}
        dt = self.position.device_type
        if dt in pro_device_types:
            self.device_type = 'pro'
        elif dt in min_device_types:
            self.device_type = 'min'
        else:
            logging.info(f'WARNING: device type {dt} not recognized. Using MINION flowcell layout as fallback')
            self.device_type = 'min'



    def grab_channels(self, run_name: str) -> None:
        """
        Wait until we find a channels TOML spec that contains the channel
        numbers assigned to the BOSS region on the flowcell.
        Result is the set of channel IDs from which data is used or empty set

        :param run_name: Experiment name from config TOML
        :return:
        """
        self.channels_toml = f'{self.out_path}/channels.toml'
        logging.info(f'looking for channels specification at : {self.channels_toml}')
        channels_found = False
        channels = []
        while not channels_found:
            if not os.path.isfile(self.channels_toml):
                logging.info("channels file does not exist (yet), waiting for 30s")
                time.sleep(30)
            else:
                channels = self._parse_channels_toml(run_name=run_name)
                channels_found = True
        # channels successfully found
        self.channels = channels



    def _parse_channels_toml(self, run_name: str) -> set:
        """
        Look into the channels toml that readfish writes
        This toml contains lists of channels assigned to each region
        Grab the channel numbers for the BOSS region

        :param run_name: experiment name of BOSS region
        :return: Set of channel numbers from which to consider data
        """
        toml_dict = rtoml.load(Path(self.channels_toml))
        # if there is only one condition, we return empty set
        # that way we can skip the regex when scanning new data
        if len(toml_dict["conditions"]) == 1:
            logging.info('Only one condition. Using all channels!')
            return set()

        # otherwise find the corresponding condition
        correct_key = ''
        for key in toml_dict["conditions"].keys():
            name = toml_dict["conditions"][key]["name"]
            if name == run_name:
                correct_key = key
                break

        if not correct_key:
            raise ValueError(f"Experiment name {run_name} in .toml not found in channel-specification toml.")

        selected_channels = set(toml_dict["conditions"][correct_key]["channels"])
        logging.info(f"found channels specification: Using {len(selected_channels)} channels.")
        return selected_channels




class LiveRun:


    @staticmethod
    def connect_sequencer(device: str, host: str = 'localhost', port: int = 9502) -> Sequencer:
        """
        Connect to the running sequencer to get the path to its output directory

        :param device: Device name on sequencing machine
        :param host: Host to connect to
        :param port: Possibility to overwrite default port
        :return: Sequencer object that holds info about the flowcellposition
        """
        LiveRun._check_minknow_api_version()
        # try connecting to the sequencing device
        try:
            flowcellposition = LiveRun._grab_target_device(device=device, host=host, port=port)
            sequencer = Sequencer(position=flowcellposition)
        except:
            raise ValueError(f"Error: issue connecting to target device {device}, at host {host} and port {port}. \n"
                             f"Please make sure to supply correct name of sequencing position in MinKNOW")
        return sequencer


    @staticmethod
    def _check_minknow_api_version() -> None:
        """
        Check compatibility of the available minknow_api version. Hard exits if not compatible
        :return:
        """
        logging.info(f"minknow API Version {minknow_api_version}")
        # minknow_api.manager supplies Manager (wrapper around MinKNOW's Manager gRPC)
        if not minknow_api_version.startswith("6"):
            raise NotImplementedError("Unsupported version of minknow_api. MinKnow <6 is not supported.")



    @staticmethod
    def _grab_target_device(device: str, host: str = 'localhost', port: int = 9502) -> FlowCellPosition:
        """
        Get the flowcell position that runs the sequencing process

        :param device: device name of the 'position' in the sequencing machine
        :param host: hostname to connect to for MinKNOW
        :param port: override default port to connect to
        :return: instance of FlowCellPosition from minknow_api
        """
        # Find a list of currently available sequencing positions.
        manager = Manager(host=host, port=int(port))
        positions = list(manager.flow_cell_positions())
        pos_dict = {pos.name: pos for pos in positions}
        # index into the dict of available devices
        target_device = pos_dict[device]
        return target_device



    @staticmethod
    def scan_dir(fastq_pass: str, processed_files: set) -> list[str]:
        """
        Periodically scan Minknow's output dir and
        return list of all NEW files

        :param fastq_pass: Path of MinKNOW's output
        :param processed_files: Set of already processed files
        :return: List of new, previously unprocessed files
        """
        patterns = ["*.fq.gz", "*.fastq.gz", "*.fastq.gzip", "*.fq.gzip", "*.fastq", "*.fq"]
        all_fq = set()
        for p in patterns:
            all_fq.update(glob.glob(f'{fastq_pass}/{p}'))
        # which files have we not seen before?
        new_fq = all_fq.difference(processed_files)
        logging.info(f"found {len(new_fq)} new fq files")
        new_fq_list = [f for f in new_fq]
        return new_fq_list



    @staticmethod
    def launch_readfish(toml: str, device: str, name: str) -> str:
        """
        Wrapper to launch readfish into the background

        :param toml: TOML config for readfish
        :param device: Flowcell position
        :param name: Name of BOSS experiment
        :return: name of the readfish log file
        """
        if device == "TEST":  # exit for testing purposes
            return ''
        if not Path(toml).exists():
            raise FileNotFoundError("Specified readfish toml does not exist.")
        # find the script to launch readfish
        module_path = inspect.getfile(LiveRun)
        logging.info(module_path)
        script_path = Path(module_path).parent / "readfish_boss.py"
        if not script_path.is_file():
            raise FileNotFoundError("readfish_boss.py not found. Something went wrong..")
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        readfish_log = f'{stamp}_readfish.log'
        readfish_comm = f'python {script_path} {toml} {device} {name} >{readfish_log} 2>&1'
        logging.info(readfish_comm)
        # launch readfish into the background
        logging.info("Launching readfish")
        subprocess.Popen(readfish_comm, shell=True)
        return readfish_log


