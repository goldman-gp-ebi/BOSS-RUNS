import logging
import subprocess
import os
import time
import sys
import glob
import inspect
from pathlib import Path

import rtoml
from minknow_api.manager import Manager
from minknow_api import __version__ as minknow_api_version



class LiveRun:

    @staticmethod
    def split_flowcell(out_path: str, run_name: str) -> set:
        """
        Perform the necessary steps if running multiple regions on a flowcell
        Wait until we find a channels TOML spec that contains the channel
        numbers assigned to the BOSS region on the flowcell.

        :param out_path: General MinKNOW output path
        :param run_name: Experiment name from config TOML
        :return: Set of channel IDs to consider data from
        """
        channel_path = f'{out_path}/channels.toml'
        logging.info(f'looking for channels specification at : {channel_path}')
        channels_found = False
        channels = []
        while not channels_found:
            if not os.path.isfile(channel_path):
                logging.info("channels file does not exist (yet), waiting for 30s")
                time.sleep(30)
            else:
                channels = LiveRun._grab_channels(channels_toml=channel_path, run_name=run_name)
                channels_found = True
        # channels successfully found
        logging.info(f"found channels specification: Using {len(channels)} channels.")
        return channels



    @staticmethod
    def connect_sequencer(device: str, host: str = 'localhost', port: int = None) -> str:
        """
        Connect to the running sequencer to get the path to its output directory

        :param device: Device name on sequencing machine
        :param host: Host to connect to
        :param port: Possibility to overwrite default port
        :return: Path of MinKNOW output directory
        """
        try:
            out_path = LiveRun._grab_output_dir(device=device, host=host, port=port)
            logging.info(f"grabbing Minknow's output path: \n{out_path}\n")
        except:
            logging.info("Minknow's output dir could not be inferred from device name. Exiting.")
            logging.info(f'\n{device}\n{host}\n{port}')
            if device == "TEST":
                out_path = "../data"
                if not os.path.exists(f'{out_path}/fastq_pass'):
                    os.mkdir(f'{out_path}/fastq_pass')
            else:
                raise ValueError("Device not found.")
        return out_path



    @staticmethod
    def _grab_output_dir(device: str, host: str = 'localhost', port: int = None) -> str:
        """
        Capture the output directory of MinKNOW,
        i.e. where fastq files are deposited during sequencing
        host and port should be optional if run on the sequencing machine

        :param device: device name of the 'position' in the sequencing machine
        :param host: hostname to connect to for MinKNOW
        :param port: override default port to connect to
        :return: String of path where sequencing data is put my MinKNOW
        """
        logging.info(f"minknow API Version {minknow_api_version}")
        # minknow_api.manager supplies Manager (wrapper around MinKNOW's Manager gRPC)
        if minknow_api_version.startswith("5"):
            if not port:
                port = 9502
            manager = Manager(host=host, port=int(port))
        else:
            raise NotImplementedError("Unsupported version of minknow_api")

        # Find a list of currently available sequencing positions.
        positions = list(manager.flow_cell_positions())
        pos_dict = {pos.name: pos for pos in positions}
        # index into the dict of available devices
        try:
            target_device = pos_dict[device]
        except KeyError:
            logging.info(f"Error: target device {device} not available")
            logging.info("Error: Please make sure to supply correct name of sequencing position in MinKNOW.")
            sys.exit()
        # connect to the device and navigate api to get output path
        device_connection = target_device.connect()
        current_run = device_connection.protocol.get_current_protocol_run()
        run_id = current_run.run_id
        logging.info(f"connected to run_id: {run_id}")
        out_path = current_run.output_path
        return out_path


    @staticmethod
    def _grab_channels(channels_toml: str, run_name: str) -> set:
        """
        Look into the channels toml that readfish writes
        This toml contains lists of channels assigned to each region
        Grab the channel numbers for the BOSS region

        :param channels_toml: Path to channels TOML from readfish
        :param run_name: experiment name of BOSS region
        :return: Set of channel numbers from which to consider data
        """
        toml_dict = rtoml.load(Path(channels_toml))
        # find the corresponding condition
        correct_key = ''
        for key in toml_dict["conditions"].keys():
            name = toml_dict["conditions"][key]["name"]
            if name == run_name:
                correct_key = key
                break

        if not correct_key:
            raise ValueError(f"Experiment name {run_name} in .toml not found in channel-specification toml.")

        selected_channels = set(toml_dict["conditions"][correct_key]["channels"])
        logging.info("grabbing channel numbers ...")
        return selected_channels



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
    def launch_readfish(toml: str, device: str, name: str) -> None:
        """
        Wrapper to launch readfish into the background

        :param toml: TOML config for readfish
        :param device: Flowcell position
        :param name: Name of BOSS experiment
        :return:
        """
        # find the script to launch readfish
        module_path = inspect.getfile(LiveRun)
        logging.info(module_path)
        script_path = Path(module_path).parent / "readfish_boss.py"
        if not script_path.is_file():
            raise FileNotFoundError("boss_readfish.py not found. Something went wrong..")
        readfish_comm = f'python {script_path} {toml} {device} {name} 2>&1 | tee -a readfish.log'
        logging.info(readfish_comm)
        if device == "TEST":  # exit for testing purposes
            return
        if not Path(toml).exists():
            raise FileNotFoundError("Specified readfish toml does not exist.")
        # launch readfish into the background
        logging.info("Launching readfish")
        subprocess.Popen(readfish_comm, shell=True)



