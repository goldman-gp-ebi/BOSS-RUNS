#!/usr/bin/env python3
from bossruns.BR_core import BossRun_live
from bossruns.BR_utils import setup_parser_live
from time import sleep
import logging

# set up argument parser
parser = setup_parser_live()
args = parser.parse_args()

# initialise constants and data processor
bossrun = BossRun_live(args=args)
logging.info("starting initialisation")
bossrun.initialise_OTUs()
bossrun.initialise_merged()

if args.ckp:
    bossrun.relaunch_checkpoint(ckp=args.ckp)

# main loop - periodically check for new data
logging.info('Initialisation completed.. waiting for sequencing data\n')
try:
    while True:
        next_update = bossrun.process_batch()
        # if processing was faster than the waiting time, sleep the rest
        if next_update > 0:
            sleep(next_update)

# stopping BR requires keyboard interrupt
except KeyboardInterrupt:
    logging.info("exiting after keyboard interrupt.. ")
