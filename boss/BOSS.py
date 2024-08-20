from time import sleep

import boss.config
from boss.runs import core as br_core
from boss.runs import simulation as br_sim
from boss.aeons import core as ba_core
from boss.aeons import simulation as ba_sim





def main(arg_list: list = None):
    # load arguments
    # DEBUG
    # arg_list = [
    # '--toml', "./config/BOSS_AEONS_SIM.toml"
    # ]
    arg_list = None if not arg_list else arg_list
    conf = boss.config.Config(parse=True, arg_list=arg_list)
    args = conf.args

    # launch entry-point
    # Real experiment mode
    if args.live_run:
        # initialise mode
        if args.ref:
            # runs initialises then launches readfish
            exp = br_core.BossRuns(args=args)
            exp.init()
            exp.launch_live_components()
            func = exp.process_batch_runs
        else:
            # aoens needs readfish running to have channel numbers available
            exp = ba_core.BossAeons(args=args)
            exp.launch_live_components()
            exp.init()
            func = exp.process_batch_aeons
        # launch main loop
        try:
            while True:
                next_update = exp.process_batch(func)
                # wait until next update
                if next_update > 0:
                    sleep(next_update)

        except KeyboardInterrupt:
            print("exiting after keyboard interrupt.. ")


    # Simulation mode
    else:
        # initialise mode
        if args.ref:
            exp = br_sim.BossRunsSim(args=args)
            func = exp.process_batch_runs_sim
        else:
            exp = ba_sim.BossAeonsSim(args=args)
            func = exp.process_batch_aeons_sim
        # launch main loop
        exp.init_sim()
        while exp.batch <= args.maxb:
            exp.process_batch_sim(func)
        exp.cleanup()



if __name__ == "__main__":
    main()


