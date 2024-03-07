import pytest
import time
import logging

import boss.core
import boss.config




@pytest.fixture
def args():
    conf = boss.config.Config()
    # assign some args since we don't load the full config
    conf.args.toml_readfish = "TEST"
    conf.args.split_flowcell = False
    conf.args.live_run = True
    return conf.args


@pytest.fixture
def Boss(args):
    return boss.core.Boss(args=args)



def test__time_to_next_update(Boss):
    tic = time.time()
    assert type(Boss._time_to_next_update(tic)) is int
    assert Boss._time_to_next_update(tic) < 1e3
    assert Boss._time_to_next_update(tic) > 0



def test_process_batch(Boss):
    # dummy function to pass into process_batch
    def dummy(new_reads, new_quals):
        logging.info(f'new_reads: {len(new_reads)}, new_quals: {len(new_quals)}')

    next_update = Boss.process_batch(main_processing_func=dummy)
    assert next_update > 10
    assert Boss.batch == 1
    # launch it again to trigger the case when there are no new files
    next_update = Boss.process_batch(main_processing_func=dummy)
    assert next_update == Boss.args.wait
    assert Boss.batch == 1



def test_process_batch_sim(Boss):
    # dummy function to pass into process_batch
    def dummy():
        logging.info(f'simulation function goes here')

    next_update = Boss.process_batch_sim(main_processing_func=dummy)
    assert next_update > 10
    assert Boss.batch == 1


