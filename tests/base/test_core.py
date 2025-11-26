import pytest
import time
import logging

import boss.core



@pytest.fixture
def Boss(args_fake_real):
    b = boss.core.Boss(args=args_fake_real)
    b.launch_live_components()
    return b



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


