import logging

from boss.sampler import Sampler
from boss.batch import ReadCache
from boss.mapper import Mapper
from boss.paf import Paf, paf_dict_type
from boss.aeons.core import BossAeons
from boss.aeons.sequences import SequencePool
from boss.aeons.repeats import RepeatFilter


class BossAeonsSim(BossAeons):



    def init_sim(self) -> None:
        """
        Initialisation wrapper for simulated experiments

        :return:
        """
        # run the init from super
        self.init()

        # initialise the sampler
        self.sampler = Sampler(
            source=self.args.fq,
            maxbatch=self.args.maxb,
            batchsize=self.args.batchsize
        )
        # initialise pseudotiming object
        self.read_cache = ReadCache(batchsize=self.args.batchsize, dumptime=self.args.dumptime)
        # run initial assembly with starting data
        self._initial_asm()
        self.strat = {}  # empty strat to accept everything




    def _initial_asm(self):
        # load some initial batches
        init_pool = SequencePool(name="init_pool", out_dir=self.out_dir)
        for i in range(self.args.binit):
            read_sequences, _, _, _ = self.sampler.sample()
            init_pool.ingest(seqs=read_sequences)
        logging.info(f"total bases in pool: {init_pool.total_bases()}")
        # increment time after preloading
        self.read_cache.update_times_aeons(read_sequences=init_pool.seqdict(),
                                           reads_decision=init_pool.seqdict())
        # set the batch counter
        self.batch = self.sampler.fq_stream.batch
        # initialise a repeat filter from the raw initial reads
        if self.args.filter_repeats:
            self.repeat_filter = RepeatFilter(name=self.args.name, seqpool=init_pool)
        # run first asm
        logging.info("Running assembly of initial data..")
        init_contigs = init_pool.initial_asm_miniasm()
        self.pool.ingest(init_contigs)
        has_contig = self.pool.has_min_one_contig(min_contig_len=self.args.min_contig_len)
        ncontigs = len(self.pool.sequences)
        logging.info(f'initial contigs: {ncontigs}')
        if len(self.pool.sequences) == 0 or not has_contig:
            raise ValueError("No contigs of sufficient length, restart simulation with more data (binit)")
        self.pool.write_seq_dict(seq_dict=self.pool.seqdict(), file=self.pool.contig_fa)



    def make_decisions(
        self,
        paf_dict: paf_dict_type,
        read_sequences: dict,
        window: int = 100,
        mu: int = 400
    ) -> dict:
        """
        Decision function for simulations only.
        In real experiments readfish performs this functionality

        :param paf_dict: Paf dict of mapping runcated reads to contigs
        :param read_sequences: Dictionary of new read sequences
        :param window: Size of downsampling window
        :return: Dict of read sequences after decisions, i.e. rejected reads are truncated
        """
        # if nothing mapped, just return. Unmapped = accept
        if len(paf_dict.items()) == 0:
            logging.info("nothing mapped")
            self.reject_count = 0
            self.accept_count = 0
            self.unmapped_count = 0
            return read_sequences

        reads_decision = dict()
        reject_count = 0
        accept_count = 0
        unmapped_count = 0
        # loop over paf dictionary
        for record_id, record_list in paf_dict.items():
            # record_id, record_list = list(gaf_dict.items())[0]
            if len(record_list) > 1:
                # should not happen often since we filter secondary mappings
                rec = Paf.choose_best_mapper(record_list)[0]
            else:
                rec = record_list[0]

            # find the start and end position relative to the whole linearised genome
            if rec.rev:
                start_pos = rec.tend - 1
            else:
                start_pos = rec.tstart

            # index into strategy to find the decision
            try:
                decision = self.strat[str(rec.tname)][start_pos // window, rec.rev]
            except:
                # if we don't have a strategy yet accept all
                decision = 1

            # ACCEPT
            if decision:
                record_seq = read_sequences[rec.qname]
                accept_count += 1

            # REJECT
            else:
                record_seq = read_sequences[rec.qname][: mu]
                reject_count += 1

            # append the read's sequence to a new dictionary of the batch after decision-making
            reads_decision[rec.qname] = record_seq

        # all unmapped reads also need to be accepted, i.e. added back into the dict
        mapped_ids = set(reads_decision.keys())

        for read_id, seq in read_sequences.items():
            if read_id in mapped_ids:
                continue
            else:
                reads_decision[read_id] = seq
                unmapped_count += 1

        logging.info(f'decisions - rejecting: {reject_count} accepting: {accept_count} unmapped: {unmapped_count}')
        self.reject_count = reject_count
        self.accept_count = accept_count
        self.unmapped_count = unmapped_count
        return reads_decision





    def process_batch_aeons_sim(self) -> None:
        """
        Process new batch of simulated data. Get new reads, make decisions,
        update times and write to files before launching the common steps of
        processing new data

        :return:
        """
        read_sequences, _, _, _ = self.sampler.sample()
        # initialise a new mapper for the current contigs
        lm = Mapper(ref=self.pool.contig_fa, default=False)
        # map the truncated version of reads to the current contigs
        paf_trunc = lm.map_sequences(sequences=read_sequences, trunc=True)
        # Dict of processed reads, i.e. rejected sequences are truncated
        reads_decision = self.make_decisions(paf_dict=paf_trunc, read_sequences=read_sequences)
        # update read length dist
        self.rl_dist.update(read_lengths=self.sampler.fq_stream.read_lengths)
        # updates to the pseudotimes
        self.read_cache.update_times_aeons(read_sequences=read_sequences, reads_decision=reads_decision)
        self.read_cache.fill_cache(read_sequences=read_sequences, reads_decision=reads_decision)
        # run the rest of the steps to update
        self.update_wrapper(new_reads=reads_decision)





