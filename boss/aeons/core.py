import logging
import time
import os

from boss.core import Boss
from boss.live import LiveRun
from boss.batch import FastqBatch
from boss.utils import execute
from boss.aeons.sequences import SequencePool, SequenceAVA, UnitigPool, ContigPool, Benefit
from boss.aeons.repeats import RepeatFilter



class BossAeons(Boss):


    def init(self):
        """
        Initialisation function to be run after the init of the superclass
        :return:
        """
        # initialise central objects
        self.pool = SequencePool(name=self.args.name,
                                 min_len=self.args.min_seq_len,
                                 out_dir=self.out_dir)
        self.ava = SequenceAVA(paf=f'{self.args.name}.ava',
                               tetra=self.args.tetra,
                               filters=self.args)
        # init scoring func
        self.score_vec = Benefit.init_scoring_vec(lowcov=self.args.lowcov)
        # launch into first assembly
        if hasattr(self.args, "live_run") and getattr(self.args, "live_run"):
            self.first_live_asm()



    def first_live_asm(self) -> None:
        """
        Construct a first assembly from some initial data.
        Scan directory for data and wait until args.data_wait megabases accumulate
        Loop until some contigs have been generated

        :return:
        """
        while True:
            new_fastq = LiveRun.scan_dir(fastq_pass=self.args.fq, processed_files=set())
            fq_batch = FastqBatch(fq_files=new_fastq, channels=self.channels)
            logging.info(f"available so far: {fq_batch.total_bases / 1e6} Mbases")
            if fq_batch.total_bases / 1e6 < self.args.data_wait:
                logging.info(f"waiting for {self.args.data_wait} Mbases of data ... \n")
                time.sleep(30)
                continue
            else:
                # try first asm
                logging.info("Attempting initial assembly")
                init_pool = SequencePool(name="init_pool",
                                         min_len=self.args.min_seq_len,
                                         out_dir=self.out_dir)
                init_pool.ingest(seqs=fq_batch.read_sequences)

                init_contigs = init_pool.initial_asm_miniasm()
                ncontigs = len(init_contigs.sequences)
                has_contig = init_pool.has_min_one_contig(min_contig_len=self.args.min_contig_len)

                if not ncontigs or not has_contig:
                    logging.info(f"Initial assembly yielded no contigs, waiting for more data ... ")
                    time.sleep(30)
                    continue
                else:
                    self.pool = SequencePool(name=self.args.name,
                                             min_len=self.args.min_seq_len,
                                             out_dir=self.out_dir)
                    self.ava = SequenceAVA(paf=f'{self.args.name}.ava',
                                           tetra=self.args.tetra,
                                           filters=self.args)
                    self.pool.ingest(init_contigs)
                    # set up a repeat filter from the initial set of reads
                    if self.args.filter_repeats:
                        self.repeat_filter = RepeatFilter(name=self.args.name, seqpool=init_pool)
                    break

        # once there are contigs, record used files
        self.processed_files = set()
        self.processed_files.update(new_fastq)
        self.n_fastq = len(new_fastq)
        logging.info("Initial asm completed\n\n")



    def assemble_unitigs(self) -> SequencePool:
        """
        Wrapper to construct new contigs from the current data.
        Removes used sequences in the process

        :return: SequencePool of the new unitigs
        """
        # write current links to paf
        self.ava.links2paf(paf_out=self.ava.paf_links)
        # write pool to file
        SequencePool.write_seq_dict(seq_dict=self.pool.seqdict(), file=self.pool.fa)
        # create gfa and unitigs
        gfa = self.ava.paf2gfa_gfatools(paf=self.ava.paf_links, fa=self.pool.fa)
        # load the new unitigs
        unitigs = SequencePool.load_unitigs(gfa=gfa)
        # put them into a collection
        unitig_pool = UnitigPool(unitigs)
        # get the coverage arrays
        unitig_pool.get_unitig_coverage_arrays(seqpool=self.pool)
        # transform into a sequence pool
        new_pool, used_sids = unitig_pool.unitigs2seqpool(
            seqpool=self.pool, min_seq_len=self.args.min_seq_len
        )
        # remove used sequences from current pool
        self.remove_seqs(used_sids)
        return new_pool


    def assemble_add_and_filter_contigs(self) -> SequencePool:
        """
        Main wrapper to update the current contigs.
        Assembles the current graph and extracts new unitigs,
        adds them to the seqpool and extracts contigs for mapping against

        :return: SequencePool of current contigs (old and new)
        """
        logging.info("assembling new unitigs.. ")
        new_pool = self.assemble_unitigs()
        # add new sequences to the dict and to the ava
        logging.info("loading and overlapping new unitigs.. ")
        self.add_new_sequences(sequences=new_pool, increment=False)
        # write the current pool to file for mapping against
        logging.info("finding contigs to map against.. ")
        contigs = self.pool.declare_contigs(min_contig_len=self.args.min_contig_len)
        SequencePool.write_seq_dict(seq_dict=contigs.seqdict(), file=self.pool.contig_fa)
        return contigs




    def remove_seqs(self, sequences: set[str]) -> None:
        """
        Wrapper to remove sequences from pool, ava, coverage etc.

        :param sequences: Set of read IDs
        :return:
        """
        if not sequences:
            return
        self.ava.remove_links(sequences=sequences)
        self.pool.remove_sequences(sequences=sequences)



    def add_new_sequences(self, sequences: SequencePool, increment: bool = True) -> None:
        """
        Wrapper to add new sequences to the main pool

        :param sequences: Sequences to add
        :param increment: Whether to increment coverage when adding new sequences
        :return:
        """
        logging.info("Adding new sequences")
        ava_new, ava_onto_pool = self.pool.add2ava(sequences)
        self.pool.ingest(seqs=sequences)
        # load new alignments
        cont_new, ovl_new = self.ava.load_ava(ava_new, seqpool=self.pool)
        if increment:
            self.pool.increment(containment=cont_new)
        cont_onto, ovl_onto = self.ava.load_ava(ava_onto_pool, seqpool=self.pool)
        if increment:
            self.pool.increment(containment=cont_onto)
        # remove contained sequences
        cont = SequenceAVA.source_union(edges0=cont_new, edges1=cont_onto)
        self.remove_seqs(sequences=cont)
        # raise temp for overlappers (new links are saved as class attr in load_ava)
        ovl = ovl_new | ovl_onto
        self.pool.reset_temperature(ovl, t=self.args.temperature)



    def overlap_pool(self) -> None:
        """
        Wrapper to run AVA for pool to find overlaps and remove contained sequences

        :return:
        """
        logging.info("Running all-versus-all of sequence pool")
        contigs = self.pool.declare_contigs(min_contig_len=self.args.min_contig_len)
        if contigs.is_empty():
            return
        pool_paf = self.pool.run_ava(sequences=contigs.seqdict(), fa=self.pool.fa, paf=self.pool.ava)
        pool_contained, pool_ovl = self.ava.load_ava(paf=pool_paf, seqpool=self.pool)
        self.pool.increment(containment=pool_contained)
        cont = SequenceAVA.source_union(edges0=pool_contained, edges1={})
        if cont:
            logging.info(f'Removing {len(cont)} contained sequences from pool')
            self.remove_seqs(sequences=cont)
        self.pool.reset_temperature(pool_ovl)



    def trim_sequences(self) -> None:
        """
        Wrapper to find reads that need trimming for potential overlaps
        done after load_ava where mappings are marked with c=6
        i.e. only trim stuff already in the pool

        :return:
        """
        logging.info('')
        trim_dict = self.ava.to_be_trimmed()
        logging.info(f"Trimming {len(trim_dict.keys())} sequences")
        # trim and ingest
        trimmed_seqs = self.pool.trim_sequences(trim_dict=trim_dict)
        trim_paf = self.pool.run_ava(sequences=trimmed_seqs,
                                     fa=f'{self.pool.fa}.trim',
                                     paf=f'{self.pool.ava}.trim')
        trim_contained, _ = self.ava.load_ava(paf=trim_paf, seqpool=self.pool)
        to_remove = self.ava.trim_success(trim_dict=trim_dict, overlaps=self.ava.overlaps)
        # remove original sequences & failed mergers
        self.remove_seqs(sequences=to_remove)





    def cleanup(self) -> None:
        """
        Move temp files after a simulation

        :return:
        """
        tmpdir = f'{self.out_dir}/tmp'
        if not os.path.exists(tmpdir):
            os.mkdir(tmpdir)
        execute(f'mv {self.args.name}.* {tmpdir}')
        with open("sim.done", 'w') as done:
            done.write('')



    def update_wrapper(self, new_reads: dict[str, str]) -> None:
        """
        Shared aspects of updates of aeons experiment
        :return:
        """
        # filter sequences with repeats at the end
        if self.args.filter_repeats:
            reads_filtered = self.repeat_filter.filter_batch(seq_dict=new_reads)
        else:
            reads_filtered = new_reads

        # load new sequences, incl length filter
        sequences = SequencePool(sequences=reads_filtered, min_len=self.args.min_seq_len)
        # add new sequences to AVA
        self.add_new_sequences(sequences=sequences)
        # check for overlaps and containment in pool
        self.overlap_pool()
        # trim sequences that might lead to overlaps
        self.trim_sequences()
        # call wrapper to update assembly
        contigs = self.assemble_add_and_filter_contigs()
        contig_pool = ContigPool(sequences=contigs.sequences)
        # write the current pool to file for mapping against
        self.pool.write_seq_dict(seq_dict=contigs.seqdict(), file=self.pool.contig_fa)
        # check if we have any frozen sequences
        frozen_ids = self.pool.decrease_temperature(lim=self.args.min_contig_len)
        self.remove_seqs(sequences=frozen_ids)
        # update the sequencing masks
        self.strat = contig_pool.process_contigs(
            score_vec=self.score_vec,
            ccl=self.rl_dist.approx_ccl,
            out_dir=self.out_dir,
            lam=self.rl_dist.lam,
            batch=self.batch)



    def process_batch_aeons(self, new_reads: dict[str, str], **kwargs) -> None:
        """
        pass-through function for superclass method
        :param new_reads: Dictionary of new sequence data
        :param kwargs: Catches other arguments
        :return:
        """
        self.update_wrapper(new_reads=new_reads)



