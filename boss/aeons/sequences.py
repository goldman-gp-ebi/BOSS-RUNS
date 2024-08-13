import logging
import os
from collections import defaultdict, Counter
from copy import deepcopy
from shutil import copy
import time
from types import SimpleNamespace
from typing import Optional, Any
from concurrent.futures import ThreadPoolExecutor as TPexe
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
import bottleneck as bn

from boss.utils import execute, random_id, load_gfa, write_logs
from boss.dependencies import Dependencies
from boss.paf import Paf, PafLine
from boss.aeons.kmer import euclidean_dist, euclidean_threshold


Edge = tuple[str, str]  # typehint for containments: tuple of (source, target)



class SequenceAVA:

    def __init__(self, paf: str, filters: SimpleNamespace, tetra: bool = False):
        """
        Initialize the SequenceAVA object.

        :param paf: Path to the PAF file.
        :param filters: filters to apply, stored as arguments.
        :param tetra: Whether to use the tetramer distance to filter overlaps. Defaults to False.
        """
        self.paf = paf
        self.gfa = f'{paf}.gfa'
        self.filters = filters
        self.tetra = tetra
        # container to keep overlaps
        # save paf lines as qname-tname with alphanum sort
        self.links = defaultdict(lambda: defaultdict(PafLine))
        self.paf_links = f"{paf}.links.paf"
        self.dep = Dependencies()


    def load_ava(self, paf: str, seqpool: "SequencePool") -> tuple[dict[Edge, PafLine], set]:
        """
        Load all entries from a PAF file as PafLines and filter the entries while loading.

        :param paf: Path to the PAF file.
        :param seqpool: SequencePool
        :return: A tuple containing containments and overlappers.
        """
        self.trims = []  # used for trimming
        self.overlaps = {}  # used for trimming
        containments = {}  # collect, used for coverage increments
        overlappers = set()  # used to track temperature of sequences
        ovl = 0
        inter = 0

        records, skip = Paf.parse_filter_classify_records(paf=paf, filters=self.filters)

        for rec in records:
            if rec.c == 2:
                # first contained
                # check if we have a containment of the two already
                if (rec.qname, rec.tname) in containments.keys():
                    if rec.s1 < containments[(rec.qname, rec.tname)].s1:
                        # previously recorded containment has higher s1
                        continue
                containments[(rec.qname, rec.tname)] = rec
            elif rec.c == 3:
                # second contained
                # check if we have a containment of the two already
                if (rec.tname, rec.qname) in containments.keys():
                    if rec.s1 < containments[(rec.tname, rec.qname)].s1:
                        # previously recorded containment has higher s1
                        continue
                containments[(rec.tname, rec.qname)] = rec
            elif rec.c in {4, 5}:
                # if applicable, check for tetramer dist
                if self.tetra:
                    intra = seqpool.is_intra(rec.qname, rec.tname)
                    if not intra:
                        inter += 1
                        continue

                if not seqpool.sequences[rec.tname].acceptor:
                    rec.c = 2
                    containments[(rec.qname, rec.tname)] = rec
                    continue
                if not seqpool.sequences[rec.qname].acceptor:
                    rec.c = 3
                    containments[(rec.tname, rec.qname)] = rec
                    continue

                ovl += 1
                self.overlaps[(rec.qname, rec.tname)] = rec
                # check if we have an overlap of the two already
                if rec.tname in self.links[rec.qname].keys():
                    if rec.s1 < self.links[rec.qname][rec.tname].s1:
                        # previously recorded overlap has higher s1
                        continue

                self.links[rec.qname][rec.tname] = rec
                self.links[rec.tname][rec.qname] = rec
                overlappers.add(rec.qname)
                overlappers.add(rec.tname)
            elif rec.c == 6:
                self.trims.append(rec)
            else:
                pass

        contained_ids = set([i for (i, j) in containments.keys()])
        # remove contained ones from further containment search
        skip_filt = [s for s in skip if s.qname not in contained_ids and s.tname not in contained_ids]
        # run multiline containment detection
        mc = MultilineContainments(records=skip_filt)
        containments.update(mc.containments)

        logging.info(f"ava load: skip {len(skip)}, cont {len(contained_ids)} cont multi: {len(mc.containments)} ovl: {ovl} inter: {inter}")
        return containments, overlappers


    def remove_links(self, sequences: set[str]) -> None:
        """
        Remove overlaps of certain sequences.

        :param sequences: Set of sequence IDs to remove overlaps from.
        """
        for sid in sequences:
            # targets for overlaps where sid is "query"
            targets = self.links[sid].keys()
            # remove the overlaps where sid is "query"
            self.links.pop(sid, None)
            # remove overlaps from targets
            for t in targets:
                self.links[t].pop(sid, None)


    def to_be_trimmed(self) -> dict[str, tuple[int, int, Any]]:
        """
        After classification, find the coordinates for trimming.

        :return: A dictionary containing sequences to be trimmed and their corresponding coordinates.
        """
        trim = self.trims
        # save the name and coordinates for trimming
        to_trim = {}
        for rec in trim:
            sid, trim_start, trim_stop, other = rec.find_trim_coords()
            if sid == '0':  # don't trim if product would be short
                continue
            to_trim[sid] = (trim_start, trim_stop, other)
        return to_trim


    @staticmethod
    def trim_success(trim_dict: dict[str, tuple[int, int, Any]], overlaps: dict[tuple[str, str], Any]) -> set[str]:
        """
        Check which trimming attempts were successful.

        :param trim_dict: Dictionary containing sequences to be trimmed and their coordinates.
        :param overlaps: Dictionary of overlaps.
        :return: A set of sequences to be removed.
        """
        success, unsuccess = set(), set()
        trim = set(trim_dict.keys())
        if not trim:
            to_remove = success | unsuccess
            return to_remove
        if not overlaps:
            unsuccess = {f'{t}%' for t in trim}
            to_remove = success | unsuccess
            return to_remove

        ovl_q, ovl_t = zip(*overlaps.keys())
        ovl = set(ovl_q) | set(ovl_t)
        trim_mod = {f'{t}%' for t in trim}
        success_raw = trim_mod & ovl
        unsuccess = trim_mod - success_raw
        # remove the percentage marker for removal
        success = {s[:-1] for s in success_raw}
        # merge the headers for removal
        to_remove = success | unsuccess
        return to_remove


    def links2paf(self, paf_out: str) -> None:
        """
        Write overlaps to a PAF file.

        :param paf_out: Output file path for the PAF file.
        """
        written = set()
        with open(paf_out, 'w') as fh:
            for node, target_dict in self.links.items():
                if not target_dict:
                    continue
                # go through all links
                for target, rec in target_dict.items():
                    if rec.line in written:
                        continue
                    else:
                        fh.write(rec.line)
                        written.add(rec.line)


    def paf2gfa_gfatools(self, paf: str, fa: str, gfa: Optional[str] = None) -> str:
        """
        Use gfatools to generate the graph and unitigs from a PAF file and a FASTA file.

        :param paf: Input PAF file path.
        :param fa: Input FASTA file path.
        :param gfa: Optional output GFA file path.
        :return: The generated GFA content as a string or empty string
        """
        if not os.path.getsize(paf):
            logging.info("no overlaps for merging")
            return ""
        paf2gfa = getattr(self.dep, "paf2gfa")

        comm = f"{paf2gfa} -i {fa} -u -c {paf}"
        stdout, stderr = execute(comm)
        # optionally write to file
        if gfa:
            with open(gfa, 'w') as fh:
                fh.write(stdout)
        return stdout


    @staticmethod
    def source_union(edges0: dict[Edge, PafLine], edges1: dict[Edge, PafLine]) -> set:
        """
        Given containment dictionaries (edges), return the union of source nodes.

        :param edges0: First containment dictionary.
        :param edges1: Second containment dictionary
        :return: Union of source nodes.
        """
        sources0 = zip(*edges0.keys())
        sources1 = zip(*edges1.keys())
        try:
            set0 = set(list(sources0)[0])
        except IndexError:
            set0 = set()

        try:
            set1 = set(list(sources1)[0])
        except IndexError:
            set1 = set()

        source_union = set0 | set1
        return source_union




class Sequence:

    def __init__(
        self,
        header: str,
        seq: str,
        cov: Optional[NDArray] = None,
        merged_components: Optional[set[str]] = None,
        merged_atoms: Optional[set[str]] = None,
        cap_l: bool = False,
        cap_r: bool = False,
    ):
        """
        Initialize a Sequence object.

        :param header: The header of the sequence.
        :param seq: The sequence data.
        :param cov: The coverage array (optional).
        :param merged_components: The list of merged components (optional).
        :param merged_atoms: The list of merged atoms (optional).
        :param cap_l: Flag indicating if the sequence has a left cap (default: False).
        :param cap_r: Flag indicating if the sequence has a right cap (default: False).
        """

        self.header = header
        self.seq = seq

        if cov is None:
            self.cov = np.ones(shape=len(seq), dtype='float')
        else:
            self.cov = cov

        # merged_headers provides the info which reads were used to create a new sequence
        if merged_components is None:
            self.atoms = set()  # includes all sequence reads that are contained in this one
            self.components = set()  # all sequences used to build this one
        else:
            # if merged info are provided
            self.components = set(merged_components)
            self.atoms = set(merged_atoms)
        # inits
        self.tetramer_zscores = 0
        # self.kmers = 0
        # temperature for ignoring reads
        # if temperature reaches 0, reads gets frozen
        self.temperature = 30
        self.cap_l = cap_l
        self.cap_r = cap_r
        # accept overlaps
        self.acceptor = True


    def is_hot(self) -> bool:
        """
        Is the sequence active?

        :return: True if temperature is >0, False otherwise.
        """
        if self.temperature > 0:
            return True
        else:
            return False


    def chunk_up_coverage(self, n: int):
        """
        Chunk up the coverage array to reduce resolution of strategies.

        :param n: The chunk size.
        """
        self.cov_chunked = np.array([np.sum(self.cov[i: i + n]) for i in range(0, len(self.cov), n)])
        # init an empty array to record nodes of interest
        self.noi = np.zeros(shape=self.cov_chunked.shape[0], dtype="bool")
        self.scores = np.zeros(shape=self.cov_chunked.shape[0], dtype="float")
        self.benefit = np.zeros(shape=self.cov_chunked.shape[0], dtype="float")


    def contig_scores(self, score_vec: NDArray, n: int = 100):
        """
        Calculate scores for the contig.

        :param score_vec: The precomputed scoring vector.
        :param n: The node size.
        """
        sc = Benefit.score_array(score_vec=score_vec, cov_arr=self.cov_chunked, node_size=n)
        self.scores = sc
        assert self.cov_chunked.shape[0] == self.scores.shape[0]


    def contig_benefits(self, mu: int, ccl: NDArray, node_size: int = 100):
        """
        Calculate benefits for the contig.

        :param mu: Length of anchor bases
        :param ccl: The ccl value.
        :param node_size: The node size.
        """
        benefit, smu_sum = Benefit.calc_fragment_benefit(
            scores=self.scores,
            mu=mu,
            approx_ccl=ccl,
            node_size=node_size,
            e1=bool(self.noi[0]),
            e2=bool(self.noi[-1])
        )
        self.benefit = benefit
        self.smu_sum = smu_sum
        assert self.cov_chunked.shape[0] == self.benefit.shape[1]


    def set_contig_ends(self, n: int, lim: int = 50):
        """
        Declare ends of contig.

        :param n: The node size.
        :param lim: The coverage limit (default: 50).
        """
        cc = self.cov_chunked
        if cc[0] > lim * n:
            pass
        elif self.cap_l:
            pass
        # only set as interesting if every check fails
        # for either model: mark as noi or give max score
        else:
            self.noi[0] = 1
            self.scores[0] = 1

        if cc[-1] > lim * n:
            pass
        elif self.cap_r:
            pass
        else:
            self.noi[-1] = 1
            self.scores[-1] = 1


    def find_strat_m0(self, threshold: float) -> NDArray:
        """
        Find the strategy for the contig using benefit threshold.

        :param threshold: threshold value.
        :return: Strategy array.
        """
        strat = np.where(self.benefit >= threshold, True, False)
        return strat.transpose()




class SequencePool:
    def __init__(
        self,
        sequences: Optional[dict[str, str | Sequence]] = None,
        name: str = "dummy",
        min_len: int = 3000,
        out_dir: str = "dummy",
        threads: int = 6
    ):
        """
        Initialize a SequencePool object. Unified pool for reads and contigs with persistent AVA

        :param sequences: Dictionary of raw sequences. Keys are headers, values are either raw seqs or Sequence objects
        :param name: Name of the sequence pool.
        :param min_len: Minimum length of sequences to be included in the pool. Default is 3000.
        :param out_dir: Output directory.
        :param threads: Number of threads to use. Default is 6.
        """
        self.min_len = min_len
        self.out_dir = out_dir
        self.sequences = dict()
        self.threads = threads
        self.name = name

        if sequences:
            input_type = type(list(sequences.values())[0])
            # raw sequences
            if input_type == str:
                self._ingest_dict(seqs=sequences)
            # Sequence objects
            elif input_type == Sequence:
                self.sequences = sequences
            else:
                print("SequencePool input type not supported")

        self.polished = {}
        # filenames
        self.fa = f'{name}.fa'  # fasta of whole pool
        self.contig_fa = f'{name}.contig.fa'  # fasta of long sequences to map against
        self.ava = f'{name}.ava'  # ava in paf
        self.gfa = f'{name}.gfa'
        self.dep = Dependencies()


    def headers(self) -> set[str]:
        """
        Get the set of sequence headers in the SequencePool.

        :return: Set of sequence headers.
        """
        return set(self.sequences.keys())


    def seqdict(self) -> dict[str, str]:
        """
        Get a dictionary representation of the raw sequences in the SequencePool.

        :return: Dictionary of sequence headers and their corresponding raw sequences.
        """
        return {header: seqo.seq for header, seqo in self.sequences.items()}


    def total_bases(self) -> int:
        """
        Get the total number of bases in the SequencePool.

        :return: Total number of bases.
        """
        return int(np.sum([len(seqo.seq) for header, seqo in self.sequences.items()]))


    def is_empty(self) -> bool:
        """
        Check if the pool is empty

        :return: True if no sequences in pool, False otherwise
        """
        empty = True if len(self.sequences) == 0 else False
        return empty


    def ingest(self, seqs: Any) -> None:
        """
        Add a pile of sequences to the SequencePool.

        :param seqs: Dictionary of raw sequences or an existing SequencePool to add.
        """
        if type(seqs) is dict:
            skipped = self._ingest_dict(seqs=seqs)
            logging.info(f"ingested: {len(seqs) - skipped} pool size: {len(self.sequences.keys())}")
        elif type(seqs) is SequencePool:
            self._ingest_pool(new_pool=seqs)
            logging.info(f"ingested: {len(seqs.sequences)} pool size: {len(self.sequences.keys())}")
        else:
            logging.info("seqs need to be dict, or SequencePool")


    def _ingest_dict(self, seqs: dict[str, str]) -> int:
        """
        Ingest a dictionary of raw sequences.

        :param seqs: Dictionary of raw sequences.
        :return: Number of skipped sequences.
        """
        skipped = 0
        for rid, seq in seqs.items():
            if len(seq) > self.min_len:
                # init sequence object without existing arrays
                seqo = Sequence(header=rid, seq=seq)
                self.sequences[rid] = seqo
            else:
                skipped += 1
        return skipped


    def _ingest_pool(self, new_pool: 'SequencePool') -> None:
        """
        Ingest sequences from an existing SequencePool, e.g. after merging Sequences

        :param new_pool: Existing SequencePool to ingest sequences from.
        """
        for rid, seqo in new_pool.sequences.items():
            if len(seqo.seq) > self.min_len:
                self.sequences[rid] = seqo


    def run_ava(self, sequences: dict[str, str], fa: str, paf: str, base_level: bool = False) -> str:
        """
        Perform All-Versus-All (AVA) comparison of sequences using minimap2.

        :param sequences: Dictionary of sequences.
        :param fa: Filename for temporary sequence file.
        :param paf: Filename for the output PAF file.
        :param base_level: Whether to perform base-level AVA. Default is False.
        :return: Filename of PAF file.
        """
        logging.info(f"Running ava for: {len(sequences)} queries")
        # write current pool to file
        self.write_seq_dict(seq_dict=sequences, file=fa)
        mm2 = getattr(self.dep, "minimap2")

        comm = f'{mm2} -x ava-ont -t{self.threads} {fa} {fa} >{paf}'
        if base_level:
            comm = f'{mm2} -cx ava-ont -t{self.threads} {fa} {fa} >{paf}'

        tic = time.time()
        stdout, stderr = execute(comm)
        toc = time.time()
        logging.info(f"t: {toc - tic}")
        if os.path.exists(f'{self.out_dir}/logs'):
            write_logs(stdout, stderr, f'{self.out_dir}/logs/ava')
        return paf


    def initial_asm_miniasm(self, c: int = 2, trds: int = 8) -> 'SequencePool':
        """
        Perform initial assembly using miniasm.

        :param c: Minimum coverage, fwd to miniasm. Default is 3.
        :param trds: Number of threads. Default is 8.
        :return: SequencePool object containing assembled unitigs.
        """
        sfile = f'{self.name}.init_reads.fa'
        self.write_seq_dict(self.seqdict(), file=sfile)
        mm2 = getattr(self.dep, "minimap2")
        miniasm = getattr(self.dep, "miniasm")

        comm0 = f'{mm2} -x ava-ont -t{trds} {sfile} {sfile} >{sfile}.ava'
        stdout, stderr = execute(comm0)
        write_logs(stdout, stderr, f'{self.out_dir}/contigs/init/init_ava')
        comm1 = f"{miniasm} -f {sfile} {sfile}.ava -c{c} >{self.name}.init.gfa"
        stdout, stderr = execute(comm1)
        write_logs(stdout, stderr, f'{self.out_dir}/contigs/init/init_asm')

        contigs = load_gfa(f'{self.name}.init.gfa')
        contig_pool = SequencePool(contigs)
        # disallow extension of circular contigs
        for header, seqo in contig_pool.sequences.items():
            if header[-1] == 'c':
                seqo.acceptor = False
        return contig_pool


    def add2ava(self, new_sequences: 'SequencePool') -> tuple[str, str]:
        """
        Add new sequences to the All-Versus-All (AVA) comparisons.
        Instead of rerunning complete ava, run ava for new sequences and map new sequences onto existing ones

        :param new_sequences: SequencePool object containing the new sequences.
        :return: Filenames of the new AVA file and the file containing mappings of new seqs to previous ones.
        """
        logging.info(f'adding to ava: {len(new_sequences.sequences)}')
        self.write_seq_dict(seq_dict=self.seqdict(), file=self.fa)
        # declare new filenames
        new_fa = f'{self.fa}.new'
        new_ava = f'{self.ava}.new'
        new_onto_pool = f'{self.fa}.new.onto_pool'
        # write new reads to file
        self.write_seq_dict(seq_dict=new_sequences.seqdict(), file=new_fa)
        # ava of new sequences
        mm2 = getattr(self.dep, "minimap2")
        comm = f'{mm2} -x ava-ont -t{self.threads} {new_fa} {new_fa} >{new_ava}'
        stdout, stderr = execute(comm)
        write_logs(stdout, stderr, f'{self.out_dir}/logs/ava_add')
        # mapping new sequences to previous pool
        comm = f'{mm2} -x map-ont -w5 -e0 -m100 -r2k -t{self.threads} {self.fa} {new_fa} >{new_onto_pool}'
        stdout, stderr = execute(comm)
        write_logs(stdout, stderr, f'{self.out_dir}/logs/map2pool')
        # return filenames to be ingested as AVA
        return new_ava, new_onto_pool


    def remove_sequences(self, sequences: set[str]) -> None:
        """
        Remove sequences from the pool based on their IDs.
        Once a read is used, it can not contribute to other alignments

        :param sequences: Set of sequence IDs to remove.
        """
        pre = len(self.sequences)
        popped = 0
        for sid in sequences:
            self.sequences.pop(sid, None)
            popped += 1
        post = len(self.sequences)
        logging.info(f'Removed: {popped} ({pre} {post})')


    def trim_sequences(self, trim_dict: dict[str, tuple[int, int, str]]) -> dict[str, str]:
        """
        Trim sequences in order to create valid alignments.

        :param trim_dict: Dictionary of [sid: (start:stop, other_name)] for trimming sequences.
        :return: Dictionary of trimmed sequences.
        """
        trimmed_seqs = {}
        other_seqs = {}
        valid_ids = set()

        for sid, (start, stop, other) in trim_dict.items():
            try:
                nsid = sid + '%'
                trimmed_seqs[nsid] = deepcopy(self.sequences[sid])
                other_seqs[other] = self.sequences[other]
                valid_ids.add(nsid)
            except KeyError:
                logging.info("key for trimming not in sequence pool")
                continue

        for sid, (start, stop, other) in trim_dict.items():
            nsid = sid + '%'
            # skip sequences that were not found in previous step
            if nsid not in valid_ids:
                continue
            seqo = trimmed_seqs[nsid]
            # deselect the trimmed bit with a boolean mask
            mask = np.ones(shape=len(seqo.seq), dtype='bool')
            mask[start: stop] = 0
            seq_arr = np.array(list(seqo.seq))
            seq_arr_trim = seq_arr[mask]
            seq_trim = ''.join(seq_arr_trim)
            # replace attributes of the Sequence Obj
            seqo.seq = seq_trim
            seqo.cov = seqo.cov[mask]
            seqo.header = nsid

        trimmed_pool = SequencePool(sequences=trimmed_seqs)
        other_pool = SequencePool(sequences=other_seqs)
        # ingest pool of trimmed sequences
        self.ingest(seqs=trimmed_pool)

        # combine sequence dicts for mapping
        seq_dict = dict(trimmed_pool.seqdict(), **other_pool.seqdict())
        return seq_dict


    @staticmethod
    def get_next_increment_edges(edges: set[Edge], previous_edges: set[Edge] = None) -> tuple[set[Edge], set[Edge]]:
        """
        Get the next increment edges based on the given edges and previous edges.
        if no previous edges are given, get the edges with in-degree of 0: This is the start of the algorithm

        :param edges: Set of edges.
        :param previous_edges: Set of previous edges.
        :return: Tuple containing the remaining edges and the next increment edges.
        """
        if not previous_edges:
            sources, targets = zip(*edges)
            next_sources = set(sources) - set(targets)
        else:
            # otherwise grab the edges starting at the previous targets
            next_sources = [t for (s, t) in previous_edges]
        # get the next edges to increment
        next_edges = {(s, t) for (s, t) in edges if s in next_sources}
        # remove the next edges
        for e in next_edges:
            edges.remove(e)
        return edges, next_edges


    def effect_increment(self, source: str, target: str, rec: PafLine, edge_multiplicity: float) -> None:
        """
        Effect the increment by adjusting the coverage and adding the source as a constituent of the target.

        :param source: Source sequence ID.
        :param target: Target sequence ID.
        :param rec: Mapping record.
        :param edge_multiplicity: Edge multiplicity.
        """
        # grab relevant coordinates of this containment
        ostart, oend, olen, cstart, cend, clen = rec.grab_increment_coords()
        # grab source coverage
        cont_cov = self.sequences[source].cov[cstart: cend].copy()
        # adjust length of coverage array
        if clen == olen:
            pass
        elif clen > olen:
            cont_cov = cont_cov[: olen]
        elif clen < olen:
            cont_cov = np.pad(cont_cov, (0, olen - clen), mode='edge')

        # account for reverse complement
        if rec.rev:
            cont_cov = cont_cov[::-1]
        else:
            pass

        # adjust for edge multiplicity
        cont_cov /= edge_multiplicity

        # add the source coverage to target coverage
        self.sequences[target].cov[ostart: oend] += cont_cov
        # limit coverage to 100 to prevent mess
        self.sequences[target].cov[np.where(self.sequences[target].cov > 100)] = 100

        # add the source as constituent of target
        if '*' not in source:
            self.sequences[target].atoms.add(source)


    def effect_increments(self, next_edges: set[Edge], containment: dict[Edge, PafLine], edge_multiplicity: dict[str, float] = None) -> None:
        """
        Effect the increments based on the next set of edges.

        :param next_edges: Set of next edges.
        :param containment: Containment records.
        :param edge_multiplicity: Edge multiplicity.
        """
        # loop over the next set of edges to effect the collected increments
        for (source, target) in next_edges:
            rec = containment[(source, target)]
            # grab edge multiplicity modifier, in case of parallel edges
            if edge_multiplicity is None:
                em = 1
            else:
                em = edge_multiplicity[source]
            self.effect_increment(source, target, rec, em)


    @staticmethod
    def find_edge_multiplicity(edges: set[Edge]):
        """
        Find the edge multiplicity for the given edges.

        :param edges: Set of edges.
        """
        sources, targets = zip(*edges)
        source_counts = Counter(sources)
        return source_counts


    def increment(self, containment: dict[Edge, PafLine]) -> set:
        """
        Increment the coverage counts based on containment records.

        :param containment: Containment records.
        :return: Set of IDs of the contained sequences.
        """
        edges = set(containment.keys())
        if not edges:
            # nothing to do
            return set()

        # debugging
        # import sys
        # sys.path.insert(0, "/home/lukas/Desktop/Aeons/code/plot")
        # from gt_plot import containment_graph
        # containment_graph(edges, 0)

        # get the first edges to increment, i.e. those with 0 in-degree
        edges, next_edges = self.get_next_increment_edges(edges, previous_edges=None)
        if not next_edges:
            return set()
        edge_multiplicity = self.find_edge_multiplicity(next_edges)
        self.effect_increments(next_edges, containment, edge_multiplicity)
        previous_edges = next_edges

        while len(edges) > 0:
            # get the next edges
            edges, next_edges = self.get_next_increment_edges(edges, previous_edges=previous_edges)
            if not next_edges:
                return set()
            edge_multiplicity = self.find_edge_multiplicity(next_edges)
            self.effect_increments(next_edges, containment, edge_multiplicity)

            # circular containment relationships could trap us here
            if len(next_edges) == len(previous_edges):
                break
            previous_edges = next_edges

        # for removal: return the ids of the contained sequences
        contained_ids = set([s for (s, t) in containment.keys()])
        return contained_ids


    def reset_temperature(self, sids: set[str], t: int = 50) -> None:
        """
        Reset the temperature of active reads to the given value.

        :param sids: Set of sequence IDs.
        :param t: Temperature value.
        """
        # give active reads a boost in temperature
        for s in sids:
            try:
                seqo = self.sequences[s]
                seqo.temperature = t
            except KeyError:
                pass


    def decrease_temperature(self, lim: int) -> set[str]:
        """
        Decrease the temperature of all reads and return the frozen sequences.
        Sequences longer than lim never freeze

        :param lim: Length limit.
        :return: Set of frozen sequence headers.
        """
        frozen_seqs = set()
        for header, seqo in self.sequences.items():
            if len(seqo.seq) < lim:
                seqo.temperature -= 1
                if not seqo.is_hot():
                    frozen_seqs.add(header)
        logging.info(f"frozen seqs: {len(frozen_seqs)}")
        return frozen_seqs


    def declare_contigs(self, min_contig_len: int) -> 'SequencePool':
        """
        Declare contigs based on a minimum contig length. These will be used for strategy generation

        :param min_contig_len: Minimum contig length.
        :return contig_pool: Pool of contigs
        """
        contigs = {header: seqo for header, seqo in self.sequences.items() if len(seqo.seq) > min_contig_len}
        contig_pool = SequencePool(sequences=contigs)
        return contig_pool


    def has_min_one_contig(self, min_contig_len: int) -> bool:
        """
        Declare contigs based on a minimum contig length. These will be used for strategy generation

        :param min_contig_len: Minimum contig length.
        :return contig: bool
        """
        contigs = {header: seqo for header, seqo in self.sequences.items() if len(seqo.seq) > min_contig_len}
        contig = True if contigs else False
        return contig


    def get_atoms(self, headers: list) -> set[str]:
        """
        Get all atomic reads based on the given headers.

        :param headers: List of headers.
        :return: Set of atomic reads.
        """
        # given a list of headers, get all atomic reads
        atoms = set()
        for h in headers:
            atm = self.sequences[h].atoms
            atoms.update(atm)
        return atoms


    def get_components(self, headers: list) -> set[str]:
        """
        Get all component reads based on the given list of headers.

        :param headers: List of headers.
        :return: Set of component reads.
        """
        components = set()
        for h in headers:
            cmp = self.sequences[h].components
            components.update(cmp)
            components.add(h)
        return components


    @staticmethod
    def write_seq_dict(seq_dict: dict[str, str], file: str) -> None:
        """
        Write a sequence dictionary to file.

        :param seq_dict: Sequence dictionary.
        :param file: File to write the dictionary to.
        """
        with open(file, 'w') as fasta:
            for sid, seq in seq_dict.items():
                fasta.write(f'>{sid}\n')
                fasta.write(f'{seq}\n')



    @staticmethod
    def load_unitigs(gfa: str) -> list['Unitig']:
        """
        Load unitigs after graph cleaning with gfatools.
        Input can be a gfa file or the content of the file

        :param gfa: Path to the GFA file or a string representation of the GFA.
        :return: A list of Unitig objects.
        """
        if not gfa:
            return []

        if gfa.startswith('S\t'):
            gfat = gfa
        else:
            with open(gfa, 'r') as fh:
                gfat = fh.read()
        # string operations
        gfatt = gfat.split('\nS\t')
        # first and last need extra treatment
        gfatt[0] = gfatt[0][2:]  # skip the first "S\t"
        gfatt[-1] = gfatt[-1].replace('\nx', '\nL').split('\nL')[0]  # exclude all L and x lines
        # get the x lines
        gfax = gfat.split('\nx\t')[1:]
        # make sure we have the same number of S_A lines and x lines
        assert len(gfatt) == len(gfax)
        # parse into unitig objects
        unitigs = []
        for sa_lines, x_line in zip(gfatt, gfax):
            unitigs.append(Unitig(sa_lines.split('\n'), x_line))
        return unitigs


    def is_intra(self, seq1: str, seq2: str) -> bool:
        """
        Check if two sequences are classified as intraspecific based on the Euclidean distance of tetramers.

        :param seq1: Name of first sequence.
        :param seq2: Name of second sequence.
        :return: True if the sequences are classified as intraspecific, False otherwise.
        """
        s1 = self.sequences[seq1]
        s2 = self.sequences[seq2]
        euc = euclidean_dist(s1, s2)
        return euc < euclidean_threshold



class ContigPool(SequencePool):


    def process_contigs(
        self,
        score_vec: NDArray,
        ccl: NDArray,
        out_dir: str,
        lam: float,
        batch: int,
        mu: int = 400,
        node_size: int = 100,
    ) -> dict[str, NDArray]:
        """
        WRAPPER to process contigs using m0 strategy.

        :param score_vec: The score vector.
        :param node_size: The node size.
        :param ccl: Read length distribution.
        :param out_dir: The output directory.
        :param mu: Length of anchor bases.
        :param lam: The lam value.
        :param batch: The batch number.
        :return: A dictionary containing the new strategies for contigs.
        """
        logging.info("finding new strategies..")
        self._chunk_up_contigs(node_size=node_size)
        self._contigs_scores(score_vec=score_vec, node_size=node_size)
        self._process_contig_ends(node_size=node_size)
        self._contigs_benefits(ccl=ccl, mu=mu, node_size=node_size)
        t = self.find_threshold(mu=mu, lam=lam, node_size=node_size)
        # find and write new strategies
        contig_strats = self._find_contig_strategies(t=t)
        logging.info("writing new strategies")
        self._write_contig_strategies(out_dir=out_dir, contig_strats=contig_strats)
        self._write_index_file(out_dir=out_dir, batch=batch)
        return contig_strats


    def _chunk_up_contigs(self, node_size: int = 100) -> None:
        """
        Chunk up collection of contigs by giving them a chunked representation to reduce resolution

        :param node_size: The node size.
        """
        n_comp = 0
        lengths = []
        for header, seqo in self.sequences.items():
            seqo.chunk_up_coverage(n=node_size)
            n_comp += 1
            lengths.append(len(seqo.seq))
        lengths_sort = np.sort(lengths)[::-1]
        logging.info(f'num components: {n_comp}')
        logging.info(f'total comp length: {lengths_sort.sum()}')
        logging.info(f'longest components: {lengths_sort[:10]}')


    def _contigs_scores(self, score_vec: NDArray, node_size: int = 100) -> None:
        """
        Get the scores for each contig.

        :param score_vec: The score vector.
        :param node_size: The node size.
        """
        for header, seqo in self.sequences.items():
            seqo.contig_scores(score_vec=score_vec, n=node_size)


    def _contigs_benefits(self, ccl: NDArray, mu: int, node_size: int = 100) -> None:
        """
        Get the benefit for each contig.

        :param ccl: Read length distribution
        :param mu: The length of anchor bases.
        :param node_size: The node size.
        """
        for header, seqo in self.sequences.items():
            seqo.contig_benefits(mu=mu, ccl=ccl, node_size=node_size)


    def find_threshold(self, mu: float, lam: float, node_size: int = 100) -> float:
        """
        Find the acceptance threshold based on all benefit values.

        :param mu: Length of anchor bases
        :param lam: The lam value.
        :param node_size: The node size.
        :return: Acceptance threshold.
        """
        # flatten all benefit values
        benefit = np.column_stack([seqo.benefit for seqo in self.sequences.values()]).ravel()
        smu_sum = np.sum([seqo.smu_sum for seqo in self.sequences.values()])
        # find acceptance threshold
        alpha = 200 // node_size
        rho = 300 // node_size
        tc = (lam - mu - 300) // node_size

        benefit_bin, counts = Benefit.benefit_bins(benefit)
        # average benefit of strategy in the case that all fragments are rejected
        ubar0 = smu_sum
        tbar0 = alpha + rho + (mu // node_size)
        # cumsum of the benefit (bins multiplied by how many sites are in the bin)
        cs_u = np.cumsum(benefit_bin * counts) + ubar0
        cs_t = np.cumsum(tc * counts) + tbar0
        peak = cs_u / cs_t
        strat_size = np.argmax(peak) + 1
        # plt.plot(cs_u)
        # plt.plot(cs_t)
        # plt.plot(peak)
        # plt.show()
        # calculate threshold exponent and where values are geq
        try:
            threshold = float(benefit_bin[strat_size])
        except IndexError:
            threshold = float(benefit_bin[-1])
        return threshold


    def _process_contig_ends(self, node_size: int = 100) -> None:
        """
        Set contig ends to reflect their status during calculation of the strategy.

        :param node_size: The node size.
        """
        for header, seqo in self.sequences.items():
            seqo.set_contig_ends(n=node_size)


    def _find_contig_strategies(self, t: float = 0.0) -> dict[str, NDArray]:
        """
        Find the strategies for all contigs.

        :param t: The score threshold. Defaults to 0.0.
        :return: A dictionary containing the new strategies for contigs.
        """
        contig_strats = {}
        for header, seqo in self.sequences.items():
            cstrat = seqo.find_strat_m0(threshold=t)
            contig_strats[header] = cstrat
        return contig_strats


    @staticmethod
    def _write_contig_strategies(out_dir: str, contig_strats: dict[str, NDArray]) -> None:
        """
        Write the strategies for all contigs to a single file.

        :param out_dir: The output directory.
        :param contig_strats: A dictionary containing the strategies for contigs.
        """
        cpath_tmp = f'{out_dir}/masks/boss_tmp.npz'
        np.savez(cpath_tmp, **contig_strats)
        # after writing to tmpfile, rename to replace the current strat
        cpath = f'{out_dir}/masks/boss.npz'
        Path(cpath_tmp).rename(cpath)
        # Example how to load these:
        # container = np.load(f'{cpath}.npz')
        # data = {key: container[key] for key in container}


    def _write_index_file(self, out_dir: str, batch: int) -> None:
        """
        Write a new contig file to map against. Readfish generates the mmi

        :param out_dir: The output directory.
        :param batch: The batch number.
        """
        fa_path_tmp = f'{out_dir}/contigs/aeons_tmp.fa'
        # save the contigs to fasta
        with open(fa_path_tmp, 'w') as fasta:
            for sid, seqo in self.sequences.items():
                fasta.write(f'>{sid}\n')
                fasta.write(f'{seqo.seq}\n')
        # rename file after being written to
        fa_path = f'{out_dir}/contigs/aeons.fa'
        Path(fa_path_tmp).rename(fa_path)
        # copy previous contigs
        if batch % 10 == 0:
            copy(fa_path, f'{out_dir}/contigs/prev/aeons_{batch}.fa')


class Unitig:

    def __init__(self, sa_line_list: list[str], x_line: str):
        """
        Initialize Unitig object.

        :param sa_line_list: A list of lines from a GFA file corresponding to the S & A lines of a single unitig.
        :param x_line: The X line from the GFA file corresponding to the unitig.
        """
        self.name = random_id()
        self._format_sa(sa_line_list)
        self._format_x(x_line)
        self.cov = None


    def _format_x(self, x_line: str) -> None:
        """
        Format the X line of the unitig.

        :param x_line: The X line from the GFA file.
        """
        xl = x_line.split("\t")
        assert xl[0].startswith('utg')
        try:
            cl = int(xl[3])
            cr = int(xl[4])
        except IndexError:
            self.cap_l = False
            self.cap_r = False
            return
        self.cap_l = True if cl > 0 else False
        self.cap_r = True if cr > 0 else False



    def _format_sa(self, line_list: list[str]) -> None:
        """
        Format the S & A lines of the unitig.

        :param line_list: A list of lines from the GFA file corresponding to the S & A lines of the unitig.
        """
        self._format_sline(sline=line_list[0])
        self._format_atoms(alines=line_list[1:])


    def _format_sline(self, sline: str) -> None:
        """
        Format the S line of the unitig.

        :param sline: The S line from the GFA file.
        """
        sls = sline.split('\t')
        self.seq = sls[1]
        self.length = int(sls[2].split(':')[-1])
        circ = sls[0][-1]
        self.circ = True if circ == 'c' else False
        assert sls[0].startswith('utg')
        assert self.length == len(self.seq)


    def _format_atoms(self, alines: list[str]) -> None:
        """
        Format the A lines (atoms) of the unitig.

        :param alines: A list of lines
        """
        atoms = []
        for line in alines:
            assert line.startswith('A')
            atom = {}
            al = line.split('\t')
            atom['pos'] = int(al[2])
            atom['strand'] = al[3]
            if atom['strand'] == '-':
                atom['rev'] = 1
            elif atom['strand'] == '+':
                atom['rev'] = 0
            else:
                print(line)
                print(al)
                print(atom)
                logging.info("wrong strand spec of unitig")
                exit(1)
            atom['name'] = al[4]
            atoms.append(atom)
        # loop a second time to get the ends
        cpos = 0
        for i in range(len(alines) - 1):
            line = alines[i + 1]
            al = line.split('\t')
            pos = int(al[2])
            to_add = pos - cpos
            atoms[i]['n'] = to_add
            cpos = pos
        # mark last atom
        atoms[-1]['n'] = -1
        self.atoms = atoms
        self.atom_headers = [a['name'] for a in atoms]


    def to_seqo(self, seqpool: SequencePool) -> Sequence:
        """
        Transform unitig into a sequence object. Needs to be done after merging coverage.

        :param seqpool: The sequence pool of the source sequence that made the unitig.
        :return: A sequence object representing the unitig.
        """
        # check that coverage merging has been done
        assert self.cov is not None
        # grab atoms of atoms
        merged_atoms = seqpool.get_atoms(headers=self.atom_headers)
        merged_components = seqpool.get_components(headers=self.atom_headers)
        seqo = Sequence(header=self.name, seq=self.seq, cov=self.cov,
                        merged_components=merged_components, merged_atoms=merged_atoms,
                        cap_l=self.cap_l, cap_r=self.cap_r)
        return seqo



class UnitigPool:

    def __init__(self, unitigs: list[Unitig]):
        """
        Initialize a UnitigPool object.

        :param unitigs: A list of Unitig objects.
        """
        self.unitigs = unitigs


    def get_unitig_coverage_arrays(self, seqpool: SequencePool) -> None:
        """
        Perform coverage array merging for each unitig in the pool.

        :param seqpool: The sequence pool object.
        """
        for u in self.unitigs:
            cm = CoverageMerger(u, seqpool.sequences)
            cov_arr = cm.cov_arr
            u.cov = cov_arr


    def unitigs2seqpool(self, seqpool: SequencePool, min_seq_len: int) -> tuple[SequencePool, set[str]]:
        """
        Transform all unitigs in the pool to sequence objects and obtain the read IDs to remove.

        :param seqpool: The sequence pool object.
        :param min_seq_len: The minimum sequence length for the new pool.
        :return: A tuple containing the new sequence pool and the set of used read IDs.
        """
        seqos = {}
        used_sids = set()
        for u in self.unitigs:
            unitig_seqo = u.to_seqo(seqpool)
            seqos[u.name] = unitig_seqo
            used_sids.update(u.atom_headers)
        # construct a new pool
        new_pool = SequencePool(sequences=seqos, min_len=min_seq_len)
        return new_pool, used_sids



class CoverageMerger:


    def __init__(self, unitig: Unitig, seqpool: dict[str, Sequence]) -> None:
        """
        Initializes a CoverageMerger object.

        :param unitig: Unitig object.
        :param seqpool: Dictionary of sequences used to create the unitigs
        :raises AssertionError: If the length of the unitig is not equal to the shape of the merged coverage array.

        """
        self.unitig = unitig
        cov_arr = self._create_merged_arr(atoms=unitig.atoms, seqpool=seqpool)
        assert unitig.length == cov_arr.shape[0]
        self.cov_arr = cov_arr


    def _create_merged_arr(self, atoms: list[dict[str, Any]], seqpool: dict[str, Sequence]) -> NDArray:
        """
        Creates merged coverage array for the unitig.

        :param atoms: List of atom dictionaries.
        :param seqpool: Dictionary of sequences used to create the unitigs
        :return: The merged coverage array.
        """
        arr_parts = []
        cpos = 0
        for a in atoms:
            assert a['pos'] >= cpos
            name = a['name']
            atom_arr = seqpool[name].cov.copy()
            atom_arr = atom_arr[::-1] if a['rev'] else atom_arr
            if not a['n'] == -1:
                atom_arr = atom_arr[: a['n']]
            else:
                # add last atom (goes until end of atom)
                if not self.unitig.circ:
                    atom_arr = atom_arr
                else:
                    diff = self.unitig.length - a['pos']
                    atom_arr = atom_arr[:diff]
            arr_parts.append(atom_arr)
            cpos = a['pos']
        return np.concatenate(arr_parts)




class MultilineContainments:


    def __init__(self, records: list[PafLine]) -> None:
        """
        Initializes a MultilineContainments object.

        :param records: A list of PafLine objects.
        """
        self.records = records
        # collect all partners with multiple mappings
        multidict = self._fill_multidict()
        # get a dict of (contained, container): [recs]
        containments = self._get_multiline_containments(multidict=multidict)
        self.containments = containments


    def _fill_multidict(self) -> dict[str, list[PafLine]]:
        """
        Fills a multidict with records that have multiple mappings.

        :return: The multidict of (key, records) pairs.
        """
        multidict = defaultdict(list)
        for rec in self.records:
            multidict[rec.keygen()].append(rec)
        multidict = {k: recs for k, recs in multidict.items() if len(recs) > 1}
        return multidict


    def _get_multiline_containments(self, multidict: dict[str, list[PafLine]]) -> dict[Edge, PafLine]:
        """
        Gets the multiline containments from the multidict using multiple internal mappings.

        :param multidict: The multidict of (key, records) pairs.
        :return: Dictionary of (contained, container) tuples mapped to PafLine objects.
        """
        containments = {}
        for k, recs in multidict.items():
            cont = self.multiline_containment(recs)
            if cont:
                containments.update(cont)
        return containments


    @staticmethod
    def multiline_containment(records: list[PafLine], n: int = 100) -> dict[Edge, PafLine]:
        """
        Checks for multiline containment from multiple internal match mappings.

        :param records: A list of PafLine objects.
        :param n: Node size, defaults to 100.
        :return: The dictionary of (contained, container) tuples mapped to PafLine objects, or empty dict
        """
        qlen = records[0].qlen // n
        tlen = records[0].tlen // n
        qarr = np.zeros(shape=qlen , dtype="bool")
        tarr = np.zeros(shape=tlen , dtype="bool")

        if len(records) > 10:
            return {}

        for r in records:
            qarr[r.qstart // n: r.qend // n] = 1
            tarr[r.tstart // n: r.tend // n] = 1

        # if more than 0.9 are covered by mappings
        if sum(qarr) > qlen * 0.9:
            # if overhang is smaller than 0.15 of len
            q_low, q_high = np.nonzero(qarr)[0][[0, -1]]
            if (q_high - q_low) > qlen * 0.85:
                # return containment tuple
                t_low, t_high = np.nonzero(tarr)[0][[0, -1]]
                ctd_low, ctd_high, ctr_low, ctr_high = q_low, q_high, t_low, t_high
                cont = MultilineContainments.generate_paf_cont(
                    records, 'q', 't', ctd_low, ctd_high, ctr_low, ctr_high, n
                )
                return cont

        if sum(tarr) > tlen * 0.9:
            t_low, t_high = np.nonzero(tarr)[0][[0, -1]]
            if (t_high - t_low) > tlen * 0.85:
                q_low, q_high = np.nonzero(qarr)[0][[0, -1]]
                ctd_low, ctd_high, ctr_low, ctr_high = t_low, t_high, q_low, q_high
                cont = MultilineContainments.generate_paf_cont(
                    records, 't', 'q', ctd_low, ctd_high, ctr_low, ctr_high, n
                )
                return cont

        # if neither q nor t are contained
        return {}


    @staticmethod
    def generate_paf_cont(
        records: list[PafLine],
        ctd: str,
        ctr: str,
        ctd_low: int,
        ctd_high: int,
        ctr_low: int,
        ctr_high: int,
        n: int
    ) -> dict[Edge, PafLine]:
        """
        Generates a PafLine object that describes the containment,
         given a list of records that make up one multicontainment

        :param records: A list of PafLine objects.
        :param ctd: The contained sequence name.
        :param ctr: The container sequence name.
        :param ctd_low: The lower bound of the contained sequence.
        :param ctd_high: The upper bound of the contained sequence.
        :param ctr_low: The lower bound of the container sequence.
        :param ctr_high: The upper bound of the container sequence.
        :param n: Node size.
        :return: The dictionary of (contained, container) tuples mapped to PafLine objects.
        """
        ctd_name = getattr(records[0], f'{ctd}name')
        ctr_name = getattr(records[0], f'{ctr}name')
        ctd_len = getattr(records[0], f'{ctd}len')
        ctr_len = getattr(records[0], f'{ctr}len')
        # use maximum span
        ctd_span = ctd_high - ctd_low
        ctr_span = ctr_high - ctr_low
        # except if span on container longer than 2.2x of contained
        if ctr_span > 2.2 * ctd_span:
            r = 0
            maplen = 0
            for i in range(len(records)):
                # then use the longest alignment between the two
                if records[i].maplen > maplen:
                    maplen = records[i].maplen
                    r = i
            ctr_low = getattr(records[r], f'{ctr}start') // n
            ctr_high = getattr(records[r], f'{ctr}end') // n
        # generate paf entry
        paf = f'{ctd_name}\t{ctd_len}\t{ctd_low * n}\t{ctd_high * n}\t+' \
              f'\t{ctr_name}\t{ctr_len}\t{ctr_low * n}\t{ctr_high * n}\t0\t0\t0'
        rec = PafLine(paf)
        # mark as first contained
        rec.c = 2
        return {(ctd_name, ctr_name): rec}




class Benefit:

    @staticmethod
    def init_scoring_vec(lowcov: float) -> NDArray:
        """
        Initialize scoring vector based on target coverage.

        :param lowcov: The target coverage value.
        :return: The scoring vector.
        """
        x = np.arange(101)
        # a = lowcov * 5
        # score_vec = -gamma.cdf(x, a=a, scale=0.2) + 1
        score_vec = 1 / (np.exp(x - lowcov) + 1)
        return score_vec


    @staticmethod
    def score_array(score_vec: NDArray, cov_arr: NDArray, node_size: int = 100) -> NDArray:
        """
        Calculate scores based on the scoring vector and coverage array.

        :param score_vec: scoring vector.
        :param cov_arr: coverage array.
        :param node_size: node size.
        :return: The calculated scores.
        """
        # grab scores using multi-indexing
        carr = cov_arr // node_size  # apply resolution reduction
        carr_int = carr.astype("int")
        scores = score_vec[carr_int]
        return scores


    @staticmethod
    def calc_fragment_benefit(
        scores: NDArray,
        mu: int,
        approx_ccl: NDArray,
        e1: bool,
        e2: bool,
        node_size: int = 100,
    ) -> tuple[NDArray, float]:
        """
        Calculate the benefit of a fragment based on scores, mu, node_size, approx_ccl, e1, and e2.

        :param scores: Fragment's position-wise scores.
        :param mu: Length of anchor bases.
        :param node_size: node size.
        :param approx_ccl: Approx of read length distribution.
        :param e1: Left-end marker.
        :param e2: Right-end marker.
        :return: The calculated benefit and smu_sum.
        """
        # expand score to account for contig ends
        mu_ds = mu // node_size
        ccl_ds = approx_ccl // node_size
        ccl_max = int(ccl_ds[-1])
        sx = Benefit._expand_scores(scores, e1, e2, ccl_max)
        smu = Benefit._calc_smu_moving(score=sx, mu_ds=mu_ds)
        benefit = Benefit._calc_benefit_moving(score=sx, ccl_ds=ccl_ds)
        smu_sum = float(np.sum(smu))
        b = benefit - smu
        b[b < 0] = 0
        b = b[:, ccl_max: -ccl_max]
        assert b.shape[1] == scores.shape[0]
        return b, smu_sum


    @staticmethod
    def _expand_scores(scores: NDArray, e1: bool, e2: bool, ccl_max: int) -> NDArray:
        """
        Expand scores to account for contig ends.

        :param scores: Fragment scores.
        :param e1: Left-end marker.
        :param e2: Right-end marker.
        :param ccl_max: Max of approx read length dist.
        :return: The expanded scores.
        """
        scoresx = np.zeros(shape=scores.shape[0] + (ccl_max * 2), dtype="float64")
        scoresx[ccl_max: -ccl_max] = scores
        scoresx[0: ccl_max] = 1 if e1 else 0
        scoresx[-ccl_max: -1] = 1 if e2 else 0
        return scoresx


    @staticmethod
    def _calc_smu_moving(score: NDArray, mu_ds: int) -> NDArray:
        """
        Calculate smu moving based on score and down-sampled mu.

        :param score: The score.
        :param mu_ds: The mu_ds.
        :return: The calculated smu.
        """
        smu_fwd = bn.move_sum(score, window=mu_ds, min_count=1)
        smu_rev = bn.move_sum(score[::-1], window=mu_ds, min_count=1)
        smu = np.stack((smu_fwd, smu_rev))
        return smu


    @staticmethod
    def _calc_benefit_moving(score: NDArray, ccl_ds: NDArray) -> NDArray:
        """
        Calculate benefit moving based on score and ccl_ds.

        :param score: Fragment scores.
        :param ccl_ds: Down-sampled read length dist array.
        :return: The calculated benefit.
        """
        score_rev = score[::-1]
        benefit = np.zeros(shape=(2, score.shape[0]), dtype="float64")
        perc = np.arange(0.1, 1.1, 0.1)[::-1]
        assert perc.shape == ccl_ds.shape
        for i in range(ccl_ds.shape[0]):
            ben_fwd = bn.move_sum(score, window=int(ccl_ds[i]), min_count=1)[ccl_ds[i]: -1]
            ben_rev = bn.move_sum(score_rev, window=int(ccl_ds[i]), min_count=1)[ccl_ds[i]: -1]
            benefit[0, 0: -ccl_ds[i] - 1] += ben_fwd * perc[i]
            benefit[1, ccl_ds[i]: -1] += ben_rev[::-1] * perc[i]
        return benefit


    @staticmethod
    def benefit_bins(benefit: NDArray) -> tuple[NDArray, NDArray]:
        """
        Group benefit into bins of similar values using binary exponent. Used to find acceptance threshold

        :param benefit: positional benefit array.
        :return: The benefit bins and counts.
        """
        benefit_nz_ind = np.nonzero(benefit)
        benefit_flat_nz = benefit[benefit_nz_ind]
        # to make binary exponents work, normalise benefit values
        normaliser = np.max(benefit_flat_nz)
        benefit_flat_norm = benefit_flat_nz / normaliser
        mantissa, benefit_exponents = np.frexp(benefit_flat_norm)
        # count how often each exponent is present
        # absolute value because counting positive integers is quicker
        benefit_exponents_pos = np.abs(benefit_exponents)
        # multi-thread counting of exponents
        exponent_arrays = np.array_split(benefit_exponents_pos, 12)
        with TPexe(max_workers=12) as executor:
            exponent_counts = executor.map(np.bincount, exponent_arrays)
        exponent_counts = list(exponent_counts)
        # aggregate results from threads
        # target array needs to have the largest shape of the thread results
        max_exp = np.max([e.shape[0] for e in exponent_counts])
        bincounts = np.zeros(shape=max_exp, dtype='int')
        # sum up results from individual threads
        for exp in exponent_counts:
            bincounts[0:exp.shape[0]] += exp
        # filter empty bins
        exponents_unique = np.nonzero(bincounts)[0]
        # counts of the existing benefit exponents
        counts = bincounts[exponents_unique]
        # use exponents to rebuild benefit values
        benefit_bin = np.power(2.0, -exponents_unique) * normaliser
        return benefit_bin, counts




