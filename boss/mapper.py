from io import StringIO
from concurrent.futures import ThreadPoolExecutor as TPE
from pathlib import Path
import logging

import mappy

from .paf import Paf, paf_dict_type



class Indexer:

    def __init__(self, fasta: str, mmi: str, t: int = 4):
        """
        Initialize simple indexing wrapper around mm2

        :param fasta: The path to the input FASTA file.
        :param mmi: The path to the output .mmi index file.
        :param t: The number of threads to use for indexing.
        """
        self.aligner = mappy.Aligner(fn_idx_in=fasta, fn_idx_out=mmi, preset="map-ont", n_threads=t)



class Mapper:

    def __init__(self, ref: str, mu: int = 400, workers: int = 4, default: bool = True):
        """
        Initialize a Mapper object; wrapper for minimap2's mappy implementation
        For the default case of mapping against some linear references

        :param ref: The path to the reference FASTA file.
        :param mu: Length of anchor bases for simulations, defaults to 400.
        :param workers: The number of worker threads to use for mapping
        :param default: Whether to use the default mappy.Aligner configuration, defaults to True.
        """
        self.mu = mu
        self.workers = workers
        # check that the given reference exists
        if not Path(ref).is_file():
            raise FileNotFoundError("Given reference file does not exist")

        if default:
            self.aligner = mappy.Aligner(fn_idx_in=ref, preset="map-ont")
        else:
            self.aligner = mappy.Aligner(fn_idx_in=ref, fn_idx_out=f'{ref}.mmi', preset="map-ont",
                                         k=13, w=5, min_cnt=2, min_chain_score=20)



    def map_sequences(self, sequences: dict[str, str], trunc: bool = False) -> paf_dict_type:
        """
        Map sequences to the reference and return the mapping results.

        :param sequences: A dictionary of read_id: sequence items.
        :param trunc: Switch to perform mapping of mu-sized fragments for aeons simulations
        :return: A dictionary of read_id: list(PafLine) objects.
        """
        # truncation option for aeons simulations
        if trunc:
            sequences = {rid: seq[: self.mu] for rid, seq in sequences.items()}
        paf_raw = self._mappy_batch(sequences=sequences)
        paf_dict = Paf.parse_PAF(StringIO(paf_raw), min_len=int(self.mu / 2))
        return paf_dict



    def _mappy_batch(self, sequences: dict[str, str], out: str = None, log: bool = True) -> str:
        """
        Map a full batch of reads and return the PAF formatted mapping hits.

        :param sequences: A dictionary of read_id: sequence items.
        :param out: The name of the file to output mappings, defaults to None.
        :param log: Whether to log the number of mapped and unmapped queries, defaults to True.
        :return: The PAF formatted mapping hits of all mapped input reads.
        """
        # container to hold the hits from all reads
        batch_alignments = []
        unmapped_count = 0
        mapped_count = 0
        # loop over all sequences and map them one by one
        with TPE(max_workers=self.workers) as executor:
            results = executor.map(self._map_query, sequences.items())

        for result in results:
            res = '\n'.join(result)
            # prevent appending an empty list if the read was not mapped
            if len(res) > 0:
                batch_alignments.append(res)
                mapped_count += 1
            else:
                unmapped_count += 1

        # transform to a single string
        alignments = '\n'.join(batch_alignments)

        # tmp write to file too
        if out:
            with open(out, 'w') as lm_out:
                lm_out.write(alignments)

        # counting the number of mapped and unmapped fragments
        if log:
            logging.info(f"MAPPY: mapped queries: {mapped_count}, unmapped queries: {unmapped_count} ")
        self.mapped_count = mapped_count
        self.unmapped_count = unmapped_count
        return alignments



    def _map_query(self, query: tuple[str, str]) -> list[str]:
        """
        Map a query and return the mapped query as a list of PAF formatted mapping hits.

        :param query: A key-value pair of read_id and query.
        :return: A list of PAF formatted mapping hits.
        """
        results = []
        read_id, seq = query

        thr_buf = mappy.ThreadBuffer()
        # For each alignment to be mapped against, returns a PAF format line
        for hit in self.aligner.map(seq, buf=thr_buf):
            # if hit.is_primary:
            results.append(f"{read_id}\t{len(seq)}\t{hit}")
        return results




