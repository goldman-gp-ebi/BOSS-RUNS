import logging
from concurrent.futures import ThreadPoolExecutor as TPexe

# non-std lib
import mappy


class Mapper:

    def __init__(self, mu):
        self.mu = mu


    def init_mappy(self, ref):
        """
        Load the index of the reference genome that we map against

        Parameters
        ----------
        ref: str
            Path to the reference fasta
        Returns
        -------
        """
        self.aligner = mappy.Aligner(fn_idx_in=ref, preset="map-ont")  # , extra_flags=0x400)


    def mappy_batch(self, sequences, truncate=False, workers=24):
        """
        Wrapper function that maps a full batch of reads.
        Also used for simulations, where reads are truncated to mu bases

        Parameters
        ----------
        sequences: dict
            dict of read_id: sequence items loaded as instance of current batch
        truncate: bool
            simulations: truncate the input reads. not needed for live application
        workers: int
            number of threads

        Returns
        -------
        alignments: str
            PAF formatted mapping hits of all mapped input reads
        """
        # container to hold the hits from all reads
        batch_alignments = []
        # for in silico experiments, truncate the reads
        if truncate:
            sequences = {read_id: seq[:self.mu] for read_id, seq in sequences.items()}

        unmapped_ids = 0
        mapped_ids = 0
        # loop over all sequences and map them one by one
        with TPexe(max_workers=workers) as executor:
            results = executor.map(self.map_sequence, sequences.items())

        for result in results:
            res = '\n'.join(result)
            # prevent appending an empty list if the read was not mapped
            if len(res) > 0:
                batch_alignments.append(res)
                mapped_ids += 1
            else:
                unmapped_ids += 1

        # write out unmapped ids for diagnosis
        # self.unmapped_ids = unmapped_ids
        # transform to a single string
        alignments = '\n'.join(batch_alignments)
        # counting the number of mapped and unmapped fragments - used for Fhat
        if not truncate:
            logging.info(f"Number of mapped reads: {mapped_ids}, unmapped reads: {unmapped_ids}")
        else:
            logging.info(f"Number of mapped reads: {mapped_ids}, unmapped reads: {unmapped_ids} (trunc)")
        return alignments


    def map_sequence(self, sequence):
        """
        Fast mapper that takes a sequence and returns the mapped sequence.
        Parameters
        ----------
        sequence: dict item
            key:value pair of read_id and sequence
        Returns
        -------
        list of str
            List of Paf formatted mapping hits
        """
        results = []
        read_id, seq = sequence

        thr_buf = mappy.ThreadBuffer()
        # For each alignment to be mapped against, returns a PAF format line
        for hit in self.aligner.map(seq, buf=thr_buf):
            if hit.is_primary:
                results.append(f"{read_id}\t{len(seq)}\t{hit}")
        return results


