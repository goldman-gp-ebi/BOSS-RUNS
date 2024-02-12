from io import StringIO
from collections import defaultdict
from types import SimpleNamespace
from pathlib import Path

import numpy as np





class PafLine:
    '''
    @DynamicAttrs
    parse a single alignment from a PAF into a flexible container
    '''

    def __init__(self, line: str, tags: bool = True):
        """
        Parse a line from a PAF file

        :param line: string representation of a line in a PAF file
        :param tags: boolean indicator whether to parse tags
        """
        self.line = line
        fields = ['qname', 'qlen', 'qstart', 'qend',
                  'strand', 'tname', 'tlen', 'tstart', 'tend',
                  'num_matches', 'alignment_block_length',
                  'mapq']
        core = 12
        record = line.strip().split("\t")

        f = PafLine.format_records(record[:core])
        for i in range(core):
            setattr(self, fields[i], f[i])

        # make sure query and target name are strings
        self.qname = str(self.qname)
        self.tname = str(self.tname)

        self.rev = 0 if self.strand == '+' else 1
        # parse the tags only if needed
        if tags:
            tags_parsed = PafLine.parse_tags(record[core:])
            self.align_score = int(tags_parsed.get("AS", 0))
            self.cigar = tags_parsed.get("cg", None)
            self.s1 = tags_parsed.get("s1", 0)
            prim = tags_parsed.get("tp", None)
            self.primary = 1 if prim == 'P' else 0
            # self.dv = tags_parsed.get('dv', 0)

        # markers for trimming
        self.qprox = False
        self.tprox = False
        # dummy inits
        self.maplen = None
        self.min_length_pair = None


    @staticmethod
    def format_records(record: list) -> list:
        """
        Helper function to make fields of a PafLine the right type

        :param record: Split string of PAFline into list
        :return: Same list but with types converted to int
        """
        return [PafLine.conv_type(x, int) for x in record]


    @staticmethod
    def parse_tags(tags: list) -> dict:
        """
        Parse tags of a PAFline into a dictionary

        :param tags: List of SAM style tags
        :return: Dict of SAM style tags
        """
        c = {"i": int, "A": str, "f": float, "Z": str}
        return {
            key: PafLine.conv_type(val, c[tag])
            for key, tag, val in (x.split(":") for x in tags)
        }


    @staticmethod
    def conv_type(s: str, func: callable):
        """
        Generic converter, to change strings to other types

        :param s: Input string to convert to a different type
        :param func: Target type of input string
        :return: Either converted or original type
        """
        try:
            return func(s)
        except ValueError:
            return s


    def filter(self, filters: SimpleNamespace) -> bool:
        """
        Check if an alignment needs to be filtered out

        :param filters: Filters stored as arguments
        :return: Indicator whether line is filtered
        """
        # like classify, pack all conditions in here
        if self._self_aligned():
            return True
        if self.map_length() < filters.min_map_len:
            return True
        if self.s1 < filters.min_s1:
            return True
        if self.min_length_in_pair() < filters.min_seq_len:
            return True
        # if none of the filters triggered
        return False


    def min_length_in_pair(self) -> float:
        """
        Get the shorter length of the two aligned sequences

        :return: Length of shorter sequence in pair
        """
        if not self.min_length_pair:
            self.min_length_pair = min(self.qlen, self.tlen)
        return self.min_length_pair


    def overhang(self) -> float:
        """
        Get the sum of the smallest overhangs on both sequences
        Used in classification of mapping

        :return: Sequence overhang in alignment
        """
        if not self.rev:
            overhang = min(self.qstart, self.tstart) +\
                       min(self.qlen - self.qend, self.tlen - self.tend)
        else:
            overhang = min(self.qstart, self.tlen - self.tend) +\
                       min(self.tstart, self.qlen - self.qend)
        return overhang


    def map_length(self) -> float:
        """
        Get the shorter alignment length. I.e. either the span on query or target

        :return: Mapping length of the paf record
        """
        if not self.maplen:
            self.maplen = min(self.qend - self.qstart, self.tend - self.tstart)
        return self.maplen


    def classify(self) -> int:
        """
        Classify the alignment according to miniasm algorithm 5
        Records that are being classified have already been filtered

        :return: Indicator of alignment type
        """
        c = -1

        if self._internal_match():
            c = 1
        elif self._first_contained():
            c = 2
        elif self._second_contained():
            c = 3
        elif self._first_contained_fallback():
            c = 2
        elif self._second_contained_fallback():
            c = 3
        else:
            pass

        # if still unclassified -> overlap
        if c < 0:
            c, qside, tside = self._overlap_orientation()
            self.qside = qside
            self.tside = tside

        # second chance for internal matches
        # consider them as ovl under special circumstances
        if c == 1:
            if self._first_contained_mostly():
                c = 2
            elif self._second_contained_mostly():
                c = 3
            elif self.internal_match_is_overlap():
                c = 6  # class 6: needs trimming
            else:
                pass

        return c



    def _internal_match(self) -> bool:
        """
        Check if the alignment is an internal match.
        We consider a maximum of 15 percent of the map length as limit

        :return: indicator if internal match
        """
        if self.overhang() > (self.map_length() * 0.15):
            return True
        else:
            return False



    def _first_contained(self) -> bool:
        """
        First contained means the query is contained in the target.
        This can be different in either forward or reverse relative orientation

        :return: indicator for query contained in target
        """
        # FORWARD
        if not self.rev:
            if (self.qstart <= self.tstart) and ((self.qlen - self.qend) < (self.tlen - self.tend)):
                return True
            else:
                return False
        # REVERSE
        else:
            if (self.qstart <= (self.tlen - self.tend)) and ((self.qlen - self.qend) < self.tstart):
                return True
            else:
                return False



    def _second_contained(self) -> bool:
        """
        Second contained means the target is contained in the query.
        This can be different in either forward or reverse relative orientation

        :return: indicator for target contained in query
        """
        # FORWARD
        if not self.rev:
            if (self.qstart >= self.tstart) and ((self.qlen - self.qend) > (self.tlen - self.tend)):
                return True
            else:
                return False
        # REVERSE
        else:
            if (self.qstart >= (self.tlen - self.tend)) and ((self.qlen - self.qend) > self.tstart):
                return True
            else:
                return False



    def _first_contained_fallback(self) -> bool:
        """
        Handling edge cases of containment: if more than 90 percent of the sequence
        is covered my a mapping, then we assume contained

        :return: Indicator for query contained in target
        """
        qcov = self.qend - self.qstart
        if (qcov / self.qlen) >= 0.90:
            return True
        else:
            return False



    def _second_contained_fallback(self) -> bool:
        """
        Handling edge cases of containment: if more than 90 percent of the sequence
        is covered my a mapping, then we assume contained

        :return: Indicator for target contained in query
        """
        tcov = self.tend - self.tstart
        if (tcov / self.tlen) >= 0.90:
            return True
        else:
            return False



    def _first_contained_mostly(self) -> bool:
        """
        This is to consider internal matches of long sequences as containments instead.
        This is done so that we can count the coverage information of such events

        :return: indicator for mostly-contained internal matches
        """
        qcov = self.qend - self.qstart
        if (qcov / self.qlen) >= 0.50 and self.qlen > 20000:
            return True
        else:
            return False



    def _second_contained_mostly(self) -> bool:
        """
        This is to consider internal matches of long sequences as containments instead.
        This is done so that we can count the coverage information of such events

        :return: indicator for mostly-contained internal matches
        """
        tcov = self.tend - self.tstart
        if (tcov / self.tlen) >= 0.50 and self.qlen > 20000:
            return True
        else:
            return False



    def _overlap_orientation(self) -> tuple[int, str, str]:
        """
        Identify details of overlaps, i.e. which overlaps which and on which end

        :return: Indicator for overlap source and their respective ends
        """
        if not self.rev:
            if self.qstart > self.tstart:
                # 'A' overlaps 'B'
                # a + b +
                return 4, 'R', 'L'
            else:
                # B overlaps A
                # b + a +
                return 5, 'L', 'R'
        elif self.qstart > (self.qlen - self.qend):
            if self.qstart > (self.tlen - self.tend):
                # 'A' overlaps 'B'
                # a + b -
                return 4, 'R', 'R'    # should this be LR?
            else:
                # B overlaps A
                # b + a -
                return 5, 'R', 'R'
        elif (self.qlen - self.qstart) > self.tend:
            # 'A' overlaps 'B'
            # a - b +
            return 4, 'L', 'L'   # should this be RL?
        else:
            # B overlaps A
            # b - a +
            return 5, 'L', 'L'



    def _self_aligned(self) -> bool:
        """
        Check if query is also target

        :return: Indicator if alignment to itself
        """
        if self.qname == self.tname:
            return True
        else:
            return False



    # def _is_merged(self) -> bool:
    #     """
    #     Check for internal match turned OVL restrictions
    #     This might be deprecated
    #
    #     :return: Indicator for already merged seqs
    #     """
    #     if '*' in self.qname and '*' in self.tname:
    #         return True
    #     else:
    #         return False



    @staticmethod
    def _is_prox(start: int, end: int, length: int, lim=1000.0) -> bool:
        """
        check if a record has a mapped region close to one of its ends
        masm definition would be lim = 1000
        a fixed basepair limit does not work for us, we need to trim off more
        i.e. if limit is given as percentage, calc a limit

        :param start: start position on seq
        :param end: end position on seq
        :param length: sequence length
        :param lim: limit of allowed overhang
        :return: indicator whether mapping in proximity to sequence end
        """
        if lim < 1:
            limit = lim * length
        else:
            limit = lim
        overhang = min(start, length - end)
        ovl = True if overhang < limit else False
        return ovl



    def _im_ovl_restrictions(self) -> bool:
        """
        check additional restrictions to count internal match as OVL

        :return: indicator for internal match to be overlap
        """
        maplen = self.map_length()
        if self.qlen > 15000:
            if self.tlen > 15000:
                if maplen > 5000:
                    return True
        return False



    def internal_match_is_overlap(self) -> bool:
        """
        due to overlapping untrimmed reads, we reconsider internal matches that might be overlaps
        if one record has a true dovetail, check if the other has a relaxed dovetail

        :return: indicator for internal match to be overlap instead
        """
        lim = 0.15
        if self._is_prox(start=self.qstart, end=self.qend, length=self.qlen):
            self.qprox = True  # mark which side the true prox is for trimming
            if self._is_prox(start=self.tstart, end=self.tend, length=self.tlen, lim=lim):
                if self._im_ovl_restrictions():
                    # check a few more restrictions
                    # relaxed dovetail
                    return True
        elif self._is_prox(start=self.tstart, end=self.tend, length=self.tlen):
            self.tprox = True  # marker for trimming
            if self._is_prox(start=self.qstart, end=self.qend, length=self.qlen, lim=lim):
                if self._im_ovl_restrictions():
                    # relaxed dovetail
                    return True
        else:
            # neither sequence has a mapping close to the end
            # i.e. true internal match
            return False
        return False



    @staticmethod
    def _find_coords(start: int, end: int, length: int) -> tuple[int, int]:
        """
        find coordinates to trim off of reads

        :param start: start of mapping on sequence to be trimmed
        :param end: end of mapping on sequence to be trimmed
        :param length: length of sequence to be trimmed
        :return: start and end coords to trim
        """
        # which side is closer to the end?
        min_overhang_idx = np.argmin([start, length - end])
        if min_overhang_idx == 0:
            # trim_start, trim_stop = start, None
            trim_start, trim_stop = 0, start
        else:
            # trim_start, trim_stop = 0, end
            trim_start, trim_stop = end, None
        return trim_start, trim_stop



    def find_trim_coords(self) -> tuple[str, int, int, str]:
        """
        if this alignment has been identified to be useful when trimmed,
        find which of the sequences we want to trim and which coordinates

        :return: Tuple of sequence name, start, end, and partner sequence
        """
        if self.qprox:
            # if q is a real prox, trim t
            sid = self.tname
            trim_start, trim_stop = self._find_coords(start=self.tstart, end=self.tend, length=self.tlen)
            other = self.qname
            other_len = self.qlen
            orig_len = self.tlen
        else:
            sid = self.qname
            trim_start, trim_stop = self._find_coords(start=self.qstart, end=self.qend, length=self.qlen)
            other = self.tname
            other_len = self.tlen
            orig_len = self.qlen

        # check that a potential trim would actually be longer than the original
        if trim_stop is None:
            stop = orig_len
        else:
            stop = trim_stop
        trimmed_bit = stop - trim_start
        # potential length of a merged sequence: original length minus trimmed bits
        # plus the length of the newly added sequence, minus the overlap length
        new_len = orig_len - trimmed_bit + other_len - self.alignment_block_length
        if new_len < orig_len:
            sid = '0'

        return sid, trim_start, trim_stop, other



    def grab_increment_coords(self) -> tuple[int, int, int, int, int, int]:
        """
        get the coordinates of a containment to increment the coverage counts

        :return: Tuple of coordinates to increment
        """
        if self.c == 2:
            ostart = self.tstart
            oend = self.tend
            cstart = self.qstart
            cend = self.qend
        elif self.c == 3:
            ostart = self.qstart
            oend = self.qend
            cstart = self.tstart
            cend = self.tend
        else:
            raise ValueError("This method should only be called on contained reads")

        olen = oend - ostart
        clen = cend - cstart
        return ostart, oend, olen, cstart, cend, clen



    def keygen(self) -> str:
        """
        New alphabetically sorted header for multiline containments

        :return: New header
        """
        if self.qname < self.tname:
            return f'{self.qname}-{self.tname}'
        else:
            return f'{self.tname}-{self.qname}'



    # for debugging
    # def plot(self, save=None):
    #     # visualize alignments as blocks
    #     import plotnine as p9
    #     import pandas as pd
    #     # this puts together the coordinates for the polygon
    #     # 5 points, starting at the lower left end
    #     if self.strand == '+':
    #         cols = ['qstart', 'qend', 'tend', 'tstart', 'qstart']
    #     else:
    #         cols = ["qend", "qstart", "tend", "tstart", "qend"]
    #
    #     pos = [getattr(self, c) for c in cols]
    #     # the query is always on top
    #     seqn = [2, 2, 1, 1, 2]
    #     cdat = pd.DataFrame({'pos': pos, 'seqn': seqn})
    #     # coordinates of the read rectangles
    #     seqlens = pd.DataFrame({'seq': [self.qname[:10], self.tname[:10]],
    #                             'start': [0, 0],
    #                             'end': [self.qlen, self.tlen],
    #                             'bottoms': [2.05, 0.8],
    #                             'tops': [2.2, 0.95]})
    #
    #     xpos = 100
    #     qpos = 1.9
    #     tpos = 1.1
    #     p = (p9.ggplot() +
    #          p9.geom_polygon(data=cdat,
    #                          mapping=p9.aes(x="pos", y="seqn"),
    #                          fill="grey", colour="black") +
    #          p9.geom_rect(data=seqlens,
    #                       mapping=p9.aes(xmin="start", xmax="end",
    #                                      ymin="bottoms", ymax="tops"),
    #                       colour="black", fill=None) +
    #          p9.annotate(x=xpos, label=self.qname[:10], y=2.12, geom='text', ha="left", va="center") +
    #          p9.annotate(x=xpos, label=self.tname[:10], y=0.87, geom='text', ha="left", va="center") +
    #          p9.annotate(x=self.qstart, label=self.qstart, y=qpos, color="darkred",
    #                      geom='text', ha="left", va="center") +
    #          p9.annotate(x=self.qend, label=self.qend, y=qpos, color="darkred",
    #                      geom='text', ha="left", va="center") +
    #          p9.annotate(x=self.tstart, label=self.tstart, y=tpos, color="darkblue",
    #                      geom='text', ha="left", va="center") +
    #          p9.annotate(x=self.tend, label=self.tend, y=tpos, color="darkblue",
    #                      geom='text', ha="left", va="center") +
    #          p9.annotate(x=xpos, label=self.strand, y=1.4, color="darkgoldenrod",
    #                      geom='text', ha="left", va="center", size=30) +
    #          p9.ylab("") +
    #          p9.xlab("") +
    #          p9.theme_minimal() +
    #          p9.theme(axis_text_y=p9.element_blank(),
    #                   axis_ticks_major_y=p9.element_blank(),
    #                   plot_background=p9.element_rect(fill="white", color="white")))
    #     if save:
    #         p.save(save)
    #     else:
    #         return p




class Paf:

    def __init__(self):
        pass


    @staticmethod
    def parse_PAF(paf_file: str | StringIO, min_len: int = 1) -> dict:
        """
        Parse the contents of a PAF file into a dictionary of records

        :param paf_file: Can be either a string or a StringIO object
        :param min_len: minimum alignment length to consider an entry
        :return: Dict of parsed PAF file
        """
        if isinstance(paf_file, str) and Path(paf_file).is_file():
            with open(paf_file, 'r') as paff:
                paf_dict = Paf._parse_content(fh=paff, min_len=min_len)
        elif isinstance(paf_file, StringIO):
            paf_dict = Paf._parse_content(fh=paf_file, min_len=min_len)
        else:
            paf_dict = dict()
            print("need file path or StringIO")
        return paf_dict



    @staticmethod
    def _parse_content(fh, min_len: int) -> dict:
        """
        Parser for PAF files into defaultdicts

        :param fh: Filehandle of PAF
        :param min_len: minimum alignment block length
        :return: parsed dict with PAF entries
        """
        paf_dict = defaultdict(list)

        for record in fh:
            paf = PafLine(record)
            # FILTERING of PAF ENTRIES
            if paf.alignment_block_length < min_len:
                continue
            if not paf.primary:
                continue

            paf_dict[str(paf.qname)].append(paf)
        return paf_dict



    @staticmethod
    def parse_filter_classify_records(paf: str, filters: SimpleNamespace) -> tuple[list, list]:
        """
        Parse a paf file by converting each line to Paf object,
        while also filtering and classifying mappings

        :param paf: String of PAF file name
        :param filters: Filters stored as Namespace
        :return: Lists of classified and skipped PAF records
        """
        records = []
        skip = []
        with open(paf, 'r') as fh:
            for record in fh:
                rec = PafLine(record)
                # check if this mapping passes the filters
                is_filtered = rec.filter(filters=filters)
                if is_filtered:
                    continue

                # classify the alignment
                rec.c = rec.classify()

                if rec.c == 1:
                    # internal match
                    skip.append(rec)
                    continue
                # if not filtered or skipped
                records.append(rec)
        return records, skip



    @staticmethod
    def choose_best_mapper(records: list) -> list:
        """
        Structured array to decide between ties, by using the score of the DP algorithm

        :param records: List of multiple mappers to decide from
        :return: Best mapper according to attributes
        """
        mapq = [(record.mapq, record.align_score) for record in records]
        custom_dtypes = [('q', int), ('dp', int)]
        mapping_qualities = np.array(mapq, dtype=custom_dtypes)
        sorted_qual = np.argsort(mapping_qualities, order=["q", "dp"])
        record = [records[sorted_qual[-1]]]
        return record



# shorthand typehint used in many places
paf_dict_type = dict[str, list[PafLine]]



