import mappy
import sys

# fasta file is read as positional argument
ref = sys.argv[1]
# name for index file
idx = f"{ref}.mmi"
# calling mappy with fn_idx_out writes index to file
aligner = mappy.Aligner(ref, preset="map-ont", n_threads=24, fn_idx_out=idx)
