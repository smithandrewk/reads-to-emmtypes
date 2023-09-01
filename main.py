from lib.utils import *

ids = get_new_sequence_ids()
download_samples(ids)
trim_raw_reads()
qc_trimmed_reads()
assemble_reads()
emmtype_assemblies()