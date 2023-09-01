from lib.utils import *

# ids = get_new_sequence_ids()
# download_samples(ids)
# trim_raw_reads()
# qc_trimmed_reads()
# assemble_reads()
# emmtype_assemblies()
from staphb_toolkit.lib import calldocker as container_engine
from staphb_toolkit.lib.autopath import path_replacer
import staphb_toolkit.lib.container_handler as container

arg_string,path_map = path_replacer(['data/4_assembled/STREP22-0001-SC-M05253-230701_S20/contigs.fasta'],os.getcwd())
command = f'emmtyper data/4_assembled/STREP22-0001-SC-M05253-230701_S20/contigs.fasta -o out.tsv'
program_object = container.Run(command=command, path=path_map, image='staphb/emmtyper', tag='latest')
program_object.run()