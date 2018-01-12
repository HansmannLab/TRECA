# configuration file for CdrExtraction 
from collections import OrderedDict

# parameters are optimized for single cell analysis 

# link to cd-hit-est binary
CD_HIT_EST = '/usr/local/bin/cd-hit-est'
# command line options for cd-hit-est
CD_HIT_EST_OPTIONS = ['-M','0','-T','0','-l','50']
# maximum number of reads to be used for generation of consensus read
MAX_READS_FOR_CONSENSUS = 1000 
# minimum reads for a sequencing cluster 
MINIMUM_READS_FOR_CLUSTER = 5
# minimum length of possible TCR reads
MINIMUM_LENGTH_OF_TCR_READ = 150
# maxmimum number of possible TCR sequences to be submitted to IMGT
MAXIMUM_NUMBER_OF_READS_PER_WELL = 5

# databases for blast as provided in subdirectory ./references - absolute path needs to be set
# Cytokines
# PATH_TO_CYTOKINE_DB='./references/cytokine_reference.fasta'
PATH_TO_CYTOKINE_DB='/Users/shaka87/bioinformatics/t_cells/cytokine_reference.fasta'
# V-Regions
# PATH_TO_TCR_DB='./references/tcr_cyt_ref_new.fa'
PATH_TO_TCR_DB='/Users/shaka87/bioinformatics/t_cells/tcr_cyt_ref_new.fa'
# available cytokines
CYTOKINE_LIST = OrderedDict([('IL2',0),('IL10',0),('IL12A',0),('IL13',0),('IL17A',0),('IL21',0),('IFNG',0),('TNF',0),('TGFB1',0),('GZMB',0),('PRF1',0),('GATA3',0),('RORC',0),('FOXP3',0),('BCL6',0),('RUNX1',0),('RUNX3',0),('TBX21',0)])
