# 2018-10-12 fix for PurgeConsensusReadAlpha

# Authors: L. Penter, L. Hansmann
# Hematology, Oncology, and Tumor Immunology, Charité - Universitätsmedizin Berlin

import ConsensusClusters
import CdrExtractionOptions

import copy
import numpy
import tempfile
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from io import StringIO

# extract cytokine from read
def CytokineExtraction(SEQ, DB):
	# blast read SEQ against local cytokine database DB
	blastx_cline = NcbiblastnCommandline(cmd='blastn', db=DB, outfmt=5, evalue=0.1, task='blastn', num_threads=8, word_size=7)
	out, err = blastx_cline(stdin=SEQ)
	out_xml = StringIO(out)
	blast_records = NCBIXML.read(out_xml)

	# if no cytokine found, return empty string
	if len(blast_records.alignments) == 0:
		return ''

	# if cytokine found, return name of cytokine
	if blast_records.alignments[0].hsps[0].expect < 1e-10:
		return  blast_records.alignments[0].title.split(' ')[1].split('_')[0] # cytokine

	return ''

# generate output of cytokines
def CytokineOutput(wellname, cytokine_list):
	output = wellname + '\t'
	for k in cytokine_list:
		output += str(cytokine_list[k]) + '\t'
	return output[:-1] + '\n'

# generate file without cytokine reads and return temporary file plus cytokine list
def FileWithoutCytokines(filename):
    # generate temporary file
    tmp_file = tempfile.NamedTemporaryFile(delete=False)
    cytokine_list = copy.deepcopy(CdrExtractionOptions.CYTOKINE_LIST)

    # read in all reads from sequence file 
    sequences = ConsensusClusters.ReadSequences(filename)
        
    # go through all reads
    for s in sequences:
        # check if read contains cytokine and count it
        cytokine = CytokineExtraction(sequences[s], CdrExtractionOptions.PATH_TO_CYTOKINE_DB)
        if cytokine != '':
            cytokine_list[cytokine]+=1
        # if it does not contain cytokine, write to file
        else:   
            tmp_file.write('>' + s + '\n' + sequences[s] + '\n')

    # close temporary file
    tmp_file.close()

    return tmp_file.name, cytokine_list

# generate input for IMGT HighV-Quest
# >wellname:index:number of reads
def HighV_QuestInput(wellname, possible_TCR_reads):
	imgt_input = ''
	counter = 0
	for p in possible_TCR_reads:
		imgt_input += '>' + wellname + ':' + str(counter) + ':' + str(p[1]) + '\n' + p[0] + '\n'
		counter += 1
	return imgt_input

# process all reads from one file, return list of possible alpha/beta chains and cytokine information
def ParseWell(filename):
	cytokine_list = copy.deepcopy(CdrExtractionOptions.CYTOKINE_LIST)
	possible_TCR_reads = []
	possible_TCR_reads_ = []

	# gets clusters
	clusters = ConsensusClusters.ConsensusClusters(filename)
	# go through clusters
	for cluster in clusters:
		# check if read contains cytokine
		cytokine = CytokineExtraction(str(cluster[0]), CdrExtractionOptions.PATH_TO_CYTOKINE_DB)
		if cytokine != '':
			cytokine_list[cytokine]+=cluster[1]
		# if it does not contain cytokine, submit to IMGT
		else:
			if cluster[1] > CdrExtractionOptions.MINIMUM_READS_FOR_CLUSTER and \
			len(cluster[0]) > CdrExtractionOptions.MINIMUM_LENGTH_OF_TCR_READ:
				possible_TCR_reads.append([str(cluster[0]), cluster[1]])
				possible_TCR_reads_.append(cluster[1])

	# find reads with most reads, return only as many as specified in MAXMIMUM_NUMBER_OF_READS_PER_WELL
	p = []
	for i in numpy.argsort(possible_TCR_reads_)[-CdrExtractionOptions.MAXIMUM_NUMBER_OF_READS_PER_WELL:][::-1]:
		p.append(possible_TCR_reads[i])

	return p, cytokine_list

# purge barcodes and primers of alpha read 
def PurgeConsensusReadAlpha(s):
    start = s.upper().find('CCAGGGTTTTCCCAGTCACGAC') + 22
    end = s.upper().find('CTGAGAGACTCTAAATCCAGTGAC')

    if (start == -1 or end == -1):
        return s

    return s[start:end]

# purge barcodes and primers of beta read
def PurgeConsensusReadBeta(s):
    start = s.upper().find('CCAGGGTTTTCCCAGTCACGAC') + 22
    end = s.upper().find('GAGCCATCAGAAGCAGAGATCTC') 
    if (start == -1 or end == -1):
        return s

    return s[start:end]
