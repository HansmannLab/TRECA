import os
import random
import subprocess
import tempfile

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align.Applications import MuscleCommandline

import CdrExtractionOptions

# return consensus sequences for each cluster in file
def ConsensusClusters(filename):
	clusters = []
	cluster_counter = -1

	# read in sequencing reads from file
	sequences = ReadSequences(filename)

	# generate temporary file
	tmp_file = tempfile.NamedTemporaryFile(delete=False)
	tmp_file.close()
	# calculate clusters of sequencing reads using cd-hit-est
	subprocess.check_output([CdrExtractionOptions.CD_HIT_EST, '-i', str(filename), '-o', str(tmp_file.name)] + CdrExtractionOptions.CD_HIT_EST_OPTIONS)
	
	# process cluster information
	f = open(tmp_file.name + '.clstr', 'rU')
	for line in f:
		# start of new cluster
		if line[0] == '>':
			clusters.append([])
			cluster_counter += 1
		else:
			sequence_number = line[line.find('>')+1:line.find('...')]
			clusters[cluster_counter].append(sequences[sequence_number])

	# delete temporary files
	os.unlink(tmp_file.name)
	os.unlink(tmp_file.name+'.clstr')

	return ConsensusSequences(clusters)

# generate consensus sequences for groups of reads      
def ConsensusSequences(clusters):
	reads = []
		
	# go through each cluster of sequencing reads provided
	for g in clusters: 
		# temporary files
		tmp_file = tempfile.NamedTemporaryFile(delete=False)
		tmp_file2 = tempfile.NamedTemporaryFile(delete=False)
		tmp_file2.close()
		# if group is larger than MAX_READS_FOR_CONSENSUS reads, select random sample of MAX_READS_FOR_CONSENSUS
		if len(g) > CdrExtractionOptions.MAX_READS_FOR_CONSENSUS:
			subsample = random.sample(g, CdrExtractionOptions.MAX_READS_FOR_CONSENSUS)
		else:
			subsample = g # use all available reads	
		# store cluster in temporary file
		for i in range(0, len(subsample)):
			tmp_file.write ('>' + str(i) + '\n')
			tmp_file.write (subsample[i] + '\n')
		tmp_file.close()
		# generate consensus read using muscle
		muscle_cline = MuscleCommandline(input=tmp_file.name, out=tmp_file2.name, maxiters=2)
		muscle_cline()
		align = AlignIO.read(tmp_file2.name, 'fasta')
		summary_align = AlignInfo.SummaryInfo(align)
		# store consensus read in list
		reads.append([summary_align.gap_consensus(ambiguous=''), len(g)])
		
		# remove tmp files
		os.unlink(tmp_file.name)
		os.unlink(tmp_file2.name)

	return reads

# read all sequences from fasta file
def ReadSequences(filename):
	sequences = {}
	identifier = ''

	f = open(filename, 'rU')
	for line in f:
		if line[0] != '>':
			sequences[identifier] = line.strip()
		else:
			identifier = line.strip()[1:]

	return sequences
