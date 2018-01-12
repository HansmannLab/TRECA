#!/usr/bin/python

# extract consensus sequences from fasta files, each containing all sequences from one individual well
# all files in current directory ending with .fasta will be processed
# 
# sequences from each fasta file are divided into groups of similar sequences using CD-HIT 
# using these groups of sequences consensus sequences are generated with Muscle
# consensus sequences are written to output file serving as input for IMGT
 
import CdrExtraction
import CdrExtractionOptions

import argparse
import glob
import multiprocessing
import os
import subprocess

# sub-thread for processing one file
def parse_file(filename, blast_cytokines):
        # if option for blasting each read to identify cytokine reads is set
        # generate temporary file without cytokine reads and proceed with it
        if blast_cytokines:
            tmp_file, cytokine_list = CdrExtraction.FileWithoutCytokines(filename)
            possible_TCR_list, empty_cytokine_list = CdrExtraction.ParseWell(tmp_file)
            os.unlink(tmp_file) 
        else:
	    possible_TCR_list, cytokine_list = CdrExtraction.ParseWell(filename)

        out_imgt.write(CdrExtraction.HighV_QuestInput(filename.split('.')[0], possible_TCR_list))
	out_cytokine.write(CdrExtraction.CytokineOutput(filename.split('.')[0], cytokine_list))

# main thread
if __name__ == '__main__':
	# parsing arguments
	parser = argparse.ArgumentParser(description='Process files containing sequencing reads.')
	parser.add_argument('--imgt_input', required=True, help='File that will contain input for IMGT High/V-Quest')
	parser.add_argument('--cytokine_output', required=True, help='File that will contain output of cytokine reads')
        parser.add_argument('-b','--blast_cytokines', help='Blast each read to identify cytokine reads', action='store_true')
	args = parser.parse_args()

	# files to be written
	out_imgt = open(args.imgt_input, 'w',0)
	out_cytokine = open(args.cytokine_output, 'w',0)
	out_cytokine.write ('Well\t' + '\t'.join(CdrExtractionOptions.CYTOKINE_LIST.keys()) + '\n')

	# starting sub-threads
	pool = multiprocessing.Pool()
	for filename in sorted(glob.glob('*.fasta')):
            pool.apply_async(parse_file, args=(filename, args.blast_cytokines))

	# clean up
	pool.close()
	pool.join()
	out_imgt.close()
	out_cytokine.close()
