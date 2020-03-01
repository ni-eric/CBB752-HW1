#!/usr/bin/python
__author__ = "Eric Ni"
__email__ = "eric.ni@yale.edu"
__copyright__ = "Copyright 2020"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage:      python hw1.py -i <input file> -s <score file> -o <open gap penalty> -e <extend gap panalty> -O <output file name>
### Example:    python hw1.py -i input.txt -s blosum62.txt
### Note:       Smith-Waterman Algorithm

import argparse
import numpy as np
### This is one way to read in arguments in Python.
### We need to read input file and score file.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
parser.add_argument('-O', '--output', help='output file', required=False, default='hw1_output.txt')
args = parser.parse_args()

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap, outFile):
	### Print input and score file names. You can comment these out.
	print ("input file : %s" % inputFile)
	print ("score file : %s" % scoreFile)
	print ("open gap penalty : %s" % openGap)
	print ("extension gap penalty : %s" % extGap)
	openGap = int(openGap)
	extGap = int(extGap)
	# read input file into separate string files
	with open(inputFile,'r') as i:
		seqs = i.read().splitlines()
	template_res = seqs[0]
	target_res = seqs[1]

	# read first line of score file for indexing
	with open(scoreFile,'r') as s:
		res = s.readlines()[0].split()

	# require there to be exactly 2 sequences, 1 per line
	if len(seqs)!=2:
		raise ValueError('input file must contain exactly 2 sequences')
	else:
	# convert sequence amino acids to index (so this should work with any matrix:DNA etc)
		template = []
		for aa in template_res:
			template.append(res.index(aa))
		target = []
		for aa in target_res:
			target.append(res.index(aa))

	# read score file as substitution matrix 
	with open(scoreFile,'r') as s:
		ncols = len(s.readline().split())
		submat = np.loadtxt(s,usecols=range(1,ncols+1),skiprows=0,dtype=int)
	#print submat

	# fill in the scoring matrix according to Smith Waterman algorithm
	rows = len(target)
	cols = len(template)
	# declare scoring matrix, add 1 to row/col for initial row/col of 0's
	SWmat = [[0 for i in range(cols+1)] for j in range(rows+1)]
	# declare arrow matrix of same size showing source of each square for each cell
	arrowmat = [[[] for i in range(cols+1)] for j in range(rows+1)]

	SWmat = np.array(SWmat)
	for n in range(1,cols+1):
		for m in range(1,rows+1):
			#print template[n]
			#print target[m]

			# diagonal by 1 cell
			diagonal = SWmat[m-1,n-1] + submat[target[m-1]][template[n-1]]
			# creating gap along m
			gap_penalty = np.array([extGap*(m-k-1)+openGap for k in range(m)])
			m_gap = SWmat[0:m,n]+gap_penalty
			# creating gap along n
			gap_penalty = np.array([extGap*(n-k-1)+openGap for k in range(n)])
			n_gap = SWmat[m,0:n]+gap_penalty

			maxscore = max(diagonal,max(m_gap),max(n_gap),0)
			SWmat[m][n] = maxscore

			# append coordinates to arrow matrix according to where score came from
			# should support any amount of arrows
			if(maxscore==diagonal):
				arrowmat[m][n].append([m-1,n-1])
			if(maxscore==max(m_gap)):
				k = np.where(m_gap==maxscore)[0]
				for i in k:
					arrowmat[m][n].append([i,n])
			if(maxscore==max(n_gap)):
				k = np.where(n_gap==maxscore)[0]
				for i in k:
					arrowmat[m][n].append([m,i])

	# write to output file 
	f = open(outFile,'w')

	f.write('************************************************************\n')
	f.write('*    Input Sequences                                       *\n')
	f.write('************************************************************\n\n')
	f.write('Sequence 1\t'+''.join(template_res)+'\n')
	f.write('Sequence 2\t'+''.join(target_res)+'\n\n')
	f.write('************************************************************\n')
	f.write('*    Score Matrix                                          *\n')
	f.write('************************************************************\n\n')
	# print scoring matrix with sequences 
	print '\t\t',				
	f.write('\t\t')
	for aa in template_res:
		print aa,'\t',
		f.write(aa+'\t')
	print
	f.write('\n')
	for num,row in enumerate(SWmat):
		if num>0:
			print target_res[num-1],'\t',
			f.write(target_res[num-1]+'\t')
		else:
			print '\t',
			f.write('\t')
		for sc in row:
			print sc,'\t',
			f.write(str(sc)+'\t')
		print
		f.write('\n')

	
	#for row in arrowmat:
	#	print row

	# declare empty alignments for the final alignment output
	template_align = ['[',']']
	lines_align = [' ',' ']
	target_align = ['[',']']


	f.write('\n************************************************************\n')
	f.write('*    Best Local Alignment                                  *\n')
	f.write('************************************************************\n\n')

	# find maximum of entire array
	[m,n] = np.where(SWmat==SWmat.max())
	end_m = m
	end_n = n
	target_align.insert(1,target_res[m-1])
	template_align.insert(1,template_res[n-1])
	print 'Maximum Alignment Score:',SWmat[m,n][0]
	f.write('Maximum Alignment Score:'+str(SWmat[m,n][0])+'\n')
	# perform traceback
	while (SWmat[m,n]!=0):
		# find where arrows are pointing from that cell
		[a,b] = arrowmat[m][n][0]
		# if multiple equivalent arrows exist, my implementation will choose diagonal>up>left
		# these will produce equal scoring alignments, however i will only choose one of them for output
		if(a==m-1 and b==n-1):	
			target_align.insert(1,target_res[a])
			template_align.insert(1,template_res[b])
			if(target_res[a]==template_res[b]):
				lines_align.insert(1,'|')
			else:
				lines_align.insert(1,' ')
		elif(a<m and b==n):
			for i in range(m-1,a-1,-1):
				target_align.insert(1,target_res[i])
				template_align.insert(1,'-')
				lines_align.insert(1,' ')
		elif(a==m and b<n):
			for i in range(n-1,b-1,-1):
				target_align.insert(1,'-')
				template_align.insert(1,template_res[i])
				lines_align.insert(1,' ')
		else:
			print('something went wrong')

		m = a
		n = b
		
	# pop out the end, since it stops at 0 score
	template_align.pop(-2)
	target_align.pop(-2)

	# output for final 
	print('Final Best Local Alignment')
	blank = ""

	n_start = template_res[0:n]
	m_start = (target_res[0:m])
	n_end =(template_res[end_n:])
	m_end = (target_res[end_m:])

	
	front_gap = len(max([n_start,m_start],key=len))
	n_start = blank.join([' ' for i in range(0,front_gap-len(n_start))])+n_start
	m_start = blank.join([' ' for i in range(0,front_gap-len(m_start))])+m_start
	blank_gap = blank.join([' ' for i in range(0,front_gap)])


	print n_start + blank.join(template_align) + n_end
	print blank_gap+ blank.join(lines_align)
	print m_start + blank.join(target_align) + m_end

	f.write(n_start + blank.join(template_align) + n_end+'\n')
	f.write(blank_gap+ blank.join(lines_align)+'\n')
	f.write(m_start + blank.join(target_align) + m_end+'\n')





### Run MY Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap, args.output)