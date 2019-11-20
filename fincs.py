import numpy as np

'''sdf
readfasta - reads file that is exactly fasta format
Parameter: 
    filename - name of the input file
Return:
	sequences - a list of 2-tuples that contain the names and sequences
'''
def readfasta(filename):
	fastafile = open(filename, 'r')
	lines = fastafile.readlines()
	i = 0
	sequences = []
	while i < len(lines) - 1:
		lines[i] = lines[i].rstrip()
		lines[i + 1] = lines[i + 1].rstrip()
		if len(lines[i]) > 0:
			if lines[i][0] == '>':
				sequences.append([lines[i][1:], lines[i+1].upper()])
				i += 1
		i += 1
	fastafile.close()
	return sequences
	
'''
align - uses dynamic programming to align two sequencs
Parameter:
	seq1 - first sequence to align
	seq2 - second sequence to align
Return:
	alnscore - the alignment score for this alignment
	alnseq1 - the first aligned sequence
	alnseq2 - the second aligned sequence
	
'''
def align(seq1, seq2):
	#the penalties for the alignment matrix are set
	gap = -7
	match = 5
	mismatch = -4
	#include the gap in the first position of the alignment matrix
	row = len(seq1) + 1
	col = len(seq2) + 1
	#initialize alignment matrix
	alnmat = np.zeros((row, col))
	
	#below is the forward pass through the matrix
	
	#score the first column
	for i in range (1, row):
		alnmat[i][0] = alnmat[i-1][0] + gap
	#score first row
	for i in range (1, col):
		alnmat[0][i] = alnmat[0][i - 1] + gap
	#score all other values from upper left to lower right
	for i in range (1, row):
		for j in range (1, col):
			diag = 0
			if seq1[i-1] == seq2[j-1]:
				diag = match
			elif seq1[i-1] == '-' or seq2[j-1] == '-':
				diag = gap
			else:
				diag = mismatch
			scores = [alnmat[i - 1][j] + gap, alnmat[i][j-1] + gap, alnmat[i - 1][j - 1] + diag]
			alnmat[i][j] = max(scores)

	#below is the backward pass through the matrix
	
	
	i = len(seq1)-1
	j = len(seq2)-1
	alnseq1 = ''
	alnseq2 = ''
	#check where current score came from and create alignment sequences
	while i > 0 or j > 0:
		if alnmat[i][j] - gap == alnmat[i-1][j]:
			alnseq1 = '-' + alnseq1
			alnseq2 = seq2[j] + alnseq2
			i -= 1
		elif alnmat[i][j] - gap == alnmat[i][j-1]:
			alnseq1 = seq1[i] + alnseq1
			alnseq2 = '-' + alnseq2
			j -= 1
		else:
			alnseq1 = seq1[i] + alnseq1
			alnseq2 = seq2[j] + alnseq2
			i -= 1
			j -= 1
	alnscore = alnmat[row-1][col-1]
	return [alnscore, alnseq1, alnseq2]

#read fasta file and store in sequencs
if __name__ == '__main__':
    infile = 'enolasefasta.txt'
    sequences = readfasta(infile)
aln = []

print 'We will find the best pairwise alignment, remove the original two sequences, '\
	'and replace them with the upper sequence of this alignment to find the next best alignment.'
count = 1
alignments = []
alignments.append(sequences[0][1])
del sequences[0]
#iterate through all sequences
while len(sequences) > 0:
	aln = []
	hiscore = 0
	alnseq1 = ''
	alnseq2 = ''
	besti = -1
	bestj = -1
	#iterate through each combination of sequences
	for i in range (0, len(alignments)):
		for j in range (0, len(sequences)):
			#store best alignment sequences
			score = 0
			aln = align(alignments[i][1], sequences[j][1])
			score = aln[0]
			if score > hiscore:
				alnseq1 = aln[1]
				alnseq2 = aln[2]
				hiscore = score
				besti = i
				bestj = j
	#remove the sequences used for the best alignment from list
	#del alignments[besti]
	del sequences[bestj]
	print
	print 'Alignment ', count
	print alnseq1
	print
	print alnseq2
	print
	#alignments.append(alnseq1)
	alignments.append(alnseq2)
	count += 1
print
print "These are the aligned sequences: "
print
for i in range(0,len(alignments)):
	print
	print alignments[i]
	print
# #align the final two alignment sequences
# aln = align(sequences[0][1], sequences[1][1])
# score = aln[0]
# top = aln[1]
# bottom = aln[2]
# print
# print "Final Alignment"
# print top
# print
# print bottom


	
	