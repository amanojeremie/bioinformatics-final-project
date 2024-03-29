import sys
from itertools import combinations

MATCH_BONUS = 5
MISMATCH_PENALTY = -4
GAP_PENALTY = -11

def valid_sequence(sequence):
	"""
		Determines if a given sequence is a valid DNA sequence
		Valid DNA sequences have at least 1 character and every character
		is a A, C, G or T.
		@param sequence The string to be checked for sequence validity
		@return True if valid False if invalid
	"""
	if len(sequence) == 0:
		return False
	for letter in sequence:
		if letter not in ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']:
			return False
	return True

def align(seq_1_tuple, seq_2_tuple, match_bonus, mismatch_penalty, gap_penalty):
	"""
		Aligns a two strings using a given match bonus mismatch penalty
		and gap penalty. Uses a forward pass to fill table of scores,
		then a backwards pass to create the string alignments
		@param seq_1 The first string to align
		@param seq_2 The second string to align
		@param match_bonus The score if the two strings match at a position
		@param mismatch_penality The score if two strings mismatch at a position
		@param gap_penality The score if a gap must be added to a string
		@return A 5-tuple with the 2 original strings, the 2 aligned strings
			and the score of the alignment
	"""
	seq_1 = seq_1_tuple[1]
	seq_2 = seq_2_tuple[1]
	alignment_table = [[] for _ in range(len(seq_2) + 1)]
	#Create the empty alignment table using the length of the second sequence for the rows plus 1 for the gap
	
	"""
		Forward pass
	"""
	for seq_1_index in range(len(seq_1) + 1):
		alignment_table[0].append(gap_penalty * seq_1_index)
		#Fills the gap penalties along the first row
	
	for seq_2_index in range(1, len(seq_2) + 1):
		alignment_table[seq_2_index].append(seq_2_index * gap_penalty)
		#Fills the gap penalty in the first column
		
		for seq_1_index in range(1, len(seq_1) + 1):
			align = [
				alignment_table[seq_2_index - 1][seq_1_index] + gap_penalty,
				alignment_table[seq_2_index][seq_1_index - 1] + gap_penalty,
				alignment_table[seq_2_index - 1][seq_1_index - 1] +
					(match_bonus if seq_1[seq_1_index - 1] == seq_2[seq_2_index - 1] else mismatch_penalty)
			]
			#Creates an array of three scores, using the column to the left, row at the top, or both
			alignment_table[seq_2_index].append(max(align))
			#Assign the best score in the score array to the cell
	
	"""
		Backwards pass
	"""
	seq_2_index = len(seq_2)
	seq_1_index = len(seq_1)
	aligned_seq_1 = ""
	aligned_seq_2 = ""
	while not (seq_2_index == 0 and seq_1_index == 0):
		#While not at top left corner
		if seq_1_index > 0 and alignment_table[seq_2_index][seq_1_index] - gap_penalty == alignment_table[seq_2_index][seq_1_index - 1]:
			#If the current cell was influenced by the cell on the left
			aligned_seq_1 = seq_1[seq_1_index - 1] + aligned_seq_1
			aligned_seq_2 = "_" + aligned_seq_2
			seq_1_index -= 1
		elif seq_2_index > 0 and alignment_table[seq_2_index][seq_1_index] - gap_penalty == alignment_table[seq_2_index - 1][seq_1_index]:
			#If the current cell was influenced by cell above
			aligned_seq_2 = seq_2[seq_2_index - 1] + aligned_seq_2
			aligned_seq_1 = "_" + aligned_seq_1
			seq_2_index -= 1
		else:
			#If the current cell was influenced by the cell on the top-left
			aligned_seq_1 = seq_1[seq_1_index - 1] + aligned_seq_1
			aligned_seq_2 = seq_2[seq_2_index - 1] + aligned_seq_2
			seq_2_index -= 1
			seq_1_index -= 1

	return (seq_1_tuple, seq_2_tuple, ((seq_1_tuple[0], seq_2_tuple[0], alignment_table[len(seq_2)][len(seq_1)]), aligned_seq_1), 
		((seq_1_tuple[0], seq_2_tuple[0], alignment_table[len(seq_2)][len(seq_1)]), aligned_seq_2), alignment_table[len(seq_2)][len(seq_1)])

def align_sequences(sequences):
	"""
		A recursive function that globally aligns multiple sequences by
		the finding the best global alignment then calling recursively using 
		the with the new alignment with the remaining sequences

		@param sequences An array of sequences to align
		@return The combined global alignment
	""" 
	combos = combinations(sequences, 2)
	#Generates combiniations of the given sequences
	alignments = []
	for combo in combos:
		alignments.append(align(combo[0], combo[1], MATCH_BONUS, MISMATCH_PENALTY, GAP_PENALTY))
	
	if len(alignments) == 1:
		return alignments[0]
	else:
		best_alignment = max(alignments, key = lambda alignment: alignment[4])
		# #Find the highest scoring alignment and removes the sequences that lead to the alignment
		# print(best_alignment[2][0])
		# print("Score: ", best_alignment[4])
		# print() #New line to separate each recursive call
		new_sequences = [best_alignment[2]]
		new_sequences_2 = [best_alignment[3]]
		sequences.remove(best_alignment[0])
		sequences.remove(best_alignment[1])
		new_sequences.extend(sequences)
		new_sequences_2.extend(sequences)
		return max([align_sequences(new_sequences), align_sequences(new_sequences_2)], key = lambda alignment: alignment[4])

def str_insert(str, insert, index):
	"""
		Inserts a string into another string at a given index
	"""
	return str[:index] + insert + str[index:]

def main(argv):
	"""
		Given a FASTA file, globally align all the sequences within
		
		Usage:
		python cs4.py fasta.txt
	"""
	if(len(argv) < 2):
		print("Globally aligns all DNA sequences within a FASTA file")
		print("Usage: python", argv[0], "fasta.txt")
		return
	
	file = open(argv[1], "r")
	sequences = []
	for line in file:
		if line[0] == ">":
			name = line[1:-1] #Omits the > and the newline character
			sequence = file.readline()[:-1].upper() #Omits the newline character and capitilizes max sequence
			if valid_sequence(sequence):
				sequences.append((name, sequence))

	print("Aligning " + str(len(sequences)) + " sequences")
	combined_alignment = align_sequences(sequences.copy())
	print(combined_alignment[2][0])
	print("Final score:", combined_alignment[4])
	print(combined_alignment[2][1])
	print(combined_alignment[3][1])
	print()


	#Realign with final alignment
	aligned_sequences = []
	for sequence in sequences:
		alignment = max(align_sequences([sequence, combined_alignment[2]]),
			align_sequences([sequence, combined_alignment[3]]), key = lambda alignment: alignment[4])
		print(sequence[0] + " Score: " + str(alignment[4]))
		print(alignment[2][1])
		print()
		aligned_sequences.append((sequence[0], alignment[2][1]))

	print("Finding conserved sequences")
	conserved_sequences = []
	latest_window = 0

	#Iterate through the combined alignment
	for window_start in range(0, len(combined_alignment[2][1])):
		
		if window_start < latest_window:
			continue

		isEqual = True
		letters = []

		#Find a index in which every realignment has a matching character 
		for sequence_tuple in aligned_sequences:
			sequence = sequence_tuple[1]
			if sequence[window_start] == "_":
				isEqual = False
				break
			if len(letters) == 0:
				letters.append(sequence[window_start])
			elif letters.count(sequence[window_start]) == 0:
				isEqual = False
				break
		

		if isEqual:
			last_alike = window_start
			#Build the matching string until characters diverge
			for window_end in range(window_start + 1, min([window_start + 25, len(combined_alignment[2][1])])):
				
				subsequence = []
				isEqual = True
				for sequence_tuple in aligned_sequences:
					sequence = sequence_tuple[1]
					if len(subsequence) == 0:
						subsequence.append(sequence[window_start:window_end])
					elif subsequence.count(sequence[window_start:window_end]) == 0:
						isEqual = False
						break
					else:
						latest_window = window_end
						last_alike = window_end
			conserved_sequences.append((combined_alignment[2][1][window_start:last_alike], window_start))

	print("Conserved sequences (sequence, alignment position):")
	print(conserved_sequences)

	#Creates an html file that highlights the conserved sequences
	html_file = open("./output.html", "w")
	for sequence_tuple in aligned_sequences:
		highlighted_sequence = sequence_tuple[1]
		char_additions = 0 #Number of additional characters that have been inserted, to realign with conserved offsets
		for conserved in conserved_sequences:
			new_sequence = str_insert(highlighted_sequence, "<span style='background-color: #ff0000'>", conserved[1] + char_additions)
			char_additions += len(new_sequence) - len(highlighted_sequence)
			highlighted_sequence = new_sequence

			new_sequence = str_insert(highlighted_sequence, "</span>", conserved[1] + len(conserved[0]) + char_additions)
			char_additions += len(new_sequence) - len(highlighted_sequence)
			highlighted_sequence = new_sequence

		print("<p>" + sequence_tuple[0] + "</p>", file=html_file)
		print("<p style='font-family: Consolas, Courier'>" + highlighted_sequence + "</p>", file=html_file)
	
	html_file.close()
		
if __name__ == "__main__":
	main(sys.argv)