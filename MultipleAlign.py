from Sequence_Align_Toolkit import *
from Sub_Matrix_Gen import *

class MultipleAlignment:
    def __init__(self, seq_list, sub_matrix, gap_penalty=-8):
        self.seq_list = seq_list
        self.sub_matrix = sub_matrix
        self.gap_penalty = gap_penalty
        self.aligned_seqs = []


    def __len__(self):
        """Return the length of aligned sequences"""
        if self.aligned_seqs:
            return len(self.aligned_seqs[0])
        return 0

    def consensus(self):
        """
        Collects frequencies of the characters in each column using dict and selects most common chars, adding them to
        a 'consensus' sequence
        """
        if not self.aligned_seqs:
            return""

        consensus_seq = ""
        for position in range(len(self.aligned_seqs[0])):
            char_freqs = {}
            for seq in self.aligned_seqs:     # iterate through all sequences
                char = seq[position]   # go through every position of all sequences
                if char in char_freqs:
                    char_freqs[char] = char_freqs[char] + 1
                else:                       # if char is found, add 1 to its count
                      char_freqs[char] = 1    # default count of existing char is 1

            maximum_freq = 0
            most_common_char = "+"
            for key, value in char_freqs.items():   # iterate through all keys (which are the characters) of the freq dict
                if key != "-" and char_freqs[key] > maximum_freq:
                    maximum_freq = char_freqs[key] # find the largest max frequency by accessing the char freqeuency counts
                    most_common_char = key
            consensus_seq = consensus_seq + most_common_char
        return consensus_seq

    def global_progressive_align(self):
        if len(self.seq_list) < 2:
            print("You need more than two sequences to run Multiple Sequence Alignment")
            return None

        # start by aligning the first two sequences
        pairwise_aligner = SequenceAlign(self.sub_matrix, self.gap_penalty)  # create SeqAlign obj w/ sub_mat and gap_pen
        pairwise_aligner.set_sequences(self.seq_list[0], self.seq_list[1])   # set seq to first two in list
        pairwise_aligner.gen_global_data()   # generate global / NW algorithm data for nw_result class function

        # initialize aligned sequences list with first pairwise alignment
        self.aligned_seqs = pairwise_aligner.nw_result.copy() # make a copy of NW align results and paste to aligned_seqs

        # progressively add each remaining sequence
        for seq_idx in range(2, len(self.seq_list)):   # iterate through remaining sequences starting from 2
            current_consensus = self.consensus()       # calculate consensus seq from current alignment
            consensus_aligner = SequenceAlign(self.sub_matrix, self.gap_penalty) # creates new SeqAlign obj
            consensus_aligner.set_sequences(current_consensus, self.seq_list[seq_idx]) # performs alignment w/ consensus
            consensus_aligner.gen_global_data() # gen global alignment data

            # get new alignment
            aligned_consensus = consensus_aligner.nw_result[0]
            aligned_new_sequences = consensus_aligner.nw_result[1]
            print("!!", aligned_consensus)
            updated_alignments = [] # create list for updated alignments
            for prev_seq in self.aligned_seqs: # for each previously aligned seq, it tracks the original pos in seq
                updated_seq = ""
                orig_pos = 0

                for cons_pos in range(len(aligned_consensus)): # for each pos in aligned consensus
                    if aligned_consensus[cons_pos] == "-":     # if it has a gap, it's added to the updated sequence
                        updated_seq += "-"
                    else:
                        updated_seq += prev_seq[orig_pos]      # otherwise, it adds the char from original sequence
                        orig_pos += 1
                updated_alignments.append(updated_seq)         # adds updated seq to the list

            updated_alignments.append(aligned_new_sequences)
            self.aligned_seqs = updated_alignments

        return self.aligned_seqs