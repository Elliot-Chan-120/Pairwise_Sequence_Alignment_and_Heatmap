import plotly.express as px

class SequenceAlign:
    def __init__(self, seq_1, seq_2, sub_matrix, gap_penalty):
        self.sub_matrix = sub_matrix
        self.gap_penalty = gap_penalty
        self.seq_1 = seq_1.upper()
        self.seq_2 = seq_2.upper()

        # Store Needleman-Wunsch data
        self.nw_s_matrix = None
        self.nw_t_matrix = None
        self.nw_optimal_alignment = None
        self.nw_result = []
        self.nw_path_matrix = None

        # Store Smith-Waterman data
        self.sw_s_matrix = None
        self.sw_t_matrix = None
        self.sw_optimal_alignment = None
        self.sw_result = []
        self.sw_path_matrix = None

# DATA GENERATION FUNCTIONS
    def gen_alignment_data(self):
        """
        Run this first.
        Generates all needed data to use the rest of the functions in this class.
        """
        self.gen_local_data()
        self.gen_global_data()

        self.global_path_matrix()
        self.local_path_matrix()

    def needleman_wunsch(self):
        """
        Builds 2 matrices:
        score_matrix - contains best score possible for aligning sequences up to each position
        traceback_matrix - contains which move contains the best score (1=diagonal, 2=up, 3=left)
        """
        score_matrix = [[0]]  # sorting matrix - contains best alignment scores for each pos
        traceback_matrix = [[0]]  # traceback matrix - records moves that gave the best scores

        # initialize gap row (all gaps in seq 1)
        for col in range(1, len(self.seq_2) + 1):
            score_matrix[0].append(self.gap_penalty * col)
            traceback_matrix[0].append(3)

        # initialize gap column (all gaps in seq 2)
        for row in range(1, len(self.seq_1) + 1):
            score_matrix.append([self.gap_penalty * row])
            traceback_matrix.append([2])

        # apply recurrence relation to fill remaining of matrix
        for row in range(0, len(self.seq_1)):
            for col in range(len(self.seq_2)):
                # calculate scores for 3 possible moves:
                diagonal_score = score_matrix[row][col] + self.score_pos(self.seq_1[row], self.seq_2[col], self.sub_matrix, self.gap_penalty)  # Diagonal
                up_score = score_matrix[row][col + 1] + self.gap_penalty  # Up
                left_score = score_matrix[row + 1][col] + self.gap_penalty  # Left
                # choose best score and move
                score_matrix[row + 1].append(max(diagonal_score, up_score, left_score))
                traceback_matrix[row + 1].append(self.max3t(diagonal_score, up_score, left_score))  # record which move it was using max3t function
        return score_matrix, traceback_matrix

    def smith_waterman(self):
        score_matrix = [[0]]
        traceback_matrix = [[0]]
        maxscore = 0

        for col in range(1, len(self.seq_2) + 1):
            score_matrix[0].append(0)
            traceback_matrix[0].append(0)

        for row in range(1, len(self.seq_1) + 1):
            score_matrix.append([0])
            traceback_matrix.append([0])

        for row in range(0, len(self.seq_1)):
            for col in range(len(self.seq_2)):
                # calculate the scores for 3 possible moves
                diagonal_score = score_matrix[row][col] + self.score_pos(self.seq_1[row], self.seq_2[col], self.sub_matrix, self.gap_penalty)
                up_score = score_matrix[row][col + 1] + self.gap_penalty
                left_score = score_matrix[row + 1][col] + self.gap_penalty
                best_score = max(diagonal_score, up_score, left_score)
                if best_score <= 0: # calculate and see if best score is greater than 0
                    score_matrix[row + 1].append(0)
                    traceback_matrix[row + 1].append(0)  # If everything is smaller than 0, then 0 is used as the score instead
                else:
                    score_matrix[row + 1].append(best_score) # Otherwise, best score will be used in matrices
                    traceback_matrix[row + 1].append(self.max3t(diagonal_score, up_score, left_score))
                    if best_score > maxscore:
                        maxscore = best_score
        return score_matrix, traceback_matrix, maxscore

    @staticmethod
    def score_pos(char_1, char_2, sub_matrix, gap_penalty):
        """
        Scores a single position in alignment
        If any character is a gap -> gap penalty
        Otherwise -> look up the score for this pair of AA/NT
        """
        if char_1 == "−" or char_2 == "−":  # if either position is a gap
            return gap_penalty  # return gap penalty
        else:
            return sub_matrix[char_1 + char_2]  # otherwise look up score in BLOSUM62 matrix

    @staticmethod
    def max3t(diagonal_score, up_score, left_score):
        """
        Returns the value which is largest out of the three, and sets it to 1, 2 or 3 depending on which direction the
        value was associated with. This will be relevant in the traceback and path matrices when we need to
        reconstruct the alignment and understand where the optimal path alignment goes towards.
        :param diagonal_score:
        :param up_score:
        :param left_score:
        :return:
        """
        if diagonal_score > up_score:
            if diagonal_score > left_score:
                return 1
            else:
                return 3
        else:
            if up_score > left_score:
                return 2
            else:
                return 3

    def score_align(self, seq1, seq2, submat, gap_penalty):
        """
        Calculates total score for entire alignment by summing them for each position
        """
        aligned_sequences = 0
        for i in range(len(seq1)):
            aligned_sequences += self.score_pos(seq1[i], seq2[i], submat, gap_penalty)
        return aligned_sequences

    @staticmethod
    def recover_align_global(traceback_matrix, seq1, seq2):
        """
        Starts at bottom right corner and follows traceback matrix to reconstruct alignment for Needleman-Wunsch
        """
        aligned_seqs = ["", ""]  # Will hold two aligned sequences
        row = len(seq1)  # Start at bottom right
        col = len(seq2)

        while row > 0 or col > 0:
            if traceback_matrix[row][col] == 1:  # Diagonal move
                aligned_seqs[0] = seq1[row - 1] + aligned_seqs[0]  # Match/Mismatch two amino acids
                aligned_seqs[1] = seq2[col - 1] + aligned_seqs[1]
                row -= 1
                col -= 1
            elif traceback_matrix[row][col] == 3:  # Left move
                aligned_seqs[0] = "-" + aligned_seqs[0]  # Insert gap in first sequence
                aligned_seqs[1] = seq2[col - 1] + aligned_seqs[1]  # Character from second sequence
                col -= 1
            else:  # Up move
                aligned_seqs[0] = seq1[row - 1] + aligned_seqs[0]  # Insert gap in second sequence
                aligned_seqs[1] = "-" + aligned_seqs[1]  # Character from first sequence
                row -= 1
        return aligned_seqs

    def recover_align_local(self, score_matrix, traceback_matrix, seq1, seq2):
        """Reconstructs alignment for Smith-Waterman starting from highest score until 0 is reached"""
        aligned_seqs = ["", ""]
        current_row, current_col = self.max_mat(score_matrix) # start at the highest score

        while traceback_matrix[current_row][current_col] > 0: # stop when 0 is reached
            move = traceback_matrix[current_row][current_col]
            if move == 1:           # diagonal move
                aligned_seqs[0] = seq1[current_row - 1] + aligned_seqs[0]
                aligned_seqs[1] = seq2[current_col - 1] + aligned_seqs[1]
                current_row -= 1
                current_col -= 1
            elif move == 3:         # left move -> insert gap in the first seq
                aligned_seqs[0] = "-" + aligned_seqs[0]
                aligned_seqs[1] = seq2[current_col - 1] + aligned_seqs[1]
                current_col -= 1
            elif move == 2:         # up move -> insert gap in the second seq
                aligned_seqs[0] = seq1[current_row - 1] + aligned_seqs[0]
                aligned_seqs[1] = "-" + aligned_seqs[1]
                current_row -= 1
        return aligned_seqs

    @staticmethod
    def max_mat(mat):
        """
        Obtain the maximum value from the matrix
        """
        maxval = mat[0][0]
        maxrow = 0
        maxcol = 0
        for i in range(0, len(mat)):
            for j in range(0, len(mat[i])):
                if mat[i][j] > maxval:
                    maxval = mat[i][j]
                    maxrow = i
                    maxcol = j
        return maxrow, maxcol

    @staticmethod
    def print_matrix(sub_matrix):
        """
        Neatly prints out any matrix
        """
        for i in range(0, len(sub_matrix)):
            print(" ".join(map(str, sub_matrix[i])))

    def gen_global_data(self):
        """
        FETCHES + GEN. RAW DATA
        stores optimal alignments, score + traceback matrices
        :return:
        """
        nw_raw_data = self.needleman_wunsch()
        self.nw_s_matrix = nw_raw_data[0]
        self.nw_t_matrix = nw_raw_data[1]
        self.nw_optimal_alignment = self.nw_s_matrix[len(self.seq_1)][len(self.seq_2)]
        align_global = self.recover_align_global(self.nw_t_matrix, self.seq_1, self.seq_2)
        self.nw_result = [align_global[0], align_global[1]]

    def gen_local_data(self):
        """
        FETCHES + GEN. RAW DATA to store in class variables
        stores optimal alignments, score + traceback matrices
        """
        sw_raw_data = self.smith_waterman()
        self.sw_s_matrix = sw_raw_data[0]
        self.sw_t_matrix = sw_raw_data[1]
        self.sw_optimal_alignment = sw_raw_data[2]
        align_local = self.recover_align_local(self.sw_s_matrix, self.sw_t_matrix, self.seq_1, self.seq_2)
        self.sw_result = [align_local[0], align_local[1]]

    def global_path_matrix(self):
        """
        Generates matrix showing the path of optimal alignment via 1's and 0's
        \n 1 indicates cell 'travelled' and 0 means not taken
        """
        trace_matrix = self.nw_t_matrix # copy traceback matrix since we're going to build on that
        row = len(self.seq_1) + 1
        col = len(self.seq_2) + 1

        self.nw_path_matrix = [[0] * col for _ in range(row)]

        # start from bottom right corner -> needleman-wunsch
        r, c = row - 1, col - 1

        while r > 0 or c > 0: # iterate through matrix
            self.nw_path_matrix[r][c] = 1 # fill target w/ a 1
            if trace_matrix[r][c] == 1: # if cell has value 1, move diagonal
                r -= 1
                c -= 1
            elif trace_matrix[r][c] == 2: # if cell has value 2, move up
                r -= 1
            else:
                c -= 1 # if cell has value 3, move left
        return self.nw_path_matrix

    def local_path_matrix(self):
        """
        Generates matrix showing the path of optimal alignment via 1's and 0's
        \n 1 indicates cell 'travelled' and 0 means not taken
        """
        trace_matrix = self.sw_t_matrix # this function needs both traceback and score
        score_matrix = self.sw_s_matrix
        row = len(self.seq_1) + 1
        col = len(self.seq_2) + 1

        self.sw_path_matrix = [[0] * col for _ in range(row)]

        r, c = self.max_mat(score_matrix) # start at highest-scoring cell as per smith-waterman algorithm

        while (r > 0 and c > 0) and trace_matrix[r][c] > 0: # iterate through the matrix and stop once 0 is found
            self.sw_path_matrix[r][c] = 1
            if trace_matrix[r][c] == 1: # move diagonal
                r -= 1
                c -= 1
            elif trace_matrix[r][c] == 2: # move up
                r -= 1
            else: # move left
                c -= 1
        return  self.sw_path_matrix


# DATA OUTPUT FUNCTIONS
    def output_nw_raw_matrices(self):
        """Outputs Score and Traceback Matrices from Needleman-Wunsch Algorithm"""
        print("\n===[GLOBAL / Needleman-Wunsch] Raw Data===")

        print("--Score Matrix--")
        self.print_matrix(self.nw_s_matrix)

        print("\n--Traceback Matrix--")
        self.print_matrix(self.nw_t_matrix)

    def output_sw_raw_matrices(self):
        """Outputs Score and Traceback Matrices from Smith-Waterman Algorithm"""
        print("\n===[Local / Smith-Waterman] Raw Data===")
        print("--Score Matrix--")
        self.print_matrix(self.sw_s_matrix)

        print("\n--Traceback Matrix--")
        self.print_matrix(self.sw_t_matrix)

    def output_nw_path(self):
        """Outputs Path Matrix from Needleman-Wunsch Algorithm"""
        print("\n===[Global / Needleman-Wunsch] Path===")
        self.print_matrix(self.nw_path_matrix)

    def output_sw_path(self):
        """Outputs Path Matrix from Smith-Waterman Algorithm"""
        print("\n===[Local / Smith-Waterman] Path===")
        self.print_matrix(self.sw_path_matrix)

    def output_nw_results(self):
        """Outputs Aligned Sequences and Optimal Alignment Score from Needleman-Wunsch Algorithm"""
        print("\n[[GLOBAL / Needleman-Wunsch Aligned Sequences]]:")
        print("Optimal Score:", self.nw_optimal_alignment)
        for item in self.nw_result:
            print(item)

    def output_sw_results(self):
        """Outputs Aligned Sequences and Optimal Alignment Score from Smith-Waterman Algorithm"""
        print("\n[[LOCAL / Smith-Waterman Aligned Sequences]]:")
        print("Optimal Score:", self.sw_optimal_alignment)
        for item in self.sw_result:
            print(item)

# HEATMAP DISPLAY FUNCTIONS
    def heatmaps_nw_raw(self):
        """Display heatmaps for Needleman-Wunsch Scoring and Traceback Matrices"""
        fig = px.imshow(self.nw_s_matrix,
                        labels=dict(x="Sequence 1", y="Sequence 2", color = "Productivity"),
                        y=[""]+[f"{char}-{i}" for i, char in enumerate(self.seq_1)],
                        x=[""]+[f"{char}-{i}" for i, char in enumerate(self.seq_2)],
                        text_auto=True,
                        color_continuous_scale='Viridis')
        fig.update_layout(title="Needleman-Wunsch Scoring Matrix Heatmap")
        fig.update_xaxes(side="top")
        fig.show()

        fig_2 = px.imshow(self.nw_t_matrix,
                          labels=dict(x="Sequence 1", y="Sequence 2", color="Productivity"),
                          y=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_1)],
                          x=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_2)],
                          text_auto=True,
                          color_continuous_scale='Viridis')
        fig_2.update_layout(title="Needleman-Wunsch Traceback Matrix Heatmap")
        fig_2.update_xaxes(side="top")
        fig_2.show()

    def heatmaps_sw_raw(self):
        """Display heatmaps for Smith-Waterman Scoring and Traceback Matrices"""
        fig = px.imshow(self.sw_s_matrix,
                        labels=dict(x="Sequence 1", y="Sequence 2", color="Productivity"),
                        y=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_1)],
                        x=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_2)],
                        text_auto=True,
                        color_continuous_scale='Viridis')
        fig.update_layout(title="Smith-Waterman Scoring Matrix Heatmap")
        fig.update_xaxes(side="top")
        fig.show()

        fig_2 = px.imshow(self.sw_t_matrix,
                          labels=dict(x="Sequence 1", y="Sequence 2", color="Productivity"),
                          y=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_1)],
                          x=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_2)],
                          text_auto=True,
                          color_continuous_scale='Viridis')
        fig_2.update_layout(title="Smith-Waterman Traceback Matrix Heatmap")
        fig_2.update_xaxes(side="top")
        fig_2.show()

    def heatmap_nw_path(self):
        """Display heatmap for Needleman-Wunsch Path Matrix"""
        fig = px.imshow(self.nw_path_matrix,
                        labels=dict(x="Sequence 1", y="Sequence 2", color="Productivity"),
                        y=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_1)],
                        x=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_2)],
                        text_auto=True,
                        color_continuous_scale='thermal')
        fig.update_layout(title="Needleman-Wunsch Path Matrix Heatmap")
        fig.update_xaxes(side="top")
        fig.show()

    def heatmap_sw_path(self):
        """Display heatmap for Smith-Waterman Path Matrix"""
        fig = px.imshow(self.sw_path_matrix,
                        labels=dict(x="Sequence 1", y="Sequence 2", color="Productivity"),
                        y=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_1)],
                        x=[""] + [f"{char}-{i}" for i, char in enumerate(self.seq_2)],
                        text_auto=True,
                        color_continuous_scale='thermal')
        fig.update_layout(title="Smith-Waterman Path Matrix Heatmap")
        fig.update_xaxes(side="top")
        fig.show()
