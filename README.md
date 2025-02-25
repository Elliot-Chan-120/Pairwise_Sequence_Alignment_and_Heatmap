# Pairwise_Sequence_Alignment_and_Heatmap
 Using the Needleman-Wunsch and Smith-Waterman algorithms, this project uses two classes to align two sequences of the same biotype (protein, R/DNA): a Substitution Matrix Generator (Sub_Matrix_Gen) and Alignment class (Sequence_Align_Toolkit), which performs the actual alignment and data processing. We can obtain the optimal alignment and the score of the two sequences using either the global (NW) or local (SW) algorithms or both simultaneously. The alignment, scores and matrices can be output neatly in your editor. All produced matrices (Scoring, Traceback, Alignment 'Travel' Path) can be visually displayed as a heatmap in this.
 
 You can import whatever protein substitution matrix such as PAM (I used blosum62), and you can edit your own gap penalties for all biotype matrices.

 Requires plotly for the heatmap visualization.
