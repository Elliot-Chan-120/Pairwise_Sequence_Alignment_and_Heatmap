from Sub_Matrix_Gen import SubstitutionMatrix
from Sequence_Align_Toolkit import SequenceAlign

# NEEDLEMAN WUNSCH IS GLOBAL ALIGNMENT
# SMITH WATERMAN IS LOCAL ALIGNMENT

# create instance of our class by choosing which two sequences we want to align, their biotype and gap penalty
global_test = SequenceAlign("HMPlGWAG", "PHAAAAADIDWG", SubstitutionMatrix("proTeIn").auto_load(), -8)

# First, generate alignment data after instantiating our class
global_test.gen_alignment_data()

# output global alignment data
global_test.output_nw_results()
global_test.output_nw_raw_matrices()
global_test.output_nw_path()

# output local alignment data
global_test.output_sw_results()
global_test.output_sw_raw_matrices()
global_test.output_sw_path()

# # output heatmaps
global_test.heatmaps_nw_raw() # 2 heatmaps
global_test.heatmaps_sw_raw() # 2 heatmaps

global_test.heatmap_nw_path()
global_test.heatmap_sw_path()
