from Sub_Matrix_Gen import SubstitutionMatrix
from Sequence_Align_Toolkit import SequenceAlign
from MultipleAlign import MultipleAlignment

# NEEDLEMAN WUNSCH IS GLOBAL ALIGNMENT
# SMITH WATERMAN IS LOCAL ALIGNMENT

# create instance of our class by choosing which two sequences we want to align, their biotype and gap penalty
pair_test = SequenceAlign(SubstitutionMatrix("proTeIn").auto_load(), -8)

def pairwise_alignment_and_heatmap_demonstration():
    """Run this to see how to use the SequenceAlign Class"""
    # create instance of our class by choosing which two sequences we want to align, their biotype and gap penalty
    pair_test = SequenceAlign(SubstitutionMatrix("proTeIn").auto_load(), -8)

    # First, generate alignment data after instantiating our class
    pair_test.set_sequences("HMPlGWAG", "PHAAAAADIDWG")
    pair_test.gen_all_data()

    # output global alignment data
    pair_test.output_nw_results()
    pair_test.output_nw_raw_matrices()
    pair_test.output_nw_path()

    # output local alignment data
    pair_test.output_sw_results()
    pair_test.output_sw_raw_matrices()
    pair_test.output_sw_path()

    # # output heatmaps
    pair_test.heatmaps_nw_raw() # 2 heatmaps
    pair_test.heatmaps_sw_raw() # 2 heatmaps

    pair_test.heatmap_nw_path()
    pair_test.heatmap_sw_path()

def multiple_sequence_demonstration():
    """Run this to see multiple sequence alignment"""
    s1 = "ATAGtccagctC"
    s2 = "AACggaattC"
    s3 = "ATGACtagcgat"
    submat = SubstitutionMatrix("DNA").auto_load()
    msa = MultipleAlignment([s1, s2, s3], submat, -8)
    result = msa.global_progressive_align()

    print("Optimal Multiple Sequence Alignment")
    for seq in result:
        print(seq)

    print("Consensus Sequence: ", msa.consensus())

multiple_sequence_demonstration()
