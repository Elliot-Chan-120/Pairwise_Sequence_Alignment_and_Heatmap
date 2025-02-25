class SubstitutionMatrix:
    """
    Generates substitution matrices depending on which biotype is chosen
    """
    def __init__(self, seq_type):
        self.seq_type = seq_type.upper()
        self.sub_matrix = {}

    def auto_load(self):
        """
        Detects what the seq_type is and outputs its default substitution matrix and parameters \n
        If you want different settings/protein matrix, use the other functions and specify your parameters
        \nProtein -> blosum62 matrix
        \nDNA -> match = 1, mismatch = 0
        \nRNA -> match = 1, mismatch = 0
        """
        if self.seq_type == "PROTEIN":
            self.load_protein_matrix("blosum62.mat")
        elif self.seq_type == "DNA":
            self.DNA_submat(1, 0)
        elif self.seq_type == "RNA":
            self.RNA_submat(1, 0)
        else:
            raise ValueError("Invalid Sequence Type - Choose from: Protein, DNA, or RNA")
        return self.sub_matrix

    def load_protein_matrix(self, filename="blosum62.mat"):
        """
        When calling any matrix (default is blosum62), will parse the .mat file into a substitution matrix
        \n substitution matrix is dict containing scores b/w diff. amino acids
        """
        if self.seq_type != "PROTEIN": # error handling in case you try to use a protein sub matrix for R/DNA
            raise ValueError("File loading is for protein sequences")
        with open(filename, 'r') as f:
            header = f.readline().strip().split()
            alphabet = header

            for line in f:
                tokens = line.strip().split()
                row_aa = tokens[0]
                scores = tokens[1:]

                for col_aa, score in zip(alphabet, scores):
                    pair = row_aa + col_aa
                    self.sub_matrix[pair] = int(score)
        return self.sub_matrix

    def DNA_submat(self, match, mismatch, alphabet="ATCG"):
        """
        Generates DNA substitution matrix
        :param match: score on match
        :param mismatch: penalty on mismatch
        :param alphabet: list of allowed characters
        :return: substitution matrix for DNA
        """
        for char_1 in alphabet:
            for char_2 in alphabet:
                if char_1 == char_2: # if it matches - score set to match value and vice versa
                    self.sub_matrix[char_1 + char_2] = match
                else:
                    self.sub_matrix[char_1 + char_2] = mismatch
        return self.sub_matrix


    def RNA_submat(self, match, mismatch, alphabet="AUCG"):
        """
        Generates RNA substitution matrix
        :param match: score on match
        :param mismatch: penalty on mismatch
        :param alphabet: list of allowed characters
        :return: substitution matrix for RNA
        """
        for char_1 in alphabet:
            for char_2 in alphabet:
                if char_1 == char_2:
                    self.sub_matrix[char_1 + char_2] = match
                else:
                    self.sub_matrix[char_1 + char_2] = mismatch
        return self.sub_matrix



def test_output():
    """
    Use this is you want to take a look at the matrix \n
    DO NOT RUN WITH main.py - THERE IS A DELIBERATE VALUEERROR TEST MATRIX
    :return: all substitution matrices (raw, not printed nicely)
    """
    print(SubstitutionMatrix("protein").auto_load())
    print(SubstitutionMatrix("dna").auto_load())
    print(SubstitutionMatrix("rna").auto_load())
    print(SubstitutionMatrix("error_test").auto_load())
# test_output()
