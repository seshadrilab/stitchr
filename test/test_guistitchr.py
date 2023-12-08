import filecmp
from Bio import SeqIO
from Stitchr import stitchrfunctions as fxn

class TestOutFiles:
    """Ensures exported FASTAs contain expected information. This is our 
    'end-to-end' test.

    To test, run the GUI with the inputs below, exporting to output FASTAs.
    Ensure Link chains=P2A and Link order=BA or GD.

    # WITH BASELINE EXAMPLE DATA (INCLUDING BASELINE CONSTANTS)
    test_alphabeta
    With TRA/TRB, Link chains=P2A, Link order=BA:
        1. 'Example data'
        2. TRAC gene name=TRAC*01
        3. TRBC gene name=TRBC1*01
        4. 'Run Stitchr'
        5. 'Export output': stitchr/test/new_ab.fasta

    test_gammadelta
    With TRG/TRD, Link chains=P2A, Link order=GD:
        1. 'Example data'
        2. 'Run Stitchr'
        3. 'Export output': stitchr/test/new_gd.fasta

    # WITH BASELINE EXAMPLE DATA PLUS NEW DEFAULT CONSTANT (MOUSE FOR ALPHA-BETA)
    test_alphabeta_mouse_c
    With TRA/TRB, Link chains=P2A, Link order=BA:
        1. 'Example data'
        2. 'Run Stitchr'
        3. 'Export output': stitchr/test/mouse_ab.fasta

    # WITH BASELINE EXAMPLE DATA AND NEW DEFAULTS AND RESTRICTION SITES ADDED
    test_alphabeta_restriction_sites
    With TRA/TRB, Link chains=P2A, Link order=BA:
        1. 'Example data'
        2. Check 'Add Restriction Sites (BamHI, SalI)'
        3. 'Run Stitchr'
        4. 'Export output': stitchr/test/restrict_ab.fasta

    test_gammadelta_restriction_sites
    With TRG/TRD, Link chains=P2A, Link order=GD:
        1. 'Example data'
        2. Check 'Add Restriction Sites (BamHI, SalI)'
        3. 'Run Stitchr'
        4. 'Export output': stitchr/test/restrict_gd.fasta
    """

    def test_alphabeta(self):
        # Does not contain restriction sites so files should be identical
        assert filecmp.cmp('baseline_ab.fasta', 'new_ab.fasta')


    def test_gammadelta(self):
        # Gamma chain contains two BamHI restriction sites
        base_nt = []
        base_aa = []
        new_nt = []
        new_aa = []

        for record in SeqIO.parse('baseline_gd.fasta', "fasta"):
            seq = record.seq
            base_nt.append(seq)
            base_aa.append(fxn.translate_nt(str(seq)))
        
        for record in SeqIO.parse('new_gd.fasta', "fasta"):
            seq = record.seq
            new_nt.append(seq)
            new_aa.append(fxn.translate_nt(str(seq)))

        # Confirm only the first record (gamma chain) and final output (both chains)
        # have changed nucleotide sequences
        assert base_nt[0] != new_nt[0]
        assert base_nt[1] == new_nt[1]
        assert base_nt[2] != new_nt[2]

        # Confirm all translations the same
        assert base_aa[0] == new_aa[0]
        assert base_aa[1] == new_aa[1]
        assert base_aa[2] == new_aa[2]


    def test_alphabeta_mouse_c(self):
        mouse_a_nt = 'gacattcagaacccggaaccggctgtataccagctgaaggacccccgatctcaggatagtactctgtgcctgttcaccgactttgatagtcagatcaatgtgcctaaaaccatggaatccggaacttttattaccgacaagtgcgtgctggatatgaaagccatggacagtaagtcaaacggcgccatcgcttggagcaatcagacatccttcacttgccaggatatcttcaaggagaccaacgcaacatacccatcctctgacgtgccctgtgatgccaccctgacagagaagtctttcgaaacagacatgaacctgaattttcagaatctgagcgtgatgggcctgagaatcctgctgctgaaggtcgctgggtttaatctgctgatgacactgcggctgtggtcctcatga'.upper()
        mouse_b_nt = 'gaagatctacgtaacgtgacaccacccaaagtctcactgtttgagcctagcaaggcagaaattgccaacaagcagaaggccaccctggtgtgcctggcaagagggttctttccagatcacgtggagctgtcctggtgggtcaacggcaaagaagtgcattctggggtctgcaccgacccccaggcttacaaggagagtaattactcatattgtctgtcaagccggctgagagtgtccgccacattctggcacaaccctaggaatcatttccgctgccaggtccagtttcacggcctgagtgaggaagataaatggccagaggggtcacctaagccagtgacacagaacatcagcgcagaagcctggggacgagcagactgtggcattactagcgcctcctatcatcagggcgtgctgagcgccactatcctgtacgagattctgctgggaaaggccaccctgtatgctgtgctggtctccggcctggtgctgatggccatggtcaagaaaaagaactct'.upper()
        mouse_stitched_nt = 'ATGAGAATCAGGCTCCTGTGCTGTGTGGCCTTTTCTCTCCTGTGGGCAGGTCCAGTGATTGCTGGGATCACCCAGGCACCAACATCTCAGATCCTGGCAGCAGGACGGCGCATGACACTGAGATGTACCCAGGATATGAGACATAATGCCATGTACTGGTATAGACAAGATCTAGGACTGGGGCTAAGGCTCATCCATTATTCAAATACTGCAGGTACCACTGGCAAAGGAGAAGTCCCTGATGGTTATAGTGTCTCCAGAGCAAACACAGATGATTTCCCCCTCACGTTGGCGTCTGCTGTACCCTCTCAGACATCTGTGTACTTCTGTGCCAGCAGTCTGAGCTTCGGCACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAGAAGATCTACGTAACGTGACACCACCCAAAGTCTCACTGTTTGAGCCTAGCAAGGCAGAAATTGCCAACAAGCAGAAGGCCACCCTGGTGTGCCTGGCAAGAGGGTTCTTTCCAGATCACGTGGAGCTGTCCTGGTGGGTCAACGGCAAAGAAGTGCATTCTGGGGTCTGCACCGACCCCCAGGCTTACAAGGAGAGTAATTACTCATATTGTCTGTCAAGCCGGCTGAGAGTGTCCGCCACATTCTGGCACAACCCTAGGAATCATTTCCGCTGCCAGGTCCAGTTTCACGGCCTGAGTGAGGAAGATAAATGGCCAGAGGGGTCACCTAAGCCAGTGACACAGAACATCAGCGCAGAAGCCTGGGGACGAGCAGACTGTGGCATTACTAGCGCCTCCTATCATCAGGGCGTGCTGAGCGCCACTATCCTGTACGAGATTCTGCTGGGAAAGGCCACCCTGTATGCTGTGCTGGTCTCCGGCCTGGTGCTGATGGCCATGGTCAAGAAAAAGAACTCTGGCAGCGGCGCCACCAACTTCAGCCTGCTGAAGCAGGCCGGCGACGTGGAGGAGAACCCCGGCCCCATGAAATCCTTGAGAGTTTTACTAGTGATCCTGTGGCTTCAGTTGAGCTGGGTTTGGAGCCAACAGAAGGAGGTGGAGCAGAATTCTGGACCCCTCAGTGTTCCAGAGGGAGCCATTGCCTCTCTCAACTGCACTTACAGTGACCGAGGTTCCCAGTCCTTCTTCTGGTACAGACAATATTCTGGGAAAAGCCCTGAGTTGATAATGTTCATATACTCCAATGGTGACAAAGAAGATGGAAGGTTTACAGCACAGCTCAATAAAGCCAGCCAGTATGTTTCTCTGCTCATCAGAGACTCCCAGCCCAGTGATTCAGCCACCTACCTCTGTGCCGTGAACTTCGGCGGAGGAAAGCTTATCTTCGGACAGGGAACGGAGTTATCTGTGAAACCCGACATTCAGAACCCGGAACCGGCTGTATACCAGCTGAAGGACCCCCGATCTCAGGATAGTACTCTGTGCCTGTTCACCGACTTTGATAGTCAGATCAATGTGCCTAAAACCATGGAATCCGGAACTTTTATTACCGACAAGTGCGTGCTGGATATGAAAGCCATGGACAGTAAGTCAAACGGCGCCATCGCTTGGAGCAATCAGACATCCTTCACTTGCCAGGATATCTTCAAGGAGACCAACGCAACATACCCATCCTCTGACGTGCCCTGTGATGCCACCCTGACAGAGAAGTCTTTCGAAACAGACATGAACCTGAATTTTCAGAATCTGAGCGTGATGGGCCTGAGAATCCTGCTGCTGAAGGTCGCTGGGTTTAATCTGCTGATGACACTGCGGCTGTGGTCCTCATGA'
        
        mouse_a_aa = 'DIQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKCVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRLWSS*'
        mouse_b_aa = 'EDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVCTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNS'
        mouse_stitched_aa = 'MRIRLLCCVAFSLLWAGPVIAGITQAPTSQILAAGRRMTLRCTQDMRHNAMYWYRQDLGLGLRLIHYSNTAGTTGKGEVPDGYSVSRANTDDFPLTLASAVPSQTSVYFCASSLSFGTEAFFGQGTRLTVVEDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVCTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNSGSGATNFSLLKQAGDVEENPGPMKSLRVLLVILWLQLSWVWSQQKEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMFIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVNFGGGKLIFGQGTELSVKPDIQNPEPAVYQLKDPRSQDSTLCLFTDFDSQINVPKTMESGTFITDKCVLDMKAMDSKSNGAIAWSNQTSFTCQDIFKETNATYPSSDVPCDATLTEKSFETDMNLNFQNLSVMGLRILLLKVAGFNLLMTLRLWSS*'
        
        # Confirm sequences are as expected based on manual checking
        base_nt = []
        base_aa = []
        mouse_nt = []
        mouse_aa = []

        for record in SeqIO.parse('baseline_ab.fasta', "fasta"):
            seq = record.seq
            base_nt.append(str(seq))
            base_aa.append(fxn.translate_nt(str(seq)))

        for record in SeqIO.parse('mouse_ab.fasta', "fasta"):
            seq = record.seq
            mouse_nt.append(str(seq))
            mouse_aa.append(fxn.translate_nt(str(seq)))

        # Nucleotide: Ensure sequences are the same up until the constant regions
        assert base_nt[0][:393] == mouse_nt[0][:393]
        assert base_nt[1][:393] == mouse_nt[1][:393]
        # Nucleotide: Ensure the mouse constant region is present in new output
        assert mouse_nt[0][393:] == mouse_a_nt
        assert mouse_nt[1][393:] == mouse_b_nt
        # Nucleotide: Ensure stitched output is as expected
        assert mouse_nt[2] == mouse_stitched_nt

        # Amino acid: Ensure sequences are the same up until the constant regions
        assert base_aa[0][:131] == mouse_aa[0][:131]
        assert base_aa[1][:131] == mouse_aa[1][:131]
        # Amino acid: Ensure the mouse constant region is present in new output
        assert mouse_aa[0][131:] == mouse_a_aa
        assert mouse_aa[1][131:] == mouse_b_aa

        # Amino acid: Ensure stitched output is as expected
        assert mouse_aa[2] == mouse_stitched_aa


    def test_alphabeta_restriction_sites(self):
        # Confirm same as mouse_ab except for restriction sites at the ends
        mouse_nt = []
        mouse_aa = []
        restrict_nt = []
        restrict_aa = []

        for record in SeqIO.parse('mouse_ab.fasta', "fasta"):
            seq = record.seq
            mouse_nt.append(str(seq))
            mouse_aa.append(fxn.translate_nt(str(seq)))

        for record in SeqIO.parse('restrict_ab.fasta', "fasta"):
            seq = record.seq
            restrict_nt.append(str(seq))
            restrict_aa.append(fxn.translate_nt(str(seq)))

        # Nucleotide: Ensure sequences the same except for restriction sites
        assert mouse_nt[2] == restrict_nt[2][6:-6]
        # Nucleotide: Ensure restriction sites added
        assert restrict_nt[2][:6] == 'GGATCC'
        assert restrict_nt[2][-6:] == 'GTCGAC'

        # Amino acid: Ensure sequences the same except for restriction sites
        assert mouse_aa[2] == restrict_aa[2][2:-2]
        # Amino acid: Ensure restriction sites added
        assert restrict_aa[2][:2] == 'GS'
        assert restrict_aa[2][-2:] == 'VD'


    def test_gammadelta_restriction_sites(self):
        # Confirm same as new_gd except for restriction sites at the ends
        new_nt = []
        new_aa = []
        restrict_nt = []
        restrict_aa = []

        for record in SeqIO.parse('new_gd.fasta', "fasta"):
            seq = record.seq
            new_nt.append(seq)
            new_aa.append(fxn.translate_nt(str(seq)))

        for record in SeqIO.parse('restrict_gd.fasta', "fasta"):
            seq = record.seq
            restrict_nt.append(str(seq))
            restrict_aa.append(fxn.translate_nt(str(seq)))

        # Nucleotide: Ensure sequences the same except for restriction sites
        assert new_nt[2] == restrict_nt[2][6:-6]
        # Nucleotide: Ensure restriction sites added
        assert restrict_nt[2][:6] == 'GGATCC'
        assert restrict_nt[2][-6:] == 'GTCGAC'

        # Amino acid: Ensure sequences the same except for restriction sites
        assert new_aa[2] == restrict_aa[2][2:-2]
        # Amino acid: Ensure restriction sites added
        assert restrict_aa[2][:2] == 'GS'
        assert restrict_aa[2][-2:] == 'VD'
