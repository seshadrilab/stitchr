from Stitchr import stitchrfunctions as fxn



class TestRestrict:
    """
    Ensures that the output of the restriction functions is the same as input.  No preparation needed.
    """
    def process(s, enzymes):
        sites = fxn.check_restricts(s, enzymes)
        sequence = fxn.wobble(s, sites, enzymes)
        return sequence
    
    def test_noSites(self):
        """
        Checks to see if two sequences are identical w/ and w/out editing with our process and no restriction sites
        """
        s = "ATGAAATCCTTGAGAGTTTTACTAGTGATCCTGTGGCTTCAGTTGAGCTGGGTTTGGAGCCAACAGAAGGAGGTGGAGCAGAATTCTGGACCCCTCAGTGTTCCAGAGGGAGCCATTGCCTCTCTCAACTGCACTTACAGTGACCGAGGTTCCCAGTCCTTCTTCTGGTACAGACAATATTCTGGGAAAAGCCCTGAGTTGATAATGTTCATATACTCCAATGGTGACAAAGAAGATGGAAGGTTTACAGCACAGCTCAATAAAGCCAGCCAGTATGTTTCTCTGCTCATCAGAGACTCCCAGCCCAGTGATTCAGCCACCTACCTCTGTGCCGTGAACTTCGGCGGAGGAAAGCTTATCTTCGGACAGGGAACGGAGTTATCTGTGAAACCCAATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGC"
        enzymes = ["BamHI", "SalI"]
        altered = TestRestrict.process(s, enzymes)
        assert altered == s
        altered = fxn.translate_nt(altered)
        unaltered = fxn.translate_nt(s)
        assert altered == unaltered

    def test_oneSite(self):
        """
        Checks to see if two sequences are identical w/ and w/out editing with our process and one restriction site
        """
        s = "ATGAAATCCTTGAGAGTTTTACTAGTGATCCTGTGGCTTCAGTTGAGCTGGGTTTGGAGCCAACAGAAGGAGGTGGAGCAGAATTCTGGACCCCTCAGTGTTCCAGAGGGAGCCATTGCCTCTCTCAACTGCACTTACAGTGACCGAGGTTCCCAGTCCTTCTTCTGGTACAGACAATATTCTGGGAAAAGCCCTGAGTTGATAATGTTCATATACTCCAATGGTGACAAAGAAGATGGAAGGTTTACAGCACAGCTCAATAAAGCCAGCCAGTATGTTTCTCTGCTCATCAGAGACTCCCAGCCCAGTGATTCAGCCACCTACCTCTGTGCCGTGAACTTCGGCGGAGGAAAGCTTATCTTCGGACAGGGAACGGAGTTATCTGTGAAACCCAATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGC"
        enzymes = ["EcoRI"]
        altered = TestRestrict.process(s, enzymes)
        assert altered != s
        altered = fxn.translate_nt(altered)
        unaltered = fxn.translate_nt(s)
        assert altered == unaltered

    #TODO tests for enzymes that have * or N in their restriction site
    #TODO tests for different length enzyme sites
    #TODO tests for different reading frames for sites
    #TODO test for an enzyme site at beginning of sequence
    #TODO test for an enzyme site at end of sequence
    #TODO test for codons that have only one nt coding option