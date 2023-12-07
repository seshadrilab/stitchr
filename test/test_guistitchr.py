import filecmp
from Bio import SeqIO
from Stitchr import stitchrfunctions as fxn

class TestOutFiles:
    """Ensures exported FASTAs are same as baseline versions.

    To test, run the GUI with the inputs below, exporting to output FASTAs.
    Ensure Link chains=P2A and Link order=BA or GD.

    test_alphabeta
    With TRA/TRB, Link chains=P2A, Link order=BA:
        1. 'Example data'
        2. 'Run Stitchr'
        3. 'Export output': stitchr/test/new_ab.fasta

    test_gammadelta
    With TRG/TRD, Link chains=P2A, Link order=GD:
        1. 'Example data'
        2. 'Run Stitchr'
        3. 'Export output': stitchr/test/new_gd.fasta

    test_template
        1. 'Find TCR input file': stitchr/templates/gui_input_example_human.tsv
        2. 'Upload TCR details'
        3. 'Run Stitchr'
        4. 'Export output': stitchr/test/new_template.fasta
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


    # TODO: Replace this test with Tori's TCR with known restriction site(s).
    # The extra sequence in the template is messing it up, and we plan to 
    # remove that functionality anyway

    # def test_template(self):
    #     # Beta chain contains one BamHI restriction site
    #     base_nt = []
    #     base_aa = []
    #     new_nt = []
    #     new_aa = []

    #     for record in SeqIO.parse('baseline_template.fasta', "fasta"):
    #         seq = record.seq
    #         base_nt.append(seq)
    #         base_aa.append(fxn.translate_nt(str(seq)))
        
    #     for record in SeqIO.parse('new_template.fasta', "fasta"):
    #         seq = record.seq
    #         new_nt.append(seq)
    #         new_aa.append(fxn.translate_nt(str(seq)))

    #     # Confirm only the second record (beta chain) and final output (both chains)
    #     # have changed nucleotide sequences
    #     assert base_nt[0] == new_nt[0]
    #     assert base_nt[1] != new_nt[1]
    #     assert base_nt[2] != new_nt[2]

    #     # Confirm all translations the same
    #     assert base_aa[0] == new_aa[0]
    #     assert base_aa[1] == new_aa[1]
    #     assert base_aa[2] == new_aa[2]
