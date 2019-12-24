# 1 : making seqence object

from Bio.Seq import Seq

tatabox_seq = Seq("tataaaggcAATATGCAGTAG")
print(type(tatabox_seq))   # <class 'Bio.Seq.Seq'>

""" print(dir(Seq))
['__add__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', 
 '__ge__','__getattribute__', '__getitem__', '__gt__', '__hash__', '__imul__', '__init__', '__init_subclass__', 
 '__le__', '__len__', '__lt__', '__module__', '__mul__', '__ne__', '__new__', '__radd__', '__reduce__', 
 '__reduce_ex__', '__repr__', '__rmul__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 
 '_get_seq_str_and_check_alphabet', 'back_transcribe', 'complement', 'count', 'count_overlap', 'endswith', 'find', 
 'index', 'join', 'lower', 'lstrip', 'reverse_complement', 'rfind', 'rindex', 'rsplit', 'rstrip', 'split', 
 'startswith', 'strip', 'tomutable', 'transcribe', 'translate', 'ungap', 'upper']"""


# 2 : adding information ~ what kind seq is this?(DNA/RNA/amino_acid)

# from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

tatabox_seq = Seq("tataaaggcAATATGCAGTAG", IUPAC.unambiguous_dna)
    # from now on , tata seq is a DNA seq


# 3 : information from simple seqence

# from Bio.Seq import Seq

exon_seq = Seq("ATGCAGTAG")
count_a = exon_seq.count("A")
print(count_a) # 3

    # exon_seq = Seq("ATGCAGTAG")
g_count = exon_seq.count("G")
c_count = exon_seq.count("C")
gc_contents = ((g_count + c_count) / 9) * 100  # 9 = len(exon_seq)
print(gc_contents) # 44.44444444444444
    # 또는
# from Bio.Seq import Seq
from Bio.SeqUtils import GC
# exon_seq = Seq("ATGCAGTAG")
print(  GC(exon_seq)  ) # 44.44444444444444

    # tatabox_seq = Seq("tataaaggcAATATGCAGTAG")
print(tatabox_seq.upper()) # TATAAAGGCAATATGCAGTAG
print(tatabox_seq.lower()) # tataaaggcaatatgcagtag

dna = Seq("AAAAAAGTCGCAGAAAATATATGAGGAAACAAAAAGCGAAGACGACAAAAAAAAAAAAAACT") # front part(62) of PHYB gene
mrna = (dna.transcribe())
ptn = mrna.translate()
print(mrna) # AAAAAAGUCGCAGAAAAUAUAUGAGGAAACAAAAAGCGAAGACGACAAAAAAAAAAAAAACU
print(ptn) # KKVAENI*GNKKRRRQKKKK
        # * = stop codon
ptn2 = mrna.translate(to_stop = True)
print(ptn2) # KKVAENI

for seq in ptn.split("*"):
    print(seq) # KKVAENI
               # GNKKRRRQKKKK
        # 실험
    """d = ptn.split("*")
    print(d) # [Seq('KKVAENI', HasStopCodon(ExtendedIUPACProtein(), '*')),
    # Seq('GNKKRRRQKKKK', HasStopCodon(ExtendedIUPACProtein(), '*'))]"""


# 4 : making reverse complement seqence

    # without biopython
        #dna2 = "CGAGAAAAGGAAAAAAAAAAATAGAAAGAGAAAACGCTTAGTATCTCCGGCGACTT" # front part(56) of FLC gene
        #comp_dic = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        #comp_seq = ""
        #for s in dna2:
        #    comp_seq += comp_dic[s]
        #revcomp_seq = comp_seq[::-1]
        #print(comp_seq)
        #print(revcomp_seq)
    # with biopython
dna2 = Seq("CGAGAAAAGGAAAAAAAAAAATAGAAAGAGAAAACGCTTAGTATCTCCGGCGACTT")
comp_dna2 = dna2.complement()
rev_comp_dna2 = dna2.reverse_complement()
print(comp_dna2) # GCTCTTTTCCTTTTTTTTTTTATCTTTCTCTTTTGCGAATCATAGAGGCCGCTGAA
print(rev_comp_dna2) # AAGTCGCCGGAGATACTAAGCGTTTTCTCTTTCTATTTTTTTTTTTCCTTTTCTCG


# 5 : codon table

from Bio.Data import CodonTable

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
print(codon_table)
"""
Table 1 Standard, SGC0

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I   | ACT T   | AAT N   | AGT S   | T
A | ATC I   | ACC T   | AAC N   | AGC S   | C
A | ATA I   | ACA T   | AAA K   | AGA R   | A
A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V   | GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--
"""
codon_table2 = CodonTable.unambiguous_dna_by_name["Plant Plastid"]
#print(codon_table2)
"""
Table 11 Bacterial, Archaeal, Plant Plastid

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
A | ATA I(s)| ACA T   | AAA K   | AGA R   | A
A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--
"""


# 6 : find ORF from sequence object

# from Bio.Seq import Seq
PHYB = Seq("AAAAAAGTCGCAGAAAATATATGAGGAAACAAAAAGCGAAGACGACAAAAAAAAAAAAAACTCTGATTTTTTTTTGTTATCTCTCTCTATCTGAGAGGCACACATTT\
TGCTTCGTCTTCTTCAATTTATTTTATTGGTTTCTCCACTTATCTCCGATCTCAATTCTCCCCATTTTCTTCTTCCTCAAGTTCAAAATTCTTGAGAATTTAGCTCTACCAGAATTCGT\
CTCCGATAACTAGTGGATGATGATTCACCCTAAATCCTTCCTTGTCTCAAGGTAATTCTGAGAAATTTCTCAAATTCAAAATCAAACGGCATGGTTTCCGGAGTCGGGGGTAGTGGCGG\
TGGCCGTGGCGGTGGCCGTGGCGGAGAAGAAGAACCGTCGTCAAGTCACACTCCTAATAACCGAAGAGGAGGAGAACAAGCTCAATCGTCGGGAACGAAATCTCTCAGACCAAGAAGCA\
ACACTGAATCAATGAGCAAAGCAATTCAACAGTACACCGTCGACGCAAGACTCCACGCCGTTTTCGAACAATCCGGCGAATCAGGGAAATCATTCGACTACTCACAATCACTCAAAAC")
    # 582 base from PHYB
start_idx = PHYB.find("ATG")
end_idx = PHYB.find("TAG", start_idx) # skipped TAA, TGA for convenience
orf = PHYB[start_idx : end_idx+3]
print(orf)
    # ATGAGGAAACAAAAAGCGAAGACGACAAAAAAAAAAAAAACTCTGATTTTTTTTTGTTATCTCTCTCTATCTGAGAGGCACACATTTTGCTTCGTCTTCTTCAATTTATTTT
    # ATTGGTTTCTCCACTTATCTCCGATCTCAATTCTCCCCATTTTCTTCTTCCTCAAGTTCAAAATTCTTGAGAATTTAG

# 7 : Weight of sequence

# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
from Bio.SeqUtils import molecular_weight

seq1 = Seq("ATGCAGTAG")
seq2 = Seq("ATGCAGTAG", IUPAC.unambiguous_dna)
seq3 = Seq("ATGCAGTAG", IUPAC.protein)
print(molecular_weight(seq1)) # 2842.8206999999993
print(molecular_weight(seq2)) # 2842.8206999999993
print(molecular_weight(seq3)) # 707.7536


# 8 : get all translation variation from  DNA sequence (가능한 모든 번역 형태 구하기)

# from Bio.Seq import Seq
from Bio.SeqUtils import six_frame_translations

seq4 = Seq("ATGCCTTGAAATGTGTAG")
print(six_frame_translations(seq4))
"""
1/1
  A  L  K  C  V
 C  L  E  M  C
M  P  *  N  V  *
atgccttgaaatgtgtag   38 %
tacggaactttacacatc
G  Q  F  T  Y 
 H  R  S  I  H  L
  A  K  F  H  T
"""

# 9 :  calculating Tm

# from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

myseq = Seq("AGTCTGGGACGCGGCAATCGCA")
print(mt.Tm_Wallace(myseq)) # 72.0


# 10 : amino acid 1 to 3

from Bio.SeqUtils import seq1

essential_amino_acid_3 = "LeuLysMetValIleThrTrpPhe"
print(seq1(essential_amino_acid_3)) # LKMVITWF

from Bio.SeqUtils import seq3

essential_amino_acid_1 = "LKMVITWF"
print(seq3(essential_amino_acid_1)) # LeuLysMetValIleThrTrpPhe
