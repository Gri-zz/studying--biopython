# 1 : sample SeqRecord

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

simple_seq = Seq("AGTC")
seqRecord = SeqRecord(simple_seq)
print(seqRecord)
""""
ID: <unknown id>
Name: <unknown name>
Description: <unknown description>
Number of features: 0
Seq('AGTC')
"""
simple_seqRecord = SeqRecord(simple_seq, id="NC_111", name="Test", description="dYEEscription")
print(simple_seqRecord)
"""
D: NC_111
Name: Test
Description: <unknown description>
Number of features: 0
Seq('AGTC')
"""
simple_seqRecord.name = "Another name"
print(simple_seqRecord)
"""
Name: Another name
"""


# 2 : making SeqRecord object

# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

    # 2-1 : directly
seq2 = Seq("ACGT")
seqRecord2 = SeqRecord(seq2)
print(seqRecord2)
print("------------")
seqRecord2.id = "NC_1111"
seqRecord2.name = "GeneA"
seqRecord2.description = "This is DD"
seqRecord2.annotations["Annotation_Key1"] = "Annotation_Value1"
seqRecord2.annotations["Annotation_Key2"] = "Annotation_Value2"
print(seqRecord2)
"""
ID: <unknown id>
Name: <unknown name>
Description: <unknown description>
Number of features: 0
Seq('ACGT')
------------
ID: NC_1111
Name: GeneA
Description: This is DD
Number of features: 0
/Annotation_Key1=Annotation_Value1
/Annotation_Key2=Annotation_Value2
Seq('ACGT')
"""

    # 2-2 : from FASTA
from Bio import SeqIO

record = SeqIO.read("J01636.1.fasta","fasta")
print(type(record))
print(record)
"""
<class 'Bio.SeqRecord.SeqRecord'>
ID: J01636.1
Name: J01636.1
Description: J01636.1 E.coli lactose operon with lacI, lacZ, lacY and lacA genes
Number of features: 0
Seq('GACACCATCGAATGGCGCAAAACCTTTCGCGGTATGGCATGATAGCGCCCGGAA...GAC', SingleLetterAlphabet())
"""

    # 2-3 : from GeneBank
# from Bio import SeqIO
record2 = SeqIO.read("J01636.1.gbk","genbank")
print(type(record2))
#print(record2)
"""
<class 'Bio.SeqRecord.SeqRecord'>
ID: J01636.1
Name: ECOLAC
Description: E.coli lactose operon with lacI, lacZ, lacY and lacA genes
Number of features: 22
/molecule_type=DNA
/topology=linear
/data_file_division=BCT
/date=05-MAY-1993
/accessions=['J01636', 'J01637', 'K01483', 'K01793']
/sequence_version=1
.
.
.
"""


# 3 : comparing seqence objects

seq11 = Seq("ACGT")
record11 = SeqRecord(seq11)
seq22 = Seq("ACGT")
record22 = SeqRecord(seq22)
print(seq11 == seq22) # True
# print(record11 == record22) # NotImplementedError:...
# SeqRecord comparison is deliberately not implemented. Explicitly compare the attributes of interest.
    # = born to be error
print(record11.seq == record22.seq) # True # comparing attributes of SeqRecords
