# 1 : read FASTA with SeqIO.parse method

from Bio import SeqIO

seq = SeqIO.parse("sample_1.fasta","fasta")
print(type(seq))
for s in seq:
    print(type(s))
    print(s)
    print("") # for interval
"""<class 'generator'>
<class 'Bio.SeqRecord.SeqRecord'>
ID: AF501235.1
Name: AF501235.1
Description: AF501235.1 Influenzavirus A (A/duck/Shanghai/1/2000) hemagglutinin gene, complete cds
Number of features: 0
Seq('ATGGAGAAAATAGTGCTTCTTCTTGCAATAGTCAGTCTTGTTAAAAGTGATCAG...AGA', SingleLetterAlphabet())
"""
seq2 = SeqIO.parse("sample_2.fasta","fasta")
print(type(seq2))
for s in seq2:
    print(type(s))
    print(s)
    print("") # for interval
"""
<class 'generator'>
<class 'Bio.SeqRecord.SeqRecord'>
ID: MH464856.1
Name: MH464856.1
Description: MH464856.1 Hepatitis B virus isolate MA134, complete genome
Number of features: 0
Seq('TTCCACAACATTCCACCAAGCTCTGCAGGATCCCAGAGTAAGAGGCCTGTATTT...GGG', SingleLetterAlphabet())

<class 'Bio.SeqRecord.SeqRecord'>
ID: CP002925.1
Name: CP002925.1
Description: CP002925.1 Streptococcus pseudopneumoniae IS7493, complete genome
Number of features: 0
Seq('TTGAAAGAAAAACAATTTTGGAATCGTATATTAGAATTTGCTCAAGAAAGACTG...ATC', SingleLetterAlphabet())
"""

# 2 : read FASTA with SeqIO.read method

seq3 = SeqIO.read("sample_1.fasta", "fasta")
print(type(seq3))
print(seq3)
"""
<class 'Bio.SeqRecord.SeqRecord'>
ID: AF501235.1
Name: AF501235.1
Description: AF501235.1 Influenzavirus A (A/duck/Shanghai/1/2000) hemagglutinin gene, complete cds
Number of features: 0
Seq('ATGGAGAAAATAGTGCTTCTTCTTGCAATAGTCAGTCTTGTTAAAAGTGATCAG...AGA', SingleLetterAlphabet())
"""
#seq4 = SeqIO.read("sample_2.fasta", "fasta")
#print(type(Seq4))
#print(seq4) # ValueError: More than one record found in handle


# 3 : read FASTQ with SeqIO.parse method

seq5 = SeqIO.parse("sample_1.fastq", "fastq")
print(type(seq5))
for s in seq5:
    print(type(s))
    print(s)
    print("") # for interval
"""
<class 'Bio.SeqRecord.SeqRecord'>
ID: SRR000982.35E745RJU01DLQBClength=53
Name: SRR000982.35E745RJU01DLQBClength=53
Description: SRR000982.35E745RJU01DLQBClength=53
Number of features: 0
Per letter annotation for: phred_quality
Seq('ATCTCTACCCAAAGATTAATGGGGATTGGTGTGATATACGGCTGAATTGTACC', SingleLetterAlphabet())
"""
seq6 = SeqIO.parse("sample_1.fastq", "fastq")
for s in seq6:
    print(s.seq)
# AAGGCACCATGCAGAGATGCAAGGCCCCTTTCTAAGCCCTAGACTTCTGGATGACACTTCTAGAAACACCCTGGGCCAGAAGTGAACCTGCTGCCTTGA
# AGGGAATAACTCAGATCTCTACCCAAAGATTAATGGGGATTGGTGTGATATACGGCTGAATTGTACC

import gzip
#from Bio import SeqIO

with gzip.open("sample_1.fastq.gz", "rt") as handle:
    seq7 = SeqIO.parse(handle, "fastq")
    for s in seq7:
        print(s.seq)
        print("")
"""
AAGGCACCATGCAGAGATGCAAGGCCCCTTTCTAAGCCCTAGACTTCTGGATGACACTTCTAGAAACACCCTGGGCCAGAAGTGAACCTGCTGCCTT
GAAGGGAATAACTCAGATCTCTACCCAAAGATTAATGGGGATTGGTGTGATATACGGCTGAATTGTACCAAGGCACCATGCAGAGATGCAAGGCCCC
TTTCTAAGCCCTAGACTTCTGGATGACACTTCTAGAAACACCCTGGGCCAGAAGTGAACCTGCTGCCTTGAAGGGAATAACTCAGATCTCTACCCAA
AGATTAATGGGGATTGGTGTGATATACGGCTGAATTGTACC
"""


# 4 : read  GenBank file

#from Bio import SeqIO

gbk = SeqIO.read("KT225476.2.gbk", "genbank")
print(type(gbk))
print(gbk)
"""
<class 'Bio.SeqRecord.SeqRecord'>
ID: KT225476.2
Name: KT225476
Description: Middle East respiratory syndrome coronavirus isolate MERS-CoV/THA/CU/17_06_2015, complete genome
Number of features: 12
/molecule_type=RNA
/topology=linear
/data_file_division=VRL
/date=22-AUG-2017
/accessions=['KT225476']
/sequence_version=2
/keywords=['']
/source=Middle East respiratory syndrome-related coronavirus (MERS-CoV)
/organism=Middle East respiratory syndrome-related coronavirus
/taxonomy=['Viruses', 'ssRNA viruses', 'ssRNA positive-strand viruses, no DNA stage', 'Nidovirales', 'Coronaviridae',\
 'Coronavirinae', 'Betacoronavirus']
/references=[Reference(title='Imported case of Middle East respiratory syndrome coronavirus (MERS-CoV) infection from\
 Oman to Thailand, June 2015', ...), Reference(title='Direct Submission', ...),\
  Reference(title='Direct Submission', ...)]
/comment=On Sep 10, 2015 this sequence version replaced KT225476.1.
/structured_comment=OrderedDict([('Assembly-Data', OrderedDict([('Sequencing Technology', \
'Sanger dideoxy sequencing')]))])
Seq('AGTGAATAGCTTGGCTATCTCACTTCCCCTCGTTCTCTTGCAGAACTTTGATTT...CTC', IUPACAmbiguousDNA())
"""
print("")

gbk2 = SeqIO.read("KT225476.2.gbk", "genbank")
print(gbk2.id)
print(gbk2.description)
print(gbk2.annotations['molecule_type'])
print(gbk2.annotations['organism'])
"""
KT225476.2
Middle East respiratory syndrome coronavirus isolate MERS-CoV/THA/CU/17_06_2015, complete genome
RNA
Middle East respiratory syndrome-related coronavirus
"""
print("")


# 5 : read from internet

from Bio import Entrez
#Entrez.email = "sample"

#handle = Entrez.efetch(db="nucleotide", rettype="gb", id="AY463215", retmode="text")

#for s in handle:
#    print(s.strip())
"""
LOCUS       AY463215                1059 bp    DNA     linear   PRI 14-DEC-2004
DEFINITION  Homo sapiens CCR5 chemokine receptor (CCR5) gene, complete cds.
ACCESSION   AY463215
VERSION     AY463215.1
KEYWORDS    .
SOURCE      Homo sapiens (human)
ORGANISM  Homo sapiens
Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
Catarrhini; Hominidae; Homo.
REFERENCE   1  (bases 1 to 1059)
AUTHORS   Capoulade-Metay,C., Ma,L., Truong,L.X., Dudoit,Y., Versmisse,P.,
Nguyen,N.V., Nguyen,M., Scott-Algara,D., Barre-Sinoussi,F.,
Debre,P., Bismuth,G., Pancino,G. and Theodorou,I.
TITLE     New CCR5 variants associated with reduced HIV coreceptor function
in southeast Asia
JOURNAL   AIDS 18 (17), 2243-2252 (2004)
PUBMED   15577536
REFERENCE   2  (bases 1 to 1059)
AUTHORS   Capoulade-Metay,C., Ma,L., Truong,L.X., Dudoit,Y., Versmisse,P.,
Nguyen,N.V., Nguyen,M., Scott-Algara,D., Barre-Sinoussi,F.,
Debre,P., Pancino,G. and Theodorou,I.
TITLE     Direct Submission
JOURNAL   Submitted (07-NOV-2003) INSERM U543, 83 Boulevard de l'Hopital,
Paris 75013, France
FEATURES             Location/Qualifiers
source          1..1059
/organism="Homo sapiens"
/mol_type="genomic DNA"
/db_xref="taxon:9606"
/country="Viet Nam"
gene            <1..>1059
/gene="CCR5"
mRNA            <1..>1059
/gene="CCR5"
/product="CCR5 chemokine receptor"
CDS             1..1059
/gene="CCR5"
/note="contains I254T variation"
/codon_start=1
/product="CCR5 chemokine receptor"
/protein_id="AAS19314.1"
/translation="MDYQVSSPIYDINYYTSEPCQKINVKQIAARLLPPLYSLVFIFG
FVGNMLVILILINCKRLKSMTDIYLLNLAISDLFFLLTVPFWAHYAAAQWDFGNTMCQ
LLTGLYFIGFFSGIFFIILLTIDRYLAVVHAVFALKARTVTFGVVTSVITWVVAVFAS
LPGIIFTRSQKEGLHYTCSSHFPYSQYQFWKNFQTLKIVILGLVLPLLVMVICYSGIL
KTLLRCRNEKKRHRAVRLIFTIMIVYFLFWAPYNTVLLLNTFQEFFGLNNCSSSNRLD
QAMQVTETLGMTHCCINPIIYAFVGEKFRNYLLVFFQKHIAKRFCKCCSIFQQEAPER
ASSVYTRSTGEQEISVGL"
variation       758
/gene="CCR5"
/replace="t"
ORIGIN
1 atggattatc aagtgtcaag tccaatctat gacatcaatt attatacatc ggagccctgc
61 caaaaaatca atgtgaagca aatcgcagcc cgcctcctgc ctccgctcta ctcactggtg
121 ttcatctttg gttttgtggg caacatgctg gtcatcctca tcctgataaa ctgcaaaagg
181 ctgaagagca tgactgacat ctacctgctc aacctggcca tctctgacct gtttttcctt
241 cttactgtcc ccttctgggc tcactatgct gccgcccagt gggactttgg aaatacaatg
301 tgtcaactct tgacagggct ctattttata ggcttcttct ctggaatctt cttcatcatc
361 ctcctgacaa tcgataggta cctggctgtc gtccatgctg tgtttgcttt aaaagccagg
421 acggtcacct ttggggtggt gacaagtgtg atcacttggg tggtggctgt gtttgcgtct
481 ctcccaggaa tcatctttac cagatctcaa aaagaaggtc ttcattacac ctgcagctct
541 cattttccat acagtcagta tcaattctgg aagaatttcc agacattaaa gatagtcatc
601 ttggggctgg tcctgccgct gcttgtcatg gtcatctgct actcgggaat cctaaaaact
661 ctgcttcggt gtcgaaatga gaagaagagg cacagggctg tgaggcttat cttcaccatc
721 atgattgttt attttctctt ctgggctccc tacaacactg tccttctcct gaacaccttc
781 caggaattct ttggcctgaa taattgcagt agctctaaca ggttggacca agctatgcag
841 gtgacagaga ctcttgggat gacgcactgc tgcatcaacc ccatcatcta tgcctttgtc
901 ggggagaagt tcagaaacta cctcttagtc ttcttccaaa agcacattgc caaacgcttc
961 tgcaaatgct gttctatttt ccagcaagag gctcccgagc gagcaagctc agtttacacc
1021 cgatccactg gggagcagga aatatctgtg ggcttgtga
//
"""

#from Bio import Entrez
#Entrez.email = "sample"
#from Bio import SeqIO

#with Entrez.efetch(db="nucleotide", rettype="fasta", remote="text", id="42540826")as handle:
#    seq9=SeqIO.read(handle, "fasta")
#print(seq9)
#print(len(seq9))
"""
ID: AY463215.1
Name: AY463215.1
Description: AY463215.1 Homo sapiens CCR5 chemokine receptor (CCR5) gene, complete cds
Number of features: 0
Seq('ATGGATTATCAAGTGTCAAGTCCAATCTATGACATCAATTATTATACATCGGAG...TGA', SingleLetterAlphabet())
1059
"""