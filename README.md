# ASAP #

**Automated Sanger Analysis Pipeline (ASAP): a tool for rapidly analyzing Sanger sequencing data with minimum user interference.
(previously AMAP)**

## Introduction ##
ASAP is a code written in Python 3.x language. It incorporates various other programs to achieve its task.

ASAP accepts AB1 chromatogram files in pairs (Forward & Reverse), reference sequence as template in FASTA format and generates final alignment.

The program converts AB1 files into FASTQ, reverse complements reverse sequence, trims low quality ends of the sequences, aligns and creates a consensus sequence of these two forward and reverse sequence, aligns with the template sequence and displays the alignment. Also, the program saves all these files for future reference and study. Please note that the program will only display the alignment if not multiplexed, i.e. if only one pair of sequence is given to the program. In case of multiplexed runs, the program only saves the files which can be viewed easily in any alignment viewer program. Optionally, if provided with a FASTA file of exonic region and a reference amino acid file, ASAP also extracts the exonic region from the input sequence, translates it into +1, +2 & +3 frames and aligns it with the test reference amino acid sequence.

## Installation ##

ASAP doesn't require any installation.

Just make sure the system requirements (listed in the requirements section) are met.

Or, use one of the following **Docker** images:

```text
ARM (Apple M1/M2 etc.): saditya88/asap:arm
x86_64: saditya88/asap:latest
```
You may simply pull one of these images and run the program.

Example on an Apple M1 Max:

```shell
docker pull saditya88/asap:arm

docker run --rm -v /home/aditya/my_files/:/data/ saditya88/asap:arm FR /data/forward1.ab1 /data/reverse1.ab1 /data/reference.fasta
```

ASAP can be executed directly without any installation.

However, for easy usage, make the code executable using "chmod +x ASAP" and copy it to your system's PATH for invoking it from any folder.

## Usage ##
For the impatient:

Typical usage

`ASAP F/R/FR/FE/RE/FRE Forward_seq.ab1 Reverse_seq.ab1 Exons.fasta AminoAcidRef.fasta Reference_seq.fasta`

F: All sequences are forward.
R: All sequences are reverse.
FR: Sequences are in pairs, Forward and Reverse per sample.
Adding 'E' to the argument will enable exon extraction based on supplied exon file, its translation and alignment with reference amino acid fasta file.

For multiplex usage for n number of samples:

Example:
n = 2; That is, Sample A and B.

Typical Usage for such a scenario:

`ASAP F/R/FR/FE/RE/FRE Forward_seq_Sample_A.ab1 Forward_seq_Sample_B.ab1 Reverse_seq_Sample_A.ab1 Reverse_seq_Sample_B.ab1 Reference_Seq.fasta`


## Requirements ##
Linux/Unix System with EMBOSS, SEQTK, Python 3.x, BioPython & CLUSTALW2 installed.

Made and tested on Mac OS X 10.11.4 with Pyhton 3.4, BioPython 1.66, EMBOSS 6.5.7, SEQTK, CLUSTALW2 2.1

**Recently tested to successfully on a MacBook Pro with ARM based M1 Max CPU running macOS Ventura, 13.1 and Python 3.10.9**

You can download these freely available softwares from the below given links:

```text
EMBOSS : http://emboss.sourceforge.net/download/
SEQTK : https://github.com/lh3/seqtk
BioPython : http://biopython.org/wiki/Main_Page
NCBI BLAST+:  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  
```

## Description of results ##
```text
*._Forward_Sequence_QC.fastq : QC Trimmed Forward Sequence in FASTQ format.
*._Reverse_Sequence_QC.fastq : QC Trimmed Reverse complemented reverse sequence in FASTQ format.
*._consensus.fasta : Consensus sequence in FASTA format.
*._result.aln : Alignemnet of consensus with reference sequecne in ALN format.
*._result.dnd : Guide tree file for alignment in DND format.
*._result.fasta_report.txt : ClustalW2 alignment report.
*.consensus_BLAST.txt: Simplified BLAST result between Exonic sequence and extracted exon.
*.EXON_*: Files related with Exon extraction.
```

## Troubleshooting ##
A.	If ASAP is failing, you can provide `check` argument to it and it will tell you exactly what is failing. Install them and make sure they are in your system’s PATH. Also, make sure you followed the installation instructions carefully.

B.	Very poor alignment:

a.	Make sure you gave the sequences to ASAP in proper order. All forward sequences for all samples first, followed by all reverse sequences of all samples second and finally the reference sequence.

b.	Make sure you/ sequencing provider did not reverse complement the reverse sequences.

c.	Although, quality trimming is performed but still if the alignment is not at all right and all above mentioned troubleshooting has been followed, manually see the chromatogram files to check for background noise.

C.	I installed all the requirements but still ASAP is unable to use them!

a.	Make sure you have installed the requirements like BioPython with Python 3. If you have installed it specifically with a sub version of Python 3 like Python 3.5 and the program is not working, then kindly edit the first line of the code using the same method described in installation instruction as:
“#!/usr/bin/env python3” to “#!/usr/bin/env python3.5” or any other version that you installed your BioPython and other dependencies with.

Also, make sure you do not use Python 2.x. The program will only work with Python 3.x.

## Citation ##
If you use the ASAP in any of your work, then kindly cite the following:

Singh A, Bhatia P. Automated Sanger Analysis Pipeline (ASAP): A Tool for Rapidly Analyzing Sanger Sequencing Data with Minimum User Interference. J Biomol Tech [Internet]. 2016 Oct 17 [cited 2016 Nov 12]; Available from: http://www.ncbi.nlm.nih.gov/pubmed/27790076

Thank you!

*Please note that `ASAP` has recenlty gone through a major overhaul, it might be that a few of the descriptions written above will not match, but the functionality remains the same*

***Hope ASAP helps you in achieving things faster and better***
