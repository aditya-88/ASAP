#!/usr/bin/env python3
__author__ = 'aditya.onco@gmail.com'
import sys
import os.path
import distutils.spawn
import multiprocessing
import re
from pathlib import Path
threads = multiprocessing.cpu_count()
check = 0
#   Setting Up  #
version = "1.0"
print("Welcome to ASAP! version {}!\n".format(version))
citation = "Singh A, Bhatia P. Automated Sanger Analysis Pipeline (ASAP): A Tool for Rapidly Analyzing Sanger Sequencing Data with Minimum User Interference. J Biomol Tech [Internet]. 2016 Oct 17\nAvailable from: http://www.ncbi.nlm.nih.gov/pubmed/27790076"
#   EDIT THESE VALUES TO MATCH YOUR INSTALLATION    #
# Also, if execution fails, edit the first line to match with your Python version.
#SEQTK
seqtk = "/usr/local/bin/seqtk"
if "darwin" in sys.platform:
    clustal = "clustalw2"
else:
    clustal = "clustalw"

#   DO NOT EDIT ANYTHING BELOW THIS LINE    #

##################################################################################################################################
## Typical usage                                                                                                                ##
## ASAP! F/R/FR/FE/RE/FRE Forward_seq.ab1 Reverse_seq.ab1 Exons.fasta AminoAcidReference.fasta Reference_seq.fasta              ##
## Where:                                                                                                                       ##
## 'F' means you are giving only forward sequences for alignment                                                                ##
## 'R' means you are giving only reverse sequences for alignment                                                                ##
## 'FR' means you are giving both Forward and Reverse sequences for alignment                                                   ##
## Adding 'E' to the argument will enable exon extraction based on supplied exon file, its translation and alignment with       ##
## reference amino acid fasta file.                                                                                             ##
## ASAP! uses all the cores available from system for running BLAST.                                                            ##
## For multiplex usage for n number of samples:                                                                                 ##
## Example:                                                                                                                     ##
## n = 2; That is, Sample A and B.                                                                                              ##
## Typical Usage for such a sceneario:                                                                                          ##
## ASAP! Forward_seq_Sample_A.ab1 Forward_seq_Sample_B.ab1 Reverse_seq_Sample_A.ab1 Reverse_seq_Sample_B.ab1 Reference_Seq.fasta##
## Reverse Sequence should not be reverse complemented, the software does it on it's own.                                       ##
## You can manually open the *._result.aln file to visualize the alignment.                                                     ##
## Requirements:                                                                                                                ##
## Linux/Unix System with EMBOSS, SEQTK, Python 3.x, BioPython, NCBI BLAST+ Suite & CLUSTALW2 installed.                        ##
## Made and tested on Mac OS X 10.11.4 with Pyhton 3.4, BioPython 1.66, EMBOSS 6.5.7, SEQTK, CLUSTALW2 2.1                      ##
##                                                                                                                              ##
## Optional:                                                                                                                    ##
## Make the code executable using "chmod +x ASAP!" and copy it to your system's PATH for easy invoking it from any folder.      ##
##                                                                                                                              ##
## You can download these freely available softwares from the below given links:                                                ##
## EMBOSS       : http://emboss.sourceforge.net/download/                                                                       ##
## SEQTK        : https://github.com/lh3/seqtk                                                                                  ##
## BioPython    : http://biopython.org/wiki/Main_Page                                                                           ##
## NCBI BLAST+  : ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/                                                   ##
##                                                                                                                              ##
## Description of results:                                                                                                      ##
## *._Forward_Sequence_QC.fastq             : QC Trimmed Forward Sequence in FASTQ format.                                      ##
## *._Reverse_Sequence_QC.fastq             : QC Trimmed Reverse complemented reverse sequence in FASTQ format.                 ##
## *._consensus.fasta                       : Consensus sequence in FASTA format.                                               ##
## *._result.aln                            : Alignemnet of consensus with reference sequecne in ALN format.                    ##
## *._result.fasta_report.txt               : CLustalW2 alignment report.                                                       ##
## *.consensus_BLAST.txt                    : Simplified BLAST result between Exonic sequence and extracted exon.               ##
## *.EXON_*                                 : Files related with Exon extraction.                                               ##
## Failed_Result.txt                        : List of sequences failing quality tests.                                          ##
##                                                                                                                              ##
## Thanks for using the program!                                                                                                ##
## For any suggestions/queries please contact:                                                                                  ##
## Aditya Singh                                                                                                                 ##
## Ph.D. Scholar                                                                                                                ##
## Advanced Paediatrics Center                                                                                                  ##
## Post Graduate Institute of Medical Education and Research                                                                    ##
## Chandigarh, India                                                                                                            ##
## aditya.onco@gmail.com; adityas@live.in                                                                                       ##
##################################################################################################################################

# Check if all the required softwares are installed
def check_system():
    check = 0
    print("\nChecking system for installation of required programs and the commandline arguments.\nTotal 5 tests to perform.")
    if os.path.isfile(seqtk):
        print("\n1. SEQTK seems installed!\t\tPASS")
    else:
        print("\n1. SEQTK was not found!\nPlease make sure it's installed and location updated at the top of this file.\t\tFAIL")
        check+=1
    try:
        from Bio import SeqIO
        from Bio.Blast import NCBIXML
        print("2. BioPython seems installed!\t\tPASS")
    except ImportError as e:
        print("2. Biopython is not installed! Please install the same and make sure the executables are available in your PATH\t\tFAIL")
        check+=1
    if distutils.spawn.find_executable('merger') is not None:
        print("3. EMBOSS seems installed.\t\tPASS")
    else:
        print("3. EMBOSS does not seems to be install. Install it and make sure it is available in System PATH\t\tFAIL")
        check+=1
    if distutils.spawn.find_executable(clustal) is not None:
        print("4. CLUSTALW2 seems installed.\t\tPASS")
    else:
        print("4. CLUSTALW2 does not seems to be install. Install it and make sure it is available in System PATH\t\tFAIL")
        check+=1
    if distutils.spawn.find_executable('blastn') is not None:
        print("5. NCBI BLAST+ seems to be intalled.\tPASS")
    else:
        print("5. NCBI BLAST+ does not seems to be installed or not in PATH. Please check.")
        check+=1
    if check > 0:
        print("\n{} out of 5 system tests failed. Please check the above mentioned errors and try again.\nThank you!".format(check))
        exit()
    else:
        pass
    print("\nSystem check: PASS\nWelcome to Automated Sanger Analysis Pipeline (ASAP!) version {0} by Aditya Singh\naditya.onco@gmail.com.\nWill utilize all {1} available threads of the system for NCBI BLAST.".format(version, threads))
    exit()
    return

# Main function
def main():
    #   Check system and command    #
    if sys.argv[1].upper() != "CHECK" and len(sys.argv) < 3:
        print("Arguments missing!\n")
        print("Usage: ASAP! F/R/FR/FE/RE/FRE Forward_seq.ab1 Reverse_seq.ab1 Exons.fasta AminoAcidReference.fasta Reference_seq.fasta")
        print("Where:")
        print("'F' means you are giving only forward sequences for alignment")
        print("'R' means you are giving only reverse sequences for alignment")
        print("'FR' means you are giving both Forward and Reverse sequences for alignment")
        print("Adding 'E' to the argument will enable exon extraction based on supplied exon file, its translation and alignment with")
        print("reference amino acid fasta file.")
        print("ASAP! uses all the cores available from system for running BLAST.")
        print("For multiplex usage for n number of samples:")
        print("Example:")   
        print("n = 2; That is, Sample A and B.")
        print("Typical Usage for such a sceneario:")
        print("ASAP! Forward_seq_Sample_A.ab1 Forward_seq_Sample_B.ab1 Reverse_seq_Sample_A.ab1 Reverse_seq_Sample_B.ab1 Reference_Seq.fasta")
        print("Reverse Sequence should not be reverse complemented, the software does it on it's own.")
        print("You can manually open the *._result.aln file to visualize the alignment.")
        exit()
    else:
        print("Arguments found!\n")
    # Run check system if user wants to
    if sys.argv[1].upper() == "CHECK":
        check_system()
        exit()
    else:
        pass
    # Import rest of the modules
    from Bio import SeqIO
    from Bio.Blast import NCBIXML

    # Assign variables
    report_file_location = Path(sys.argv[2]).parent
    report_file = open(str(report_file_location)+"/Failed_Report.txt", "a")
    report_file.write("List of failed Sequencing runs:\nNumber\tSequence name\n")

    ## Importing modules
    import os.path
    from Bio.Align.Applications import ClustalwCommandline
    from Bio import AlignIO
    #_______________Multiplexing
    exons = ""
    all_sequences = sys.argv[2:-1]
    mode = sys.argv[1]
    number_of_sequences = int(len(all_sequences))
    if mode.upper() == "FR":
        print("Running in dual mode per sample (Forward and Reverse sequence per sample)\n")
        half = int(number_of_sequences/2)
        forward_sequences = all_sequences[:half]
        reverse_sequences = all_sequences[half:]
        count = 0
        reference_seq = sys.argv[-1]

    #   WORKING Dual mode#
        while (count <= half-1):
            try:
                Forward_seq = open(forward_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found.\nProgram will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(forward_sequences[count]))
                exit()
            try:
                Reverse_seq = open(reverse_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found. Program will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(reverse_sequences[count]))
                exit()
            forward_file = open(Forward_seq.name[:-4]+"_Forward_Sequence.fastq", 'w+b')
            reverse_file_not_rev = open(Forward_seq.name[:-4]+"_Reverse_Sequence_Not_Reversed.fastq", 'w+b')
            reverse_file = open(Forward_seq.name[:-4]+"_Reverse_Sequence.fastq", 'w+b')
            Forward_fastq = SeqIO.convert(Forward_seq, "abi",forward_file.name, "fastq" )
            Reverse_fastq = SeqIO.convert(Reverse_seq, "abi",reverse_file_not_rev.name, "fastq" )
            os.system(seqtk+" seq -r "+reverse_file_not_rev.name+" > "+reverse_file.name)
            print("Performing Quality control and building consensus sequence....")
            os.system(seqtk+" trimfq -q 0.05 "+forward_file.name+" > "+forward_file.name[:-6] + "_QC.fastq")
            os.system(seqtk+" trimfq -q 0.05 "+reverse_file.name+" > "+reverse_file.name[:-6] + "_QC.fastq")
            os.system(seqtk + " seq -q23 -n N " + reverse_file.name[:-6] + "_QC.fastq" + " >" + reverse_file.name[:-23] + "_N_.fasta")
            os.system(seqtk + " seq -q23 -n N " + forward_file.name[:-6] + "_QC.fastq" + " >" + forward_file.name[:-23] + "_N_.fasta")
            os.system("merger "+"-verbose "+ "-asequence "+forward_file.name[:-23] + "_N_.fasta"+" -bsequence "+reverse_file.name[:-23] + "_N_.fasta"+" -outfile "+reverse_file.name[:-23]+"_consensus_alignment.aln "+"-outseq "+reverse_file.name[:-23]+"_consensus.fasta")
            count = count+1
            size = os.stat(reverse_file.name[:-23]+"_consensus.fasta")
            if size.st_size > 100:
                print("Performing Alignment with reference sequence....")
                os.system("cat "+reference_seq+" "+reverse_file.name[:-23]+"_consensus.fasta "+"> "+reverse_file.name[:-23]+"_result.fasta")
                in_file = reverse_file.name[:-23]+"_result.fasta"
                print("Chosen file:\n",in_file)
                clustalw_cline = ClustalwCommandline(clustal, infile=in_file)
                print("Running program....\n",clustalw_cline)
                os.system(str(clustalw_cline)+">"+in_file+"_report.txt")
                alignment = AlignIO.read(in_file[:-6]+".aln", "clustal")
            else:
                print("\n\nA very short consensus sequence of {} bytes made. Not using it. Check sequence quality.\n\n".format(size.st_size))
                report_file.write(str(count)+"\t"+reverse_file.name[:-23]+".ab1\n")
    #   Renaming and Purging of unimportant files    #
            os.system("rm "+reverse_file.name[:-23]+"_consensus_alignment.aln "+Forward_seq.name[:-4]+"_Reverse_Sequence_Not_Reversed.fastq "+forward_file.name+" "+reverse_file.name+" ")

    # EXON Extraction___________________________________________

    elif mode.upper() == "FRE":
        print("Running in dual mode per sample (Forward and Reverse sequence per sample) and also extracting Exonic region\n")
        all_sequences = sys.argv[2:-3]
        exons = sys.argv[-3]
        amino_acid = sys.argv[-2]
        number_of_sequences = int(len(all_sequences))
        half = int(number_of_sequences/2)
        forward_sequences = all_sequences[:half]
        reverse_sequences = all_sequences[half:]
        count = 0
        reference_seq = sys.argv[-1]
        os.system("makeblastdb -in "+exons+" -out "+exons+"_DB"+" -dbtype nucl")
    #   WORKING Dual mode#
        while (count <= half-1):
            try:
                Forward_seq = open(forward_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found.\nProgram will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(forward_sequences[count]))
                exit()
            try:
                Reverse_seq = open(reverse_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found. Program will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(reverse_sequences[count]))
                exit()
            forward_file = open(Forward_seq.name[:-4]+"_Forward_Sequence.fastq", 'w+b')
            reverse_file_not_rev = open(Forward_seq.name[:-4]+"_Reverse_Sequence_Not_Reversed.fastq", 'w+b')
            reverse_file = open(Forward_seq.name[:-4]+"_Reverse_Sequence.fastq", 'w+b')
            Forward_fastq = SeqIO.convert(Forward_seq, "abi",forward_file.name, "fastq" )
            Reverse_fastq = SeqIO.convert(Reverse_seq, "abi",reverse_file_not_rev.name, "fastq" )
            os.system(seqtk+" seq -r "+reverse_file_not_rev.name+" > "+reverse_file.name)
            print("Performing Quality control and building consensus sequence....")
            os.system(seqtk+" trimfq -q 0.05 "+forward_file.name+" > "+forward_file.name[:-6] + "_QC.fastq")
            os.system(seqtk+" trimfq -q 0.05 "+reverse_file.name+" > "+reverse_file.name[:-6] + "_QC.fastq")
            os.system(seqtk + " seq -q23 -n N " + reverse_file.name[:-6] + "_QC.fastq" + " >" + reverse_file.name[:-23] + "_N_.fasta")
            os.system(seqtk + " seq -q23 -n N " + forward_file.name[:-6] + "_QC.fastq" + " >" + forward_file.name[:-23] + "_N_.fasta")
            os.system("merger "+"-verbose "+ "-asequence "+forward_file.name[:-23] + "_N_.fasta"+" -bsequence "+reverse_file.name[:-23] + "_N_.fasta"+" -outfile "+reverse_file.name[:-23]+"_consensus_alignment.aln "+"-outseq "+reverse_file.name[:-23]+"_consensus.fasta")
            count = count+1
            size = os.stat(reverse_file.name[:-23]+"_consensus.fasta")
            if size.st_size > 100:
                #print(os.stat(reverse_file.name[:-23]+"_consensus.fasta".st_size))
                os.system("blastn -query "+reverse_file.name[:-23]+"_consensus.fasta -db "+exons+"_DB -out "+reverse_file.name[:-23]+"_consensus_BLAST.xml -max_hsps 1 -num_alignments 1 -outfmt 5 -num_threads "+str(threads))
                result = open(reverse_file.name[:-23]+"_consensus_BLAST.xml", "r")
                outfile = open(reverse_file.name[:-23]+"_consensus_BLAST.txt", "w")
                blast_records = NCBIXML.parse(result)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            extracted_exon = hsp.query
                            extracted_exon = re.sub('-', "", extracted_exon)
                            os.system("echo '>"+os.path.basename(reverse_file.name)[:-23]+ "\n"+extracted_exon+"' > "+reverse_file.name[:-23]+"_EXON.fasta")
                            pi = float(hsp.identities)/float(alignment.length)*100
                            pi_deci = round(pi,ndigits=2)
                            print("Percent Identiy with the Exonic region",pi_deci)
                            print ('sequence:', alignment.title,"\t", len(hsp.sbjct))
                            outfile.write(os.path.basename(reverse_file.name)[:-23]+"\n")
                            outfile.write('\n****Alignment****\n\n')
                            outfile.write("Sequence: "+ str(alignment.title) + "\n")
                            outfile.write("Percent Identity: "+ str(pi_deci)+"%"+ "\n")
                            outfile.write("Length: "+ str(alignment.length) + "\n")
                            outfile.write("e-value: "+ str(hsp.expect) + "\n")
                            outfile.write("Gaps: "+ str(hsp.gaps) + "\n")
                            outfile.write("Query:   "+ str(hsp.query[0:75]) + "\n")
                            outfile.write("Match:   "+ str(hsp.match[0:75]) + "\n")
                            outfile.write("Subject: "+ str(hsp.sbjct[0:75]) + "\n")
                os.system("degapseq -sequence "+reverse_file.name[:-23]+"_EXON.fasta"+" -outseq "+reverse_file.name[:-23])
                print("Translating input sequence at frames +1, +2 & +3...\n")
                os.system("transeq -sequence "+reverse_file.name[:-23]+"_EXON.fasta"+" -outseq "+reverse_file.name[:-23]+"_EXON_AminoAcid.fasta"+" -frame 6")
                os.system("cat "+amino_acid+" "+reverse_file.name[:-23]+"_EXON_AminoAcid.fasta"+" "+" > "+reverse_file.name[:-23]+"_EXON_AminoAcid_with_Reference.fasta")
                in_file_amino_acid = reverse_file.name[:-23]+"_EXON_AminoAcid_with_Reference.fasta"
                clustalw_cline_amino_acid = ClustalwCommandline(clustal, infile=in_file_amino_acid)
                print("Aligning translated amino acid with reference protein sequence '{}'\n".format(amino_acid))
                os.system(str(clustalw_cline_amino_acid)+">"+in_file_amino_acid+"_report.txt")
                print("Performing nucleic acid alignment with reference sequence....")
                os.system("cat "+reference_seq+" "+reverse_file.name[:-23]+"_consensus.fasta "+"> "+reverse_file.name[:-23]+"_result.fasta")
                in_file = reverse_file.name[:-23]+"_result.fasta"
                print("Chosen file:\n",in_file)
                clustalw_cline = ClustalwCommandline(clustal, infile=in_file)
                print("Running program....\n",clustalw_cline)
                os.system(str(clustalw_cline)+">"+in_file+"_report.txt")
                alignment = AlignIO.read(in_file[:-6]+".aln", "clustal")
            else:
                print("\n\nA very short consensus sequence of {} bytes made. Not using it. Check sequence quality.\n\n".format(size.st_size))
                report_file.write(str(count)+"\t"+reverse_file.name[:-23]+".ab1\n")
            #   Renaming and Purging of unimportant files    #
                os.system("rm "+reverse_file.name[:-23]+"_consensus_alignment.aln "+Forward_seq.name[:-4]+"_Reverse_Sequence_Not_Reversed.fastq "+forward_file.name+" "+reverse_file.name+" ")
    elif mode.upper() == "F":
        print("Running in Forward ONLY mode per sample (Only Forward sequence per sample)\n")
        forward_sequences = all_sequences
        count = 0
        reference_seq = sys.argv[-1]
    #   WORKING forward only mode#
        while (count <= number_of_sequences-1):
            try:
                Forward_seq = open(forward_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found.\nProgram will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(forward_sequences[count]))
                exit()
            forward_file = open(Forward_seq.name[:-4]+"_Forward_Sequence.fastq", 'w+b')
            Forward_fastq = SeqIO.convert(Forward_seq, "abi",forward_file.name, "fastq" )
            print("Performing Quality control....")
            os.system(seqtk+" trimfq -q 0.05 "+forward_file.name+" > "+forward_file.name[:-6] + "_QC.fastq")
            os.system(seqtk+" seq -q23 -n N "+forward_file.name[:-6] + "_QC.fastq"+" >"+forward_file.name[:-23]+"_consensus.fasta")
            #os.system(seqtk +" seq -a "+forward_file.name[:-6] + "_QC.fastq >"+forward_file.name[:-23]+"_consensus.fasta")
            count = count+1
            size = os.stat(forward_file.name[:-23]+"_consensus.fasta")
            if size.st_size > 100:
                print("Performing Alignment with reference sequence....")
                os.system("cat "+reference_seq+" "+forward_file.name[:-23]+"_consensus.fasta "+"> "+forward_file.name[:-23]+"_result.fasta")
                in_file = forward_file.name[:-23]+"_result.fasta"
                print("Chosen file:\n",in_file)
                clustalw_cline = ClustalwCommandline(clustal, infile=in_file)
                print("Running program....\n",clustalw_cline)
                os.system(str(clustalw_cline)+">"+in_file+"_report.txt")
                alignment = AlignIO.read(in_file[:-6]+".aln", "clustal")
            else:
                print("\n\nA very short consensus sequence of {} bytes made. Not using it. Check sequence quality.\n\n".format(size.st_size))
                report_file.write(str(count)+"\t"+forward_file.name[:-23]+".ab1"+".ab1\n")
    # EXON Extraction___________________________________________
    elif mode.upper() == "FE":
        all_sequences = sys.argv[2:-3]
        exons = sys.argv[-3]
        amino_acid = sys.argv[-2]
        print("Running in Forward ONLY mode per sample (Only Forward sequence per sample) with Exon extraction\n")
        forward_sequences = all_sequences
        count = 0
        reference_seq = sys.argv[-1]
        os.system("makeblastdb -in "+exons+" -out "+exons+"_DB"+" -dbtype nucl")
    #   WORKING forward only mode#
        while (count <= number_of_sequences-3):
            try:
                Forward_seq = open(forward_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found.\nProgram will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(forward_sequences[count]))
                exit()
            forward_file = open(Forward_seq.name[:-4]+"_Forward_Sequence.fastq", 'w+b')
            Forward_fastq = SeqIO.convert(Forward_seq, "abi",forward_file.name, "fastq" )
            print("Performing Quality control....")
            os.system(seqtk+" trimfq -q 0.05 "+forward_file.name+" > "+forward_file.name[:-6] + "_QC.fastq")
            os.system(seqtk + " seq -q23 -n N " + forward_file.name[:-6] + "_QC.fastq" + " >" + forward_file.name[:-23] + "_consensus.fasta")
            #os.system(seqtk +" seq -a "+forward_file.name[:-6] + "_QC.fastq >"+forward_file.name[:-23]+"_consensus.fasta")
            count = count+1
            size = os.stat(forward_file.name[:-23]+"_consensus.fasta")
            if size.st_size > 100:
                os.system("blastn -query "+forward_file.name[:-23]+"_consensus.fasta -db "+exons+"_DB -out "+forward_file.name[:-23]+"_consensus_BLAST.xml -max_hsps 1 -num_alignments 1 -outfmt 5 -num_threads "+str(threads))
                result = open(forward_file.name[:-23]+"_consensus_BLAST.xml", "r")
                outfile = open(forward_file.name[:-23]+"_consensus_BLAST.txt", "w")
                blast_records = NCBIXML.parse(result)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            extracted_exon = hsp.query
                            extracted_exon = re.sub('-', "", extracted_exon)
                            os.system("echo '>"+os.path.basename(forward_file.name)[:-23]+ "\n"+extracted_exon+"' > "+forward_file.name[:-23]+"_EXON.fasta")
                            pi = float(hsp.identities)/float(alignment.length)*100
                            pi_deci = round(pi,ndigits=2)
                            print("Percent Identiy with the Exonic region",pi_deci)
                            print ('sequence:', alignment.title,"\t", len(hsp.sbjct))
                            outfile.write(os.path.basename(forward_file.name)[:-23]+"\n")
                            outfile.write('\n****Alignment****\n\n')
                            outfile.write("Sequence: "+ str(alignment.title) + "\n")
                            outfile.write("Percent Identity: "+ str(pi_deci)+"%"+ "\n")
                            outfile.write("Length: "+ str(alignment.length) + "\n")
                            outfile.write("e-value: "+ str(hsp.expect) + "\n")
                            outfile.write("Gaps: "+ str(hsp.gaps) + "\n")
                            outfile.write("Query:   "+ str(hsp.query[0:75]) + "\n")
                            outfile.write("Match:   "+ str(hsp.match[0:75]) + "\n")
                            outfile.write("Subject: "+ str(hsp.sbjct[0:75]) + "\n")
                os.system("degapseq -sequence "+forward_file.name[:-23]+"_EXON.fasta"+" -outseq "+forward_file.name[:-23])
                os.system("transeq -sequence "+forward_file.name[:-23]+"_EXON.fasta"+" -outseq "+forward_file.name[:-23]+"_EXON_AminoAcid.fasta"+" -frame 6")
                print("Translating input sequence at frames +1, +2 & +3...\n")
                os.system("transeq -sequence "+forward_file.name[:-23]+"_EXON.fasta"+" -outseq "+forward_file.name[:-23]+"_EXON_AminoAcid.fasta"+" -frame 6")
                os.system("cat "+amino_acid+" "+forward_file.name[:-23]+"_EXON_AminoAcid.fasta"+" "+" > "+forward_file.name[:-23]+"_EXON_AminoAcid_with_Reference.fasta")
                in_file_amino_acid = forward_file.name[:-23]+"_EXON_AminoAcid_with_Reference.fasta"
                clustalw_cline_amino_acid = ClustalwCommandline(clustal, infile=in_file_amino_acid)
                print("Aligning translated amino acid with reference protein sequence '{}'\n".format(amino_acid))
                os.system(str(clustalw_cline_amino_acid)+">"+in_file_amino_acid+"_report.txt")
                print("Performing Alignment with reference sequence....")
                os.system("cat "+reference_seq+" "+forward_file.name[:-23]+"_consensus.fasta "+"> "+forward_file.name[:-23]+"_result.fasta")
                in_file = forward_file.name[:-23]+"_result.fasta"
                print("Chosen file:\n",in_file)
                clustalw_cline = ClustalwCommandline(clustal, infile=in_file)
                print("Running program....\n",clustalw_cline)
                os.system(str(clustalw_cline)+">"+in_file+"_report.txt")
                alignment = AlignIO.read(in_file[:-6]+".aln", "clustal")
            else:
                print("\n\nA very short consensus sequence of {} bytes made. Not using it. Check sequence quality.\n\n".format(size.st_size))
                report_file.write(str(count)+"\t"+forward_file.name[:-23]+".ab1\n")
    elif mode.upper() == "R":
        print("Running in Reverse ONLY mode per sample (Only Reverse sequence per sample)\n")
        reverse_sequences = all_sequences
        count = 0
        reference_seq = sys.argv[-1]
    #   WORKING reverse only mode#
        while (count <= number_of_sequences-1):
            try:
                Reverse_seq = open(reverse_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found. Program will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(reverse_sequences[count]))
                exit()
            reverse_file_not_rev = open(Reverse_seq.name[:-4]+"_Reverse_Sequence_Not_Reversed.fastq", 'w+b')
            reverse_file = open(Reverse_seq.name[:-4]+"_Reverse_Sequence.fastq", 'w+b')
            Reverse_fastq = SeqIO.convert(Reverse_seq, "abi",reverse_file_not_rev.name, "fastq" )
            os.system(seqtk+" seq -r "+reverse_file_not_rev.name+" > "+reverse_file.name)
            print("Performing Quality control and building consensus sequence....")
            os.system(seqtk+" trimfq -q 0.05 "+reverse_file.name+" > "+reverse_file.name[:-6] + "_QC.fastq")
            os.system(seqtk + " seq -q23 -n N " + reverse_file.name[:-6] + "_QC.fastq" + " >" + reverse_file.name[:-23] + "_consensus.fasta")
            #os.system(seqtk+" seq -a "+reverse_file.name[:-6] + "_QC.fastq >"+reverse_file.name[:-23]+"_consensus.fasta")
            count = count+1
            size = os.stat(reverse_file.name[:-23]+"_consensus.fasta")
            if size.st_size > 100:
                print("Performing Alignment with reference sequence....")
                os.system("cat "+reference_seq+" "+reverse_file.name[:-23]+"_consensus.fasta "+"> "+reverse_file.name[:-23]+"_result.fasta")
                in_file = reverse_file.name[:-23]+"_result.fasta"
                print("Chosen file:\n",in_file)
                clustalw_cline = ClustalwCommandline(clustal, infile=in_file)
                print("Running program....\n",clustalw_cline)
                os.system(str(clustalw_cline)+">"+in_file+"_report.txt")
                alignment = AlignIO.read(in_file[:-6]+".aln", "clustal")
            else:
                print("\n\nA very short consensus sequence of {} bytes made. Not using it. Check sequence quality.\n\n".format(size.st_size))
                report_file.write(str(count)+"\t"+reverse_file.name[:-23]+".ab1\n")

    # EXON Extraction___________________________________________

    elif mode.upper() == "RE":
        print("Running in Reverse ONLY mode per sample (Only Reverse sequence per sample) with Exon extraction\n")
        all_sequences = sys.argv[2:-3]
        exons = sys.argv[-3]
        amino_acid = sys.argv[-2]
        reverse_sequences = all_sequences
        count = 0
        reference_seq = sys.argv[-1]
        os.system("makeblastdb -in "+exons+" -out "+exons+"_DB"+" -dbtype nucl")
    #   WORKING reverse only mode#
        while (count <= number_of_sequences-3):
            try:
                Reverse_seq = open(reverse_sequences[count], "rb")
            except FileNotFoundError:
                print("\nCheck the file address.\n'{}' : File not found. Program will exit now.\nResults might be incomplete.\nPlease verify the file addresses and try again.".format(reverse_sequences[count]))
                exit()
            reverse_file_not_rev = open(Reverse_seq.name[:-4]+"_Reverse_Sequence_Not_Reversed.fastq", 'w+b')
            reverse_file = open(Reverse_seq.name[:-4]+"_Reverse_Sequence.fastq", 'w+b')
            Reverse_fastq = SeqIO.convert(Reverse_seq, "abi",reverse_file_not_rev.name, "fastq" )
            os.system(seqtk+" seq -r "+reverse_file_not_rev.name+" > "+reverse_file.name)
            print("Performing Quality control and building consensus sequence....")
            os.system(seqtk+" trimfq -q 0.05 "+reverse_file.name+" > "+reverse_file.name[:-6] + "_QC.fastq")
            os.system(seqtk + " seq -q23 -n N " + reverse_file.name[:-6] + "_QC.fastq" + " >" + reverse_file.name[:-23] + "_consensus.fasta")
            #os.system(seqtk+" seq -a "+reverse_file.name[:-6] + "_QC.fastq >"+reverse_file.name[:-23]+"_consensus.fasta")
            count = count+1
            size = os.stat(reverse_file.name[:-23]+"_consensus.fasta")
            if size.st_size > 100:
                os.system("blastn -query "+reverse_file.name[:-23]+"_consensus.fasta -db "+exons+"_DB -out "+reverse_file.name[:-23]+"_consensus_BLAST.xml -max_hsps 1 -num_alignments 1 -outfmt 5 -num_threads "+str(threads))
                result = open(reverse_file.name[:-23]+"_consensus_BLAST.xml", "r")
                outfile = open(reverse_file.name[:-23]+"_consensus_BLAST.txt", "w")
                blast_records = NCBIXML.parse(result)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            extracted_exon = hsp.query
                            extracted_exon = re.sub('-', "", extracted_exon)
                            os.system("echo '>"+os.path.basename(reverse_file.name)[:-23]+ "\n"+extracted_exon+"' > "+reverse_file.name[:-23]+"_EXON.fasta")
                            pi = float(hsp.identities)/float(alignment.length)*100
                            pi_deci = round(pi,ndigits=2)
                            print("Percent Identiy with the Exonic region",pi_deci)
                            print ('sequence:', alignment.title,"\t", len(hsp.sbjct))
                            outfile.write(os.path.basename(reverse_file.name)[:-23]+"\n")
                            outfile.write('\n****Alignment****\n\n')
                            outfile.write("Sequence: "+ str(alignment.title) + "\n")
                            outfile.write("Percent Identity: "+ str(pi_deci)+"%"+ "\n")
                            outfile.write("Length: "+ str(alignment.length) + "\n")
                            outfile.write("e-value: "+ str(hsp.expect) + "\n")
                            outfile.write("Gaps: "+ str(hsp.gaps) + "\n")
                            outfile.write("Query:   "+ str(hsp.query[0:75]) + "\n")
                            outfile.write("Match:   "+ str(hsp.match[0:75]) + "\n")
                            outfile.write("Subject: "+ str(hsp.sbjct[0:75]) + "\n")
                os.system("degapseq -sequence "+reverse_file.name[:-23]+"_EXON.fasta"+" -outseq "+reverse_file.name[:-23])
                os.system("transeq -sequence "+reverse_file.name[:-23]+"_EXON.fasta"+" -outseq "+reverse_file.name[:-23]+"_EXON_AminoAcid.fasta"+" -frame 6")
                print("Translating input sequence at frames +1, +2 & +3...\n")
                os.system("transeq -sequence "+reverse_file.name[:-23]+"_EXON.fasta"+" -outseq "+reverse_file.name[:-23]+"_EXON_AminoAcid.fasta"+" -frame 6")
                os.system("cat "+amino_acid+" "+reverse_file.name[:-23]+"_EXON_AminoAcid.fasta"+" "+" > "+reverse_file.name[:-23]+"_EXON_AminoAcid_with_Reference.fasta")
                in_file_amino_acid = reverse_file.name[:-23]+"_EXON_AminoAcid_with_Reference.fasta"
                clustalw_cline_amino_acid = ClustalwCommandline(clustal, infile=in_file_amino_acid)
                print("Aligning translated amino acid with reference protein sequence '{}'\n".format(amino_acid))
                os.system(str(clustalw_cline_amino_acid)+">"+in_file_amino_acid+"_report.txt")
                print("\nPerforming Alignment with reference sequence....")
                os.system("cat "+reference_seq+" "+reverse_file.name[:-23]+"_consensus.fasta "+"> "+reverse_file.name[:-23]+"_result.fasta")
                in_file = reverse_file.name[:-23]+"_result.fasta"
                print("Chosen file:\n",in_file)
                clustalw_cline = ClustalwCommandline(clustal, infile=in_file)
                print("Running program....\n",clustalw_cline)
                os.system(str(clustalw_cline)+">"+in_file+"_report.txt")
                alignment = AlignIO.read(in_file[:-6]+".aln", "clustal")
            else:
                print("\n\nA very short consensus sequence of {} bytes made. Not using it. Check sequence quality.\n\n".format(size.st_size))
                report_file.write(str(count)+"\t"+reverse_file.name[:-23]+".ab1\n")
    else:
        print("'{}' is not a valid mode argument.\nMake sure you gave the mode argument correct (F/FE/FR/FRE/R/RE)\nPlease try Again".format(mode))
        exit()
    #   End message #
    print("Thank you for using ASAP! version {0}!\nHope it was useful to you!\n\nIf you use ASAP! in your work, please cite:\n\n{1}".format(version,citation))
    return

# Run ASAP! #
if __name__ == "__main__":
    main()