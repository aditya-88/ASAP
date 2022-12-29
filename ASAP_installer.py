#!/usr/bin/env python3
__author__ = 'Aditya Singh'
import os
import os.path
import easygui
import distutils.spawn
details = []
directory = (os.path.dirname(os.path.realpath(__file__)))
choice = easygui.boolbox(msg='Do you wish to install ASAP ?', title=' ', choices=('Continue', 'Exit'), image=directory+"/asap.gif")
if choice == 1:
    easygui.msgbox("Please choose your SeaView install directory in the next window.", ok_button="Choose SeaView Directory")
    seaview_directory = (easygui.diropenbox(msg="Choose the SeaView installation folder.", title="SeaView Directory"))+'/seaview'
    seaview = '"'+seaview_directory+'"'
else:
    exit()

easygui.msgbox("Please choose your SEQTK install directory in the next window", ok_button="Choose SEQTK Directory")
seqtk_directory = easygui.diropenbox(msg="Choose the SEQTK installation folder.", title="SEQTK Directory")+'/seqtk'
seqtk = '"'+seqtk_directory+'"'
check = 0
if os.path.isfile(seqtk_directory):
    print("\n1. SEQTK seems installed!")
    details.append("\n1. SEQTK seems installed!\n")
else:
    print("\n1. SEQTK was not found!\nPlease make sure it's installed properly.\n")
    check=+1
    details.append("\n1. SEQTK was not found!\nPlease make sure it's installed properly.\n")
if os.path.isfile(seaview_directory):
    print("\n2. SeaView seems installed!\n")
    details.append("2. SeaView seems installed!\n")
else:
    print("\n2. SeaView was not found!\nPlease make sure it's installed properly.\n")
    check=+1
    details.append("2. SeaView was not found!\nPlease make sure it's installed properly.\n")
try:
    from Bio import SeqIO
    from Bio.Blast import NCBIXML
    print("\n3. BioPython seems installed!\n")
    details.append("3. BioPython seems installed!\n")
except ImportError as e:
    print("\n3. Biopython is not installed! Please install the same and make sure the executables are available in your PATH\n")
    check=+1
    details.append("3. Biopython is not installed! Please install the same and make sure the executables are available in your PATH\n")
if distutils.spawn.find_executable('merger') is not None:
    print("\n4. EMBOSS seems installed.\n")
    details.append("4. EMBOSS seems installed.\n")
else:
    print("\n4. EMBOSS does not seems to be install. Install it and make sure it is available in System PATH\n")
    check=+1
    details.append("4. EMBOSS does not seems to be install. Install it and make sure it is available in System PATH\n")
if distutils.spawn.find_executable('clustalw2') is not None:
    print("\n5. CLUSTALW2 seems installed.\n")
    details.append("5. CLUSTALW2 seems installed.\n")
else:
    print("\n5. CLUSTALW2 does not seems to be install. Install it and make sure it is available in System PATH\n")
    check=+1
    details.append("5. CLUSTALW2 does not seems to be install. Install it and make sure it is available in System PATH\n")
if distutils.spawn.find_executable('blastn') is not None:
    print("\n6. NCBI BLAST+ seems to be intalled.\n")
    details.append("6. NCBI BLAST+ seems to be intalled.\n")
else:
    print("\nNCBI BLAST+ does not seems to be installed or not in PATH. Please check.\n")
    check=+1
    details.append("\nNCBI BLAST+ does not seems to be installed or not in PATH. Please check.\n")
if check > 0:
    print("\n{} out of 6 system tests failed. Please check the above mentioned errors and try again.\nThank you!".format(check))
    details.append("\n{} out of 6 system tests failed.\nPlease check the above mentioned errors and try again.\nThank you!".format(check))
    easygui.msgbox(''.join(details), ok_button="Exit")
    exit()
else:
    pass
with open(directory+"/ASAP.py", 'r') as infile:
    with open(directory+"/ASAP", 'w') as outfile:
        for line in infile:
            if line.startswith('seaview'):
                outfile.write("seaview = "+seaview+"\n")
                continue
            if line.startswith('seqtk'):
                outfile.write("seqtk = "+seqtk+"\n")
                continue
            if not line.startswith('seaview'):
                outfile.write(line)
                continue
            if not line.startswith('seqtk'):
                outfile.write(line)
                continue
os.system("chmod +x "+directory+"/ASAP")
easygui.msgbox("Installation done.\n\n"+"Will be using SeaView located at:\n"+seaview+"\n\nWill be using SEQTK located at:\n"+seqtk+"\n\nPlease find ASAP installed at:\n'{}/ASAP'".format(directory)+"\n\nYou might wish to move the ASAP executable in your PATH if you have super user access.\n\nThank you for installing ASAP!\nHope it helps you achieve better in lesser time.\n\nPlease cite ASAP if you use it in your work.")