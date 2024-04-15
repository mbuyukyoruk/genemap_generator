import argparse
import sys
import os
import subprocess
import textwrap
import shutil

try:
    from Bio import SeqIO
except:
    print("SeqIO module is not installed! Please install SeqIO and try again.")
    sys.exit()

try:
    import tqdm
except:
    print("tqdm module is not installed! Please install tqdm and try again.")
    sys.exit()

trimal_path = shutil.which("trimal")

if trimal_path is None:
    print("Trimal is required! Please install trimal and and make sure it is in $PATH.")
    sys.exit()

parser = argparse.ArgumentParser(prog='python genemap_generator.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=textwrap.dedent('''\

      	Author: Murat Buyukyoruk

        genemap_generator help:

This script is developed to generate a dataframe that can be used to create alignment panel with phylogenetic trees (i.e., to use with Pylo2genemap.py script)
tqdm is required to provide a progress bar.
trimal is required for -a option.

Syntax:

        genemap_generator.py -i demo_aligned.fasta -o demo_out.txt -d demo_hmmsearch_domtblout.txt -p genomic -a auto

genemap_generator dependencies:
	tqdm                                refer to https://pypi.org/project/tqdm/
	trimal                              refer to http://trimal.cgenomics.org/downloads

Input Paramaters (REQUIRED):
----------------------------
	-i/--input		FASTA			    Specify a fasta file used for generating phylogenetic tree and HMMER search. FASTA file requires headers starting with accession number. (i.e. >NZ_CP006019 [Prodigal Description])

	-o/--output		Output file	        Specify a output file name that should contain genemap info.

	-d/--data		HMMER domtblout		Specify a hmmsearch or hmmscan domtblout file to parse the positions of predicted hmm hits.

Parameters [optional]:
----------------------

	-p/--position	genomic or ORF		Option to generate ORF-specific or genomic-based genemap positions. It is recommended to use "genomic" option if you are planning to perform Neighbourhood analysis.

	-a/--anchor	    no or auto		    Option to generate ORF-specific or genomic-based genemap positions. It is recommended to use "genomic" option if you are planning to perform Neighbourhood analysis.

Basic Options:
--------------
	-h/--help		HELP			    Shows this help text and exits the run.

      	'''))
parser.add_argument('-i', '--input', required=True, type=str, dest='filename',
                    help='Specify a original fasta file.\n')
parser.add_argument('-d', '--data', nargs='+', required=True, type=str, dest='data',
                    help='Specify HMMSearch dombtlout data files. This flag accepts multiple file (i.e., '
                         '-d domtblout1.txt domtblout2.txt).\n')
parser.add_argument('-o', '--output', required=True, dest='out',
                    help='Specify a output file name.\n')
parser.add_argument('-p', '--position', required=False, dest='pos', default="ORF",
                    help='Specify type of position (i.e., genomic, ORF).\n')
parser.add_argument('-a', '--anchor', nargs='+', required=False, type=str, dest='anchor', default="no",
                    help='Option to set anchoring point for the genemap. Auto: This script will asign an auto '
                         'residue.\nNo: Will skip this step and Phylo_2_genemap script will use the HMMsearch '
                         'Domtblout to center align the target gene.\n<Accession> <residue> (i.e., '
                         '-a NZ_CP017479.1_2600 450) to provide a reference ORF and position of the residue.\n')

results = parser.parse_args()
filename = results.filename
data = results.data
out = results.out
pos = results.pos
anchor = results.anchor

orig_stdout = sys.stdout

os.system('> ' + out)

seq_id = []
seq_start = []
seq_end = []
seq_strand = []

f = open(out, 'a')
sys.stdout = f

print("molecule\tgenome\tgene\tstart\tend\torientation")

if anchor[0].lower() != "no":
    if anchor[0].lower() != "auto":
        ref_seq_id = anchor[0]
        ref_pos = len(anchor[1])
        record_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        ref_seq = str(record_dict[ref_seq_id].seq)
        for i in range(len(ref_seq)):
            if len(ref_seq[i].replace("-", "")) == ref_pos - 1:
                ali_pos = i
                break
    elif anchor[0].lower() == "auto":
        residue = []
        gap_prec = []
        timal_proc = subprocess.Popen("trimal -sgc -in " + filename, shell=True, stdout=subprocess.PIPE,text=True)
        trimal = timal_proc.communicate()[0].split('\n')
        for i in range(len(trimal)):
            if trimal[i] != "":
                if trimal[i][0] not in ["|", "+"]:
                    residue.append((trimal[i].split("\t")[0]).replace(" ", ""))
                    gap_prec.append(float(trimal[i].split("\t")[2]))
        min_gap = min(gap_prec)
        ind = gap_prec.index(min_gap)
        ali_pos = int(residue[ind])
    else:
        sys.stdout = orig_stdout
        print("Invalid option for anchor. Use 'no' or 'auto' or 'Accession position' options.")
        sys.exit()

proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length + 1)) as pbar:
    pbar.set_description('Getting sequence info...')
    for record in SeqIO.parse(filename, "fasta"):
        pbar.update()
        if pos.lower() == "orf":
            seq_len = len(str(record.seq).replace("-", ""))
            acc = record.id
            genome = acc.rsplit("_", 1)[0]
            print(acc + "\t" + acc + "\t" + genome + "\tlength\t1\t" + str(seq_len) + "\t1")
            if anchor[0].lower() != "no":
                ank = len(str(record.seq[:ali_pos]).replace("-", ""))
                print(acc + "\t" + acc + "\t" + genome + "\tAnchor\t" + str(ank) + "\t" + str(ank + 1) + "\t1")

        elif pos.lower() == "genomic":
            acc = record.id
            genome = acc.rsplit("_", 1)[0]
            start = int(record.description.split(" # ")[1])
            stop = int(record.description.split(" # ")[2])
            strand = str(record.description.split(" # ")[3]).replace("-1", "0")
            print(acc + "\t" + acc + "\t" + genome + "\tlength\t" + str(start) + "\t" + str(stop) + "\t" + strand)

            if anchor[0].lower() != "no":
                ank = len(str(record.seq[:ali_pos]).replace("-", ""))
                if strand == "1":
                    start_h = start + (ank * 3)
                    end_h = stop + ((ank + 1) * 3)
                elif strand == "0":
                    start_h = stop - (ank * 3)
                    end_h = stop - (ank + 1 * 3)
                print(
                    acc + "\t" + acc + "\t" + genome + "\tAnchor\t" + str(start_h) + "\t" + str(end_h) + "\t" + strand)

            seq_id.append(acc)
            seq_start.append(start)
            seq_end.append(stop)
            seq_strand.append(strand)
        else:
            sys.stdout = orig_stdout
            print("Invalid option for position argument. Use 'genomic' or 'ORF' options.")
            sys.exit()

for i in range(len(data)):
    proc = subprocess.Popen("wc -l < " + data[i], shell=True, stdout=subprocess.PIPE, text=True)
    length = int(proc.communicate()[0].split('\n')[0])

    proc = subprocess.Popen("grep 'Program' " + data[i], shell=True, stdout=subprocess.PIPE, text=True)
    proc_out = proc.communicate()[0].split('\n')[0]
    program = proc_out.split()[-1]

    with tqdm.tqdm(range(length + 1)) as pbar:
        pbar.set_description('Adding HMMsearch/HMMscan data... ' + data[i])
        with open(data[i], "r") as file:
            for line in file:
                pbar.update()
                if program == "hmmsearch":
                    if line[0] != "#":
                        if pos.lower() == "orf":
                            arr = line.split()
                            accession = arr[0]
                            genome_acc = accession.rsplit("_", 1)[0]
                            gene = arr[3]
                            start_h = arr[19]
                            end_h = arr[20]
                            print(accession + "\t" + accession + "\t" + genome_acc + "\t" + gene + "\t" + str(
                                start_h) + "\t" + str(end_h) + "\t1")
                        elif pos.lower() == "genomic":
                            arr = line.split()
                            accession = arr[0]
                            genome_acc = accession.rsplit("_", 1)[0]
                            gene = arr[3]
                            start_i = int(arr[19])
                            end_i = int(arr[20])
                            gene_start = int(line.split(" # ")[1])
                            gene_stop = int(line.split(" # ")[2])
                            gene_strand = line.split(" # ")[3].replace("-1", "0")
                            if gene_strand == "1":
                                start_h = gene_start + (start_i * 3)
                                end_h = gene_start + (end_i * 3)
                            elif gene_strand == "0":
                                start_h = gene_stop - (end_i * 3)
                                end_h = gene_stop - (start_i * 3)
                            print(accession + "\t" + accession + "\t" + genome_acc + "\t" + gene + "\t" + str(
                                start_h) + "\t" + str(end_h) + "\t" + gene_strand)
                elif program == "hmmscan":
                    if line[0] != "#":
                        if pos.lower() == "orf":
                            arr = line.split()
                            accession = arr[3]
                            genome_acc = accession.rsplit("_", 1)[0]
                            gene = arr[0]
                            start_h = arr[19]
                            end_h = arr[20]
                            print(accession + "\t" + accession + "\t" + genome_acc + "\t" + gene + "\t" + str(
                                start_h) + "\t" + str(end_h) + "\t1")
                        elif pos.lower() == "genomic":
                            arr = line.split()
                            accession = arr[3]
                            genome_acc = accession.rsplit("_", 1)[0]
                            gene = arr[0]
                            start_i = int(arr[19])
                            end_i = int(arr[20])
                            ind = seq_id.index(accession)
                            gene_start = int(seq_start[ind])
                            gene_stop = int(seq_end[ind])
                            gene_strand = seq_strand[ind].replace("-1", "0")
                            if gene_strand == "1":
                                start_h = gene_start + (start_i * 3)
                                end_h = gene_start + (end_i * 3)
                            elif gene_strand == "0":
                                start_h = gene_stop - (end_i * 3)
                                end_h = gene_stop - (start_i * 3)
                            print(accession + "\t" + accession + "\t" + genome_acc + "\t" + gene + "\t" + str(
                                start_h) + "\t" + str(end_h) + "\t" + gene_strand)
