Author: Murat Buyukyoruk
       
    genemap_generator help:

This script is developed to generate a dataframe that can be used to create alignment panel with phylogenetic trees (i.e., to use with Pylo2genemap.py script)
tqdm is required to provide a progress bar.
trimal is required for -a option.

Syntax:

    genemap_generator.py -i demo_aligned.fasta -o demo_out.txt -d demo_hmmsearch_domtblout.txt -p genomic -a auto

genemap_generator dependencies:
	tqdm                            refer to https://pypi.org/project/tqdm/
	trimal                          refer to http://trimal.cgenomics.org/downloads

Input Paramaters (REQUIRED):
----------------------------
	-i/--input		FASTA			        Specify a fasta file used for generating phylogenetic tree and HMMER search. FASTA file requires headers starting with accession number. (i.e. >NZ_CP006019 [Prodigal Description])

	-o/--output		Output file	      Specify a output file name that should contain genemap info.

	-d/--data     HMMER domtblout		Specify a hmmsearch or hmmscan domtblout file to parse the positions of predicted hmm hits.

Parameters [optional]:
----------------------

	-p/--position  genomic or ORF		Option to generate ORF-specific or genomic-based genemap positions. It is recommended to use "genomic" option if you are planning to perform Neighbourhood analysis.

	-a/--anchor    no or auto		    Option to generate ORF-specific or genomic-based genemap positions. It is recommended to use "genomic" option if you are planning to perform Neighbourhood analysis.

Basic Options:
--------------
	-h/--help		HELP			          Shows this help text and exits the run.
