base_path = "queries/"
import pdb
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import math
import time
import os
from typing import Optional
from Bio import Entrez
handle = Entrez.esearch(db="nucleotide", retmax=300, term="virus[All Fields] AND ('unidentified plasmid'[Organism] OR plasmid[All Fields]) AND vector[All Fields]")
record = Entrez.read(handle)
gi_list = record["IdList"]
print(gi_list)
print("I found a gi list of length " + str(len(gi_list)))
gi_str = ",".join(gi_list)
handle = Entrez.efetch(db="nucleotide", id=gi_str, rettype="gb", retmode="text")
records = SeqIO.parse(handle, "gb")
window_size = 200
stride = 40
for record in records:
    sequence_length = len(record.seq)
    if (sequence_length < window_size) or (sequence_length > 1000*window_size):
        continue
    id = record.id
    os.mkdir(base_path + "analysis_sequences/windowed_sequences/" + str(id) + "/")
    os.mkdir("results/analysis_sequences/windowed_sequences/" + str(id) + "/")
    start = 0
    end = start + window_size
    i = 0
    while (end < sequence_length):
        window_fragment = record.seq[start:end]
        window_record = SeqRecord(window_fragment, str(id) + "_fragment_%i" % (i + 1), "", "")
        start += stride
        end += stride
        window_save_directory = "analysis_sequences/windowed_sequences/" + str(id)
        fragment_filename = "fragment_%i" % (i + 1)
        window_save_file_path = window_save_directory + "/fragment_%i.fsa" % (i + 1)
        SeqIO.write(window_record, base_path + window_save_file_path, "fasta")
        os.system("docker run --rm     -v $HOME/blastdb:/blast/blastdb:ro     -v $HOME/blastdb_custom:/blast/blastdb_custom:ro     -v $HOME/queries/%s:/blast/queries:ro     -v $HOME/results/%s:/blast/results:rw  ncbi/blast     blastn -query /blast/queries/%s.fsa -db blastdb/ref_viruses_rep_genomes     -out /blast/results/%s.out" % (window_save_directory, window_save_directory, fragment_filename, fragment_filename))
        i += 1
