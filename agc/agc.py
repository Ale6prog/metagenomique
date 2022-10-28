#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import textwrap
import nwalign3 as nw
from tqdm import tqdm
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/


__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as  monfich:
        ephem = ""
        for line in monfich:
            if line.startswith(">"):
                if len(ephem) >= minseqlen:
                    yield ephem
                ephem = ""
                continue
            ephem += line.strip()
        yield ephem


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
 
    sequences = read_fasta(amplicon_file, minseqlen)
    for seq, nbr in sorted({seq: nbr for seq, nbr in Counter(sequences).items() if nbr >= mincount}.items(), 
                            key=lambda x: x[1] ,reverse=True):
        yield(seq, nbr)





def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    return sum([1 for i in range(len(alignment_list[0])) if alignment_list[0][i] == alignment_list[1][i]])/len(alignment_list[0]) * 100

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    ref = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    mem = []
    print("Comparing")
    [mem.append(seq) for seq in ref for seq_b in dereplication_fulllength(amplicon_file, minseqlen, mincount)
     if  get_identity(nw.global_align(seq[0], seq_b[0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath
     (os.path.join(os.path.dirname(__file__), "MATCH")))) <= 97]
    return mem

def write_OTU(otu_list, output_file):
    with open(output_file, "w") as file:
        for i, seq in enumerate(otu_list):
            file.write(f">OTU_{i+1} occurrence:{seq[1]}\n{textwrap.fill(seq[0], width=80)}\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,args.mincount, args.chunk_size, args.kmer_size)
    print("saving")
    write_OTU(otu, args.output_file)
    # Votre programme ici

#==============================================================
# Chimera removal section
#==============================================================
"""
def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def get_chunks(sequence, chunk_size):

    pass

def cut_kmer(sequence, kmer_size):

    pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass

def detect_chimera(perc_identity_matrix):
    pass

def search_mates(kmer_dict, sequence, kmer_size):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass
"""

if __name__ == '__main__':
    main()
