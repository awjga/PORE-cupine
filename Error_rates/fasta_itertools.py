#!python

__doc__="""
    A itertool for fasta
    """
__author__="noah_yu_zhang"
__version__="0.1 beta"
__last_modify__="27-05-2019"

import argparse
from sys import argv

def fasta_itertools(fapath):
    with open(fapath) as fafile:
        head = ""
        seq = ""
        for l in fafile:
            if l[0] == ">":
                if head == "":
                    pass
                else:
                    yield head,seq
                head = l.strip('\n')[1:]
                seq = ""
            else:
                seq = seq + l.strip('\n')
        # one time again for last one
        yield head,seq

if __name__=="__main__":
    for fa in fasta_itertools(argv[1]):
        print (fa)
