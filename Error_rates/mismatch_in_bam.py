#!python
__doc__="""
    find all mismatches and class them into types (indel, different kinds of subsititution)
    Only one .fasta file in the prefix fold that contains the gene seq for the gene length
    make sure the gene name is exactly same with the one in the fasta file
    (no space or tab followed)
    """
__version__="1.0"
__author__="noah_yu_zhang"
__last_modify__="7_Aug_2019"

import argparse
import pysam
import glob
import os
from collections import OrderedDict
import numpy as np

from fasta_itertools import fasta_itertools

def make_pos_dict(path,gene):
    genefafile = glob.glob(path+"*.fa")[0]
    dseq = {}
    for fa in fasta_itertools(genefafile):
        dseq[fa[0]] = fa[1]

    #
    genelength = len(dseq[gene])
    a = range(0,genelength)
    dnt = OrderedDict()
    for i in a:
        dnt[i] = dict(zip(['A','T','C','G','ref','D','I'],[0]*7))

    for n,refnt in enumerate(list(dseq[gene])):
        #print (n,refnt,dnt[n])
        dnt[n]['ref'] = refnt.upper()
        #print (n,refnt,dnt[n])
    return dnt

def is_clip(align_pairs,n):
    downstream_clipper = False
    for tempreadpos,temprefpos,temprefnt in align_pairs[n:]:
        #print ("is_clip: ",tempreadpos,temprefpos,temprefnt)
        if temprefpos != None:
            downstream_clipper = False
            break
            #return False
        else:
            downstream_clipper = True
    if downstream_clipper == True:
        return True
    else:
        return False

def find_mismatches(bampath,prefix,outputfile,gene):
    dnt = make_pos_dict(prefix,gene)
    bam = pysam.AlignmentFile(bampath, "rb")
    readcount = 0
    for read in bam.fetch():
        #readcount += 1
        #print (readcount)
        #positions = read.get_reference_positions(full_length=False)
        #get_reference_sequence(self)
        align_pairs = read.get_aligned_pairs(with_seq=True)
        readseq = read.query_sequence
        #print (len(readseq), align_pairs[-1][0])
        previous_refpos = 0
        last_insert = 0
        for num,(readpos,refpos,refnt) in enumerate(align_pairs):

            #print ("input: ",num,readpos,refpos,refnt)
            if refpos == None: ## it may be a insertion or soft clip
                if previous_refpos == 0: ## soft clip at start terminal
                    pass
                else:
                    #print (is_clip(align_pairs,num))
                    if not is_clip(align_pairs,num):  ## to identify if it's a soft clip
                        if last_insert != previous_refpos+1:
                            #print ("before: ",(readpos,refpos,refnt),previous_refpos+1,dnt[previous_refpos],dnt[previous_refpos+1])
                            dnt[previous_refpos+1]['I'] += 1
                            #print ("double check insert: ",(readpos,refpos,refnt),previous_refpos+1,dnt[previous_refpos], dnt[previous_refpos+1])
                            last_insert = previous_refpos+1
                        else:
                            #print ("multiple_insert: ",previous_refpos+1,(readpos,refpos,refnt),dnt[previous_refpos], dnt[previous_refpos+1])
                            pass
                    else: ## soft clip at end terminal
                        #print ("the read end with soft clip")
                        break
            elif readpos == None: ## it is a deletion
                #print ("before: ",(readpos,refpos,refnt),previous_refpos,refpos,dnt[refpos])
                dnt[refpos]['D'] += 1
                #print ("double check deletion: ",(readpos,refpos,refnt),previous_refpos,refpos,dnt[refpos])
                previous_refpos = refpos
            elif refnt.islower:
                #print (readpos)
                readnt = readseq[readpos].upper()
                #print ("before: ",(readpos,refpos,readnt,refnt),previous_refpos,refpos,dnt[refpos])
                dnt[refpos][readnt] += 1
                #print ("double check sub: ",(readpos,refpos,readnt,refnt),previous_refpos,refpos,dnt[refpos])
                #print (readnt,dnt[refpos])
                previous_refpos = refpos
            #else:
                #previous_refpos = refpos
    bam.close()

    if os.path.exists(outputfile):
        os.remove(outputfile)
    with open(outputfile,'a') as output:
        output.write("\t".join(['pos','ref','A','T','C','G','D','I','sum','A%','T%','C%','G%','D%','I%','Match%'])+"\n")
        for k,v in dnt.items():
            counts = [v['A'],v['T'],v['C'],v['G'],v['D'],v['I']]
            total_count = sum(counts)
            rate = list(np.array(counts)/total_count)
            match_rate = v[v['ref']]/total_count
            ol = list(map(str, [k,v['ref']]+counts+[total_count]+rate+[match_rate]))
            #print (ol)
            output.write("\t".join(ol)+"\n")

if __name__=="__main__":

    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-b","--bam",help="input of nanopore bam file")
    parser.add_argument("-p","--prefix", help="path prefix for project")
    parser.add_argument("-o","--output",help="output file")
    parser.add_argument("-g","--gene",help="the gene to statistics")
    args=parser.parse_args()

    find_mismatches(args.bam,args.prefix,args.output,args.gene)
    #make_pos_dict(args.prefix)
