#!/usr/bin/env python

import argparse
import array
import gzip
import logging
import os.path
import sys
import numpy
import matplotlib
from matplotlib.mlab import PCA
from numpy.testing import assert_array_almost_equal
import statsmodels.api as sm
from statsmodels.sandbox.tools import pca
from matplotlib import pyplot as plt

#from statsmodels.sandbox.tools.cross_val import LeaveOneOut



########################
########################
### Global constants ###
########################
########################

HOST_CODE  = 1
SYMB_CODE  = 2

BINARIES_DIR    = 'binaries'

LETTERS = ('A', 'T', 'G', 'C')

####################################################################
####################################################################
### Class to store the various paths used throughout the program ###
####################################################################
####################################################################

class pcaOptions:
    """This class consists of attributes to allow database and file paths to be
    obtained conveniently and consistently."""

    def __init__(self, args):
        self.args              = args
        self.suffix            = None


    def _getNumberOfSequences(self):
        """Return the length part of the kmer file name"""
        if self.args.numberOfSeq == 0:
            length = 'all'
        else:
            length = self.args.numberOfSeq
        return str(length)

    def _getSuffix(self):
        """Create the suffix of the SVM input files"""
        if not self.suffix:
            suffix = self._getNumberOfSequences()
            self.mink = str(self.args.minWordSize)
            self.maxk = str(self.args.maxWordSize)
            self.suffix = suffix + '_c' + self.mink + '_k' + self.maxk
        return self.suffix

######################
######################
### Misc utilities ###
######################
######################

def iterFasta(path):
    """Iterates over the sequences of a fasta file"""
    logging.info("Loading fasta files from %s" % path)
    name = None
    seq = []
    if path.endswith('.gz') or path.endswith('.gz"'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))
    handle.close()

############################
############################
### Compute Kmer vectors ###
############################
############################

def prepareMaps(k, maxk, kmers):
    """Prepares the kmer maps for the specified kmer range"""
    if k == maxk:
        n        = 0
        kmer2int = {}
        for kmer in kmers:
            kmer2int[kmer] = n
            n += 1
        return kmer2int
    newKmers = []
    for kmer in kmers:
        for letter in LETTERS:
            newKmers.append(kmer + letter)
    kmers = newKmers
    return prepareMaps(k + 1, maxk, kmers)


def computerKmers(args, path, label):
    """Compute the kmer counts throughout the kmer range for each sequence, and
    write the output to a file.  Each kmer counts will be scaled accordingly
    with the sequence size."""
    length = args.numberOfSeq
    # Prepare all maps
    kMin    = args.minWordSize
    kMax    = args.maxWordSize
    maps    = []
    logging.info('Preparing kmer maps')
    for i in xrange(kMin, kMax + 1):
        maps.append(prepareMaps(0, i, ['']))
    #obtain sequence size
    if length != 0:
        numberOfSeq = length
    else:
        l = 0
        for name,seq in iterFasta(path):
            l += 1
    	numberOfSeq = l	    
    # Initialise counts
    counts     = {}
    kmerArray = []
    labelList  = []
    kmerSize = 0
    for i in xrange(kMin, kMax + 1):
        counts[i] = array.array('d', [0 for x in xrange(4 ** i)])
    	kmerRange = 4 ** i
    	kmerSize  = kmerSize + kmerRange
    #Initialise data Matrix	
    kmerMatrix = numpy.zeros(shape=(numberOfSeq,kmerSize))
    # Iterate over sequences
    nSeqs   = 0
    logging.info('Computing kmers for %s' % path)
    for name, seq in iterFasta(path):
        if length > 0 and nSeqs >= length:
            break
        size   = len(seq)
        n      = 0
        labelList.append(label)
        #handle.write('%d' % label)
        # For each kmer value
        for i in xrange(kMin, kMax + 1):
            kCounts  = counts[i]
            # For word in the sequence
            for j in xrange(size - i + 1):
                word = seq[j:j + i]
                kMap = maps[i - kMin]
                idx  = kMap[word]
                kCounts[idx] += 1
            kCountsSum = sum(kCounts)
            for j in xrange(len(kCounts)):
                kCounts[j] /= kCountsSum
            kmerArray.append(kCounts)
        #Write to Matrix    
        t = 0
        for i in xrange(kMin, kMax + 1):
            for j in xrange(len(counts[i])):
	        kmerMatrix[nSeqs,t] = kmerArray[i - kMin][j] 
	        t += 1     
        #Reset counts
        for i in xrange(kMin, kMax + 1):
            for j in xrange(len(counts[i])):
	            counts[i][j] = 0
        nSeqs += 1
    # Trace
    logging.info('Processed %d sequences' % nSeqs)
    return kmerMatrix, labelList


def doPCA(dataMatrix):
    a,b,c,d= pca(dataMatrix, keepdim=2, normalize=True, demean=True)
    return a,b,c,d

	    
def plotPCA(args, kmerMatrix, labelList):
    xreduced, factors, evals, evecs = pca(kmerMatrix.T, keepdim=2, normalize=True, demean=True)
    pcaMatrix = evecs
    fig = matplotlib.pyplot.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel('Principal Component 1')
    ax.set_ylabel('Principal Component 2')
    ax.set_title('PCA Plot')
    for i in xrange(len(kmerMatrix)):
    	xpoint = pcaMatrix[i,0]
    	ypoint = pcaMatrix[i,1]
    	if labelList[i] == 1:
    		plt.scatter(xpoint, ypoint, color = 'green')
    	elif labelList[i] == 2:
    		plt.scatter(xpoint, ypoint, color = 'red')

    fig.savefig(args.outputFileName, format='png')
    matplotlib.pyplot.show()   

    
def checkExecutable(program):
    """Check whether a program is installed and executable"""

    # First check in $PATH
    path = os.getenv('PATH')
    for d in path.split(os.path.pathsep):
        exe = os.path.join(d, program)
        if os.path.exists(exe) and os.access(exe, os.X_OK):
            return exe
    # Then check in the subdirectory
    root = os.path.dirname(os.path.abspath(sys.argv[0]))
    exe  = os.path.join(root, BINARIES_DIR, program)
    if os.path.exists(exe) and os.access(exe, os.X_OK):
        return exe


def mainArgs():
    """Process command-line arguments"""

    parser = argparse.ArgumentParser(description='Perform SVM Classification of Host and Symbiont (or Parasite) Sequences')
    parser.add_argument('hostFile',
                        help='The 1st input fasta sequences')
    parser.add_argument('symbFile',
                        help='The 2nd input fasta sequences')
    #subparsers = parser.add_subparsers(help='Choose between option_1 or option_2 input format')
    #group = parser.add_mutually_exclusive_group()
    #parser_1 = subparsers.add_parser('option_1', help='Provide raw protein sequences, and perform blast before preparation for SVM')
    #parser_2 = subparsers.add_parser('option_2', help='Provide blast results as an input, directly start the preparation for SVM')
    parser.add_argument('-n',
                        '--numberOfSeq',
                        type=int,
                        default='0',
                        help='Maximum number of training & testing sequences')
    parser.add_argument('-c',
                        '--minWordSize',
                        type=int,
                        default='1',
                        help='Minimum value of DNA word length')
    parser.add_argument('-k',
                        '--maxWordSize',
                        type=int,
                        default='4',
                        help='Maxmimum value of DNA word length')
    parser.add_argument('-o',
                        '--outputFileName',
                        type=str,
                        default='output',
                        help='Name of the output file')
    args = parser.parse_args()
    if args.minWordSize > args.maxWordSize:
        logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
        sys.exit(1)
    return args

def main():
    """Banzai !!!"""
    logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
    args    = mainArgs()
    options = pcaOptions(args)
    logging.info("Arguments parsed. Starting...")
    hostPath = args.hostFile
    symbPath = args.symbFile
    hostKmer, hostList = computerKmers(args, hostPath, 1)
    symbKmer, symbList = computerKmers(args, symbPath, 2)
    
    kmerData = numpy.vstack([hostKmer, symbKmer])
    seqList  = numpy.hstack([hostList, symbList])
    plotPCA(args, kmerData, seqList)
    #blastClassification = loadBlastClassification(options)
    #predictSVM(args, blastClassification, kmerTrain, kmerTest)
    #print seqList
    logging.info("PCA analysis completed successfully.")

if __name__ == '__main__':
    main()
# vim:ts=4:sw=4:sts=4:et:ai:



