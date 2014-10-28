#!/usr/bin/env python

#################################################################################
#################################################################################
############PSYTRANS V3.0 Using Reverse Complementary Strand for SVM Training####
#################################################################################
#################################################################################

import argparse
import array
import gzip
import logging
import os.path
import subprocess
import sys
import shutil
import os
import threading
import traceback
import math
import random

random.seed(1234)
if(sys.hexversion < 0x03000000):
    import Queue
else:
    import queue as Queue


########################
########################
### Global constants ###
########################
########################

HOST_NAME  = 'host'
SYMB_NAME  = 'symb'
DB_NAME    = 'HostSymbDB'
DB_FASTA   = DB_NAME + '.fasta'
BLAST_FILE = HOST_NAME + SYMB_NAME + '_blastResults.txt'
BLAST_SORT = HOST_NAME + SYMB_NAME + '_blastClassification.txt'
HOST_CODE  = 1
SYMB_CODE  = 2

#BLAST CLASSIFICATION VARIABLES
MIN_BIT_RATIO = 2
MIN_BIT_DELTA = 100

# SVM GLOBAL VARIABLES
SVM_FOLD   = 5
SVM_CSTART = -5
SVM_CEND   = 15
SVM_CSTEP  = 2
SVM_GSTART = 3
SVM_GEND   = -15
SVM_GSTEP  = -2

HOST_TRAINING = HOST_NAME + '_training.fasta'
HOST_TESTING  = HOST_NAME + '_testing.fasta'
SYMB_TRAINING = SYMB_NAME + '_training.fasta'
SYMB_TESTING  = SYMB_NAME + '_testing.fasta'
BINARIES_DIR    = 'binaries'

LETTERS = ('A', 'T', 'G', 'C')

####################################################################
####################################################################
### Class to store the various paths used throughout the program ###
####################################################################
####################################################################

class PsyTransOptions:
    """This class consists of attributes to allow database and file paths to be
    obtained conveniently and consistently."""

    def __init__(self, args):
        self.args              = args
        self.dbPath            = None
        self.fastaDbPath       = None
        self.blastResultsPath  = None
        self.suffix            = None
        self.inputFile         = None
        self.trainFile         = None
        self.testFile          = None
        self.hostTrainPath     = None
        self.hostTestPath      = None
        self.symbTrainPath     = None
        self.symbTestPath      = None
        self.blastSortPath     = None
        self.SVMOutPath        = None
        self.chunkList         = []
        self.threadBlastList   = []


    def getDbPath(self):
        """Return the path to the blast database"""
        if not self.dbPath:
            self.dbPath = os.path.join(self.args.tempDir, DB_NAME)
        return self.dbPath

    def getFastaDbPath(self):
        """Return the path of the fasta file for the blast database"""
        if not self.fastaDbPath:
            self.fastaDbPath = os.path.join(self.args.tempDir, DB_FASTA)
        return self.fastaDbPath

    def getChunkList(self):
        """Return the list of path to the fasta chunks"""
        if not self.chunkList:
            blastInput  = self.args.queries
            fastaPrefix = os.path.basename(blastInput)
            for i in xrange(self.args.nbThreads):
                fastaChunkName = '%s_chunk_%06d' % (fastaPrefix, i)
                fastaChunkName = os.path.join(self.args.tempDir, fastaChunkName)
                self.chunkList.append(fastaChunkName)
        return self.chunkList

    def getThreadBlastList(self):
        """Return the list of output paths of the multi threaded blast"""
        if not self.threadBlastList:
            for i in xrange(self.args.nbThreads):
                blastThreadFile = '%s.%06d' % (BLAST_FILE, i)
                blastThreadPath = os.path.join(self.args.tempDir, blastThreadFile)
                self.threadBlastList.append(blastThreadPath)
        return self.threadBlastList

    def getBlastResultsPath(self):
        """Return the path to the blast results"""
        if self.args.blastResults:
            return self.args.blastResults
        if not self.blastResultsPath:
            self.blastResultsPath = os.path.join(self.args.tempDir, BLAST_FILE)
        return self.blastResultsPath

    def getCheckPointPath(self, dFile):
        """Return the path for the chekpoint (.done) file"""
        return os.path.join(self.args.tempDir, dFile)

    def createCheckPoint(self, cpFile):
        """Create the checkpoint file"""
        path = self.getCheckPointPath(cpFile)
        open(path, 'w').close()

    def checkPoint(self, dFile):
        """Check if a particular checkpoint has been created"""
        path = self.getCheckPointPath(dFile)
        return os.path.exists(path)

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

    def getTrainPath(self):
        """Return the path of the SVM training file"""
        if not self.trainFile:
            fName = self._getSuffix()
            self.trainFile = 'Training' + '_' + fName + '.txt'
        return self.trainFile

    def getTestPath(self):
        """Return the path of the SVM testing file"""
        if not self.testFile:
            fName = self._getSuffix()
            self.testFile = 'Testing' + '_' + fName + '.txt'
        return str(self.testFile)

    def getHostTrainPath(self):
        """Return the path of the host training sequences"""
        if not self.hostTrainPath:
            self.hostTrainPath = os.path.join(self.args.tempDir, HOST_TRAINING)
        return self.hostTrainPath

    def getHostTestPath(self):
        """Return the path of the host testing sequences"""
        if not self.hostTestPath:
            self.hostTestPath = os.path.join(self.args.tempDir, HOST_TESTING)
        return self.hostTestPath

    def getSymbTrainPath(self):
        """Return the path of the symbiont training sequences"""
        if not self.symbTrainPath:
            self.symbTrainPath = os.path.join(self.args.tempDir, SYMB_TRAINING)
        return self.symbTrainPath

    def getSymbTestPath(self):
        """Return the path of the symbiont testing sequences"""
        if not self.symbTestPath:
            self.symbTestPath = os.path.join(self.args.tempDir, SYMB_TESTING)
        return self.symbTestPath

    def getBlastSortPath(self):
        """Get the path of the sequences sorted using blast"""
        if not self.blastSortPath:
            self.blastSortPath = os.path.join(self.args.tempDir, BLAST_SORT)
        return self.blastSortPath

    def getSVMOutPath(self):
        """Get the SVM output path"""
        if not self.SVMOutPath:
            fName = self._getSuffix()
            self.SVMOutPath = fName + '.out'
            self.SVMOutPath = os.path.join(self.args.tempDir, self.SVMOutPath)
        return self.SVMOutPath

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
                yield (name, ''.join(i for i in seq if i.isalpha()))
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(i for i in seq if i.isalpha()))
    handle.close()

def seqCount(path):
    """Counts the number of sequences of a fasta file"""
    c = 0
    if path.endswith('.gz') or path.endswith('.gz"'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        if line.startswith(">"): 
            c += 1
    return c
    
#####################################
#####################################
### Make training set using BLAST ###
#####################################
#####################################

def writeDatabase(args, options, fastaPath):
    """Write Host and Symbiont sequences with standardised names to a new fasta file"""

    logging.info('Creating Database.')
    hostPath      = args.hostSeq
    symbPath      = args.symbSeq
    targetPath = open(fastaPath, 'w')
    #Writing Host Sequences to target database
    i = 0
    for name, seq in iterFasta(hostPath):
        i += 1
        targetPath.write('>%s_%d\n%s\n' % (HOST_NAME, i, seq))
    i = 0
    for name, seq in iterFasta(symbPath):
        i += 1
        targetPath.write('>%s_%d\n%s\n' % (SYMB_NAME, i, seq))
    targetPath.close()
    options.createCheckPoint('writedatabase.done')

def makeDB(args, options):
    """Build the blast database in the temporary folder"""

    dbPath      = options.getDbPath()
    fastaPath   = options.getFastaDbPath()
    logPath     = dbPath + '.log'
    makeblastdb = checkExecutable('makeblastdb')
    makeDBCmd   = [makeblastdb,
                  '-title',
                  DB_NAME,
                  '-in',
                  fastaPath,
                  '-dbtype nucl', ### TODO enable tblastx in the future
                  '-out ',
                  dbPath,
                  '-logfile',
                  logPath]
    makeDBCmd = ' '.join(makeDBCmd)
    submakeDB = subprocess.call(makeDBCmd, shell=True)
    if not submakeDB == 0:
        logging.error('[ERROR] Failed to create blast database')
        sys.exit(1)
    options.createCheckPoint('makeDB.done')

def splitBlastInput(args, options):
    """Split the input fasta file into chunks for parallel blast searches"""

    logging.info('Splitting sequences into %d chunks' % args.nbThreads)
    chunkList = options.getChunkList()
    handles   = []
    for i in xrange(args.nbThreads):
        handle = open(chunkList[i], 'w')
        handles.append(handle)
    #writing to each chunk .fasta
    i = 0
    for name, seq in iterFasta(args.queries):
        handles[i % args.nbThreads].write('>%s\n%s\n' % (name, seq))
        i += 1
    for i in xrange(args.nbThreads):
        handles[i].close()

def runBlast(args, options, threadId):
    """Invoke the blast command. The output format of the result by default is
    set to '6' (tab-seaparated without comments)."""

    #Define BLAST variables
    logging.info('Performing Blast search with thread %d' % threadId)
    eVal         = '%.2e' % args.maxBestEvalue
    dbPath       = options.getDbPath()
    blastOutput  = options.getThreadBlastList()[threadId]
    blast        = checkExecutable(args.blastType)
    blastCmd     = [blast,
                    '-evalue',
                    eVal,
                    '-query',
                    options.getChunkList()[threadId],
                    '-db',
                    dbPath,
                    '-outfmt 6',
                    '-out',
                    blastOutput]
    blastCmd     = ' '.join(blastCmd)
    retCode      = subprocess.call(blastCmd, shell=True)
    if not retCode == 0:
        logging.error('[ERROR] Failed to excecute blast command')
        sys.exit(1)

def mergeBlastOutput(args, options):
    """Merge the output from the blast searches"""

    logging.info('Merging Blast results')
    blastOut       = options.getBlastResultsPath()
    blastOutHandle = open(blastOut, 'w')
    for i in xrange(args.nbThreads):
        threadPath   = options.getThreadBlastList()[i]
        threadHandle = open(threadPath)
        for line in threadHandle:
            blastOutHandle.write(line)
        threadHandle.close()
    blastOutHandle.close()

def runBlastThreads(args, options):
    """Split the queries into chunks, run blast and merge the results"""

    logging.info('Launching threaded Blast search')
    splitBlastInput(args, options)
    threads = []
    for i in xrange(args.nbThreads):
        t = threading.Thread(target=runBlast, args=(args, options, i))
        threads.append(t)
        t.start()
    for i in xrange(args.nbThreads):
        threads[i].join()
    mergeBlastOutput(args, options)
    options.createCheckPoint('runBlast.done')


def parseBlast(args, options):
    """Parse the blast results to be used later to prepare training and testing
    set with unambiguously classified sequences"""

    logging.info('Parsing blast results')
    if not args.blastResults:
        path = options.getBlastResultsPath()
        if path.endswith('.gz'):
            handle = gzip.open(path)
        else:
            handle = open(path,'r')
    else:
        path = options.getBlastResultsPath()
        if path.endswith('.gz'):
            handle = gzip.open(path)
        else:
            handle = open(path,'r')
    querries = {}
    n        = 0
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '#':
            continue
        fields = line.split()
        qName  = fields[0]
        hName  = fields[1]
        evalue = float(fields[10])
        bitscore = float(fields[11])
        if not qName in querries:
            querries[qName] = []
            n += 1
        hit = (hName, evalue, bitscore)
        querries[qName].append(hit)
    logging.info('Parsed %d blast records' % n)
    logging.info('Found %d queries hits' % len(querries))
    handle.close()
    return querries

def classifyFromBlast(querries, args):
    """Classify blast results into ambiguous and unambiguous sequences from
    host and symbiont"""

    def sortHits(h1, h2):
        if h1[1] > h2[1]:
            return 1
        if h1[1] < h2[1]:
            return -1
        return 0

    logging.info('Classifying using Blast results')
    trainingClassification = {}
    blastClassification    = {}
    hostTrained    = 0
    symbTrained    = 0
    hostClassified = 0
    symbClassified = 0
    for qName in querries:
        hits = querries[qName]
        hits.sort(sortHits)
        hasCoral        = False
        hasZoox         = False
        coralBestEvalue = -1
        zooxBestEvalue  = -1
        coralBestBit    = 0
        zooxBestBit     = 0
        for hName, evalue, bitscore in hits:
            if hName.startswith(HOST_NAME):
                if not hasCoral:
                    hasCoral        = True
                    coralBestEvalue = evalue
                    coralBestBit    = bitscore
            else :
                if not hasZoox:
                    hasZoox        = True
                    zooxBestEvalue = evalue
                    zooxBestBit   = bitscore
        if hasCoral and not hasZoox and coralBestEvalue <= args.maxBestEvalue:
            trainingClassification[qName] = HOST_CODE
            blastClassification[qName]    = HOST_CODE
            hostTrained                  += 1
            hostClassified               += 1
        elif hasZoox and not hasCoral and zooxBestEvalue <= args.maxBestEvalue:
            trainingClassification[qName] = SYMB_CODE
            blastClassification[qName]    = SYMB_CODE
            symbTrained                  += 1
            symbClassified               += 1
        if hasZoox and hasCoral:
            bitRatio = float(coralBestBit)/float(zooxBestBit)
            bitDelta = coralBestBit - zooxBestBit
            if bitRatio > MIN_BIT_RATIO and bitDelta > MIN_BIT_DELTA:   
                blastClassification[qName] = HOST_CODE
                hostClassified            += 1
            elif bitRatio <= 0.5 and bitDelta <= -100:
                blastClassification[qName] = SYMB_CODE
                symbClassified            += 1
    logging.info('Found %d unambiguous hits' % len(trainingClassification))
    logging.info('Found %d host only hits' % hostTrained)
    logging.info('Found %d symbiont only hits' % symbTrained)
    logging.info('Found %d likely host hits' % hostClassified)
    logging.info('Found %d likely symbiont hits' % symbClassified)
    return trainingClassification, blastClassification


def trainingSplit(args, options, hCode, sCode):
    """Write the unambiguously classified sequences into four fasta files:
    training.fasta for host sequences, testing.fasta for host sequences,
    training.fasta for symb sequences and testing.fasta for symb sequences."""

    logging.info('Splitting training and testing sequences')
    m = 0
    n = 0
    j = 0
    k = 0
    length = args.numberOfSeq
    handle    = open(args.queries)
    hostTrain = open(options.getHostTrainPath(), 'w')
    hostTest  = open(options.getHostTestPath(), 'w')
    symbTrain = open(options.getSymbTrainPath(), 'w')
    symbTest  = open(options.getSymbTestPath(), 'w')
    blastSort = open(options.getBlastSortPath(), 'w')
    spec1size = seqCount(args.species1) 
    spec2size = seqCount(args.species2)
    rand1List = random.sample(xrange(spec1size), min(length*4,spec1size))
    rand2List = random.sample(xrange(spec2size), min(length*4,spec2size))
    rand1List.sort()
    rand2List.sort()    
    for name, seq in iterFasta(args.species1):
        identity = (name.split(' ')[0])
        if m in rand1List and n < length:
            hostTrain.write('>%s\n%s\n' % (identity, seq))
            n += 1
        else :
            hostTest.write('>%s\n%s\n' % (identity, seq))
        blastSort.write('%s\t%d\n' %(identity, hCode))
        m += 1 
    for name, seq in iterFasta(args.species2):
        identity = (name.split(' ')[0])
        if j in rand2List and k < length:
            symbTrain.write('>%s\n%s\n' % (identity, seq))
            k += 1
        else :
            symbTest.write('>%s\n%s\n' % (identity, seq))
        blastSort.write('%s\t%d\n' %(identity, sCode))
        j += 1    
    
    handle.close()
    hostTest.close()
    hostTrain.close()
    symbTest.close()
    symbTrain.close()
    blastSort.close()
    options.createCheckPoint('parseBlast.done')


def seqSplit(args, options, trainingClassification, blastClassification):
    """Write the unambiguously classified sequences into four fasta files:
    training.fasta for host sequences, testing.fasta for host sequences,
    training.fasta for symb sequences and testing.fasta for symb sequences."""

    logging.info('Splitting training and testing sequences')
    m = 0
    n = 0
    j = 0
    k = 0
    length = args.numberOfSeq
    minSeqSize = args.minSeqSize
    handle    = open(args.queries)
    spec1Path = args.species1
    spec2Path = args.species2
    hostTrain = open(options.getHostTrainPath(), 'w')
    hostTest  = open(options.getHostTestPath(), 'w')
    symbTrain = open(options.getSymbTrainPath(), 'w')
    symbTest  = open(options.getSymbTestPath(), 'w')
    blastSort = open(options.getBlastSortPath(), 'w')
    for name, seq in iterFasta(args.queries):
        size = len(seq)
        identity = (name.split(' ')[0])
        seqClass = trainingClassification.get(identity, 0)
        if seqClass == HOST_CODE:
            if size < minSeqSize:
                n += 1
                continue
            if m < length:
                hostTrain.write('>%s\n%s\n' % (identity, seq))
            else :
                hostTest.write('>%s\n%s\n' % (identity, seq))
            m += 1
    for name, seq in iterFasta(args.queries):
        size = len(seq)
        identity = (name.split(' ')[0])
        seqClass = trainingClassification.get(identity, 0)
        if seqClass == SYMB_CODE:
            if size < minSeqSize:
                k += 1
                continue
            if j < length:
                symbTrain.write('>%s\n%s\n' % (identity, seq))
            else :
                symbTest.write('>%s\n%s\n' % (identity, seq))
            j += 1        
    for blastId in blastClassification:
        blastCode = blastClassification[blastId]
        blastSort.write('%s\t%d\n' % (blastId, blastCode))
    total = n+k
    handle.close()
    hostTest.close()
    hostTrain.close()
    symbTest.close()
    symbTrain.close()
    blastSort.close()
    logging.info('%d training sequences are not parsed due to failure to satisfy minimum seq size requirement.'  %total)
    options.createCheckPoint('parseBlast.done')

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


def computerKmers(args, path, outfile, code, sLength, mode, computeAll):
    """Compute the kmer counts throughout the kmer range for each sequence, and
    write the output to a file.  Each kmer counts will be scaled accordingly
    with the sequence size."""

    logging.info('Computing kmers for %s' % path)
    label = int(code)
    sCounts = seqCount(path)
    if computeAll:
        length = 0
        randList = []
    else:
        length  = args.numberOfSeq
        if length == 0:
            randList = range(0, sCounts)
        else:
            #Check fasta size and create sorted random sequence list        
            randList = random.sample(xrange(sCounts), sLength)
            randList.sort()
        randList = dict.fromkeys(randList)
    # Prepare all maps
    kMin    = args.minWordSize
    kMax    = args.maxWordSize
    maps    = []
    logging.info('Preparing kmer maps')
    for i in xrange(kMin, kMax + 1):
        maps.append(prepareMaps(0, i, ['']))
    # Initialise output
    out     = outfile
    outPath = os.path.join(args.tempDir,out)
    handle  = open(outPath, mode)
    # Initialise counts
    counts  = {}
    for i in xrange(kMin, kMax + 1):
        counts[i] = array.array('d', [0 for x in xrange(4 ** i)])
    # Iterate over sequences
    nSeqs    = 0
    position = 0
    for name, seq in iterFasta(path):
        size = len(seq)
        if length > 0 and nSeqs >= sLength:
            break
        if size < args.maxWordSize:
            position += 1
            continue
        if not position in randList and not computeAll:
            position += 1
            continue
        n      = 0
        handle.write('%d' % label)
        # Obtain Reverse complementary strand
        ###TODO Decide if complementary strand kmers are needed ###
        #seqC = seq[::-1]
        #seqC = seqC.replace('A','3').replace('C','4').replace('G','C').replace('T','A').replace('3','T').replace('4','G')
        # For each kmer value
        for i in xrange(kMin, kMax + 1):
            kCounts  = counts[i]
            # For word in the sequence
            for j in xrange(size - i + 1):
                word  = seq[j:j + i]
                #wordC = seqC[j:j + i]
                kMap  = maps[i - kMin]
                idx   = kMap.get(word,None)
                #idxC  = kMap.get(wordC,None)
                #if idx is None or idxC is None:
                if idx is None:
                    continue
                kCounts[idx]  += 1
                #kCounts[idxC] += 1
            kCountsSum = sum(kCounts)
            for j in xrange(len(kCounts)):
                kCounts[j] /= kCountsSum
            for j in kCounts:
                n += 1
                if j != 0:
                    handle.write(' %d:%.3e' % (n, j))
        handle.write('\n')
        # Reset counts
        for i in xrange(kMin, kMax + 1):
            for j in xrange(len(counts[i])):
                counts[i][j] = 0
        nSeqs    += 1
        position += 1
    # Trace
    logging.info('From %d total sequences, found %d sequences within target range' % (sCounts, nSeqs))
    logging.info('Processed %d sequences' % nSeqs)
    handle.close()

def prepareTrainingKmers(args, options, kmerTrain, kmerTest):
    """Compute the kmer counts for the training and testing sequences.  The
    function outputs two files: a training file and a testing file to be used
    as inputs for the SVM training."""

    logging.info('Preparing kmers for training')
    hostTrainPath = options.getHostTrainPath()
    hostTestPath  = options.getHostTestPath()
    symbTrainPath = options.getSymbTrainPath()
    symbTestPath  = options.getSymbTestPath()
    sLength = min(seqCount(hostTestPath),seqCount(symbTestPath))
    sLength = int(math.floor(sLength/10)*10)
    Length = args.numberOfSeq
    computerKmers(args, hostTrainPath, kmerTrain, HOST_CODE, Length, "w", True)
    computerKmers(args, hostTestPath, kmerTest, HOST_CODE, sLength, "w", False)
    computerKmers(args, symbTrainPath, kmerTrain, SYMB_CODE, Length, "a", True)
    computerKmers(args, symbTestPath, kmerTest, SYMB_CODE, sLength, "a", False)
    options.createCheckPoint('kmers.done')

##################################################################
##################################################################
### SVM computations, based on svm-easy / svm-grid from libsvm ###
##################################################################
##################################################################

def doSVMEasy(args, options, kmerTrain, kmerTest):
    """Scale the input, optimise the SVM parameters and test the predictions.
    the This is roughly the equivalent of the svm-easy script from libsvm."""

    logging.info('Starting SVM training')
    kmerTrain       = os.path.join(args.tempDir, kmerTrain)
    kmerTest        = os.path.join(args.tempDir, kmerTest)
    svmTrain        = checkExecutable('svm-train')
    svmPredict      = checkExecutable('svm-predict')
    svmScale        = checkExecutable('svm-scale')
    scaledFile      = kmerTrain + '.scale'
    modelFile       = kmerTrain + '.model'
    rangeFile       = kmerTrain + '.range'
    scaledTestFile  = kmerTest  + '.scale'
    predictTestFile = kmerTest  + '.predict'
    resultLog       = kmerTrain + '_accuracy.log'
    cmdScale        = [svmScale,
                       '-s',
                       rangeFile,
                       kmerTrain,
                       '>',
                       scaledFile]
    cmdScale        = ' '.join(cmdScale)
    subprocess.call(cmdScale, shell=True)
    c, g, rate      = doSVMGrid(args, options, scaledFile)
    cmdTrain        = [svmTrain,
                       '-c',
                       str(c),
                       '-g',
                       str(g),
                       scaledFile,
                       modelFile]
    cmdTrain        = ' '.join(cmdTrain)
    subprocess.call(cmdTrain, shell=True)
    cmdScale        = [svmScale,
                       '-r',
                       rangeFile,
                       kmerTest,
                       '>',
                       scaledTestFile]
    cmdScale        = ' '.join(cmdScale)
    subprocess.call(cmdScale, shell=True)
    resultHandle    = open(resultLog, 'w')
    cmdPredict      = [svmPredict,
                       scaledTestFile,
                       modelFile,
                       predictTestFile]
    cmdPredict      = ' '.join(cmdPredict)
    subprocess.call(cmdPredict, shell=True, stdout = resultHandle)
    #Adding classification-result to logger
    for line in open(resultLog, 'r'):
        if line.startswith('Accuracy'):
            logging.info('Summary: %s' % line)
    logging.info('Prediction in: %s' % predictTestFile)
    options.createCheckPoint('svm.done')

def calculateSVMGridJobs():
    """Calculate the coordinates of the search space"""

    def rangeF(begin, end, step):
        seq = []
        while True:
            if step > 0 and begin > end:
                break
            if step < 0 and begin < end:
                break
            seq.append(begin)
            begin = begin + step
        return seq

    def permuteSequence(seq):
        n = len(seq)
        if n <= 1:
            return seq
        mid   = int(n / 2)
        left  = permuteSequence(seq[:mid])
        right = permuteSequence(seq[mid+1:])
        ret   = [seq[mid]]
        while left or right:
            if left:
                ret.append(left.pop(0))
            if right:
                ret.append(right.pop(0))
        return ret

    logging.info('Calculating grid coordinates of SVM parameter')
    cSeq = permuteSequence(rangeF(SVM_CSTART, SVM_CEND, SVM_CSTEP))
    gSeq = permuteSequence(rangeF(SVM_GSTART, SVM_GEND, SVM_GSTEP))

    nC   = float(len(cSeq))
    nG   = float(len(gSeq))
    i    = 0
    j    = 0
    jobs = []
    while i < nC or j < nG:
        if i / nC < j / nG:
            # increase C resolution
            line = []
            for k in xrange(0, j):
                line.append((cSeq[i], gSeq[k]))
            i = i + 1
            jobs.append(line)
        else:
            # increase g resolution
            line = []
            for k in xrange(0, i):
                line.append((cSeq[k], gSeq[j]))
            j = j + 1
            jobs.append(line)
    return jobs

class SVMGridWorkerStopToken:
    """Notify the worker to stop"""
    pass

class SVMGridWorker(threading.Thread):
    """Worker thread that calls successive svm-train commands"""

    def __init__(self, name, jobQueue, resultQueue, dataPath):
        threading.Thread.__init__(self)
        self.name        = name
        self.jobQueue    = jobQueue
        self.resultQueue = resultQueue
        self.dataPath    = dataPath

    def run(self):
        """Call svm-train jobs until all the grid coordinates have been explored"""

        while True:
            (c, g) = self.jobQueue.get()
            if c is SVMGridWorkerStopToken:
                self.jobQueue.put((c, g))
                break
            try:
                rate = self.runeOne(2.0 ** c, 2.0 ** g)
                if rate is None:
                    raise RuntimeError(RuntimeError("Got no rate"))
            except:
                # We failed, let others do that and we just quit
                excInfo = sys.exc_info()
                msg     = traceback.format_exception(excInfo[0], excInfo[1], excInfo[2])
                msg     = ''.join(msg)
                logging.warning('[WARNING] Worker %s failed:' % self.name)
                logging.warning(msg)
                self.jobQueue.put((c, g))
                break
            else:
                # TODO do we really need the worker's name?
                self.resultQueue.put((self.name, c, g, rate))

    def runeOne(self, c, g):
        """Call a single svm-train job"""

        svmTrain = checkExecutable('svm-train')
        cmd      = [svmTrain,
                    '-c',
                    str(c),
                    '-g',
                    str(g),
                    '-v',
                    str(SVM_FOLD),
                    self.dataPath]
        cmd      = ' '.join(cmd)
        proc     = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        result   = proc.stdout.readlines()
        for line in result:
            if line.find("Cross") != -1:
                return float(line.split()[-1][:-1])

def doSVMGrid(args, options, dataPath):
    """Search a parameter grid to optimise the SVM parameters.  This is roughly
    equivalent to the svm-grid script from libsvm."""

    logging.info('Optimising SVM parameters')
    jobs        = calculateSVMGridJobs()
    jobQueue    = Queue.Queue(0)
    resultQueue = Queue.Queue(0)
    for line in jobs:
        for (c, g) in line:
            jobQueue.put((c, g))
    jobQueue._put = jobQueue.queue.appendleft

    for i in xrange(args.nbThreads):
        worker = SVMGridWorker('Worker%03d' % i, jobQueue, resultQueue, dataPath)
        worker.start()

    doneJobs   = {}
    svmOutPath = options.getSVMOutPath()
    resultFile = open(svmOutPath, 'w')
    bestRate   = -1
    bestC1     = 0
    bestG1     = 0
    bestC      = 1
    bestG      = 1

    for line in jobs:
        for (c, g) in line:
            while (c, g) not in doneJobs:
                (workerName, c1, g1, rate) = resultQueue.get()
                doneJobs[(c1, g1)]         = rate
                resultFile.write('%f %f %f\n' % (c1, g1, rate))
                if (rate > bestRate) or (rate == bestRate and g1 == bestG1 and c1 < bestC1):
                    bestRate = rate
                    bestC1   = c1
                    bestG1   = g1
                    bestC    = 2.0 ** c1
                    bestG    = 2.0 ** g1
    jobQueue.put((SVMGridWorkerStopToken, None))
    ### TODO check if we need to keep track of the threads and call join
    resultFile.close()
    logging.info('Optimal SVM parameters: c=%f, g=%f, rate=%f' % (bestC, bestG, bestRate))
    return bestC, bestG, bestRate

#Prediction SVM

def loadSVMPredictions(path):
    handle      = open(path)
    content     = handle.read()
    predictions = content.strip().split('\n')
    handle.close()
    return predictions


def loadBlastClassification(options):
    blastSort      = options.getBlastSortPath()
    handle         = open(blastSort)
    classification = {}
    n                   = 0
    for line in handle:
        line = line.strip()
        if not line:
            continue
        fields                 = line.split()
        seqId                  = fields[0]
        seqCode                = fields[1]
        classification[seqId]  = seqCode
        n                     += 1
    logging.info('Parsed %d blast classifications' % n)
    return classification

def writeOutput(args, predictions, blastClassification, fastaPath, fastaName, prefix1, prefix2):
    """Write the final results"""

    logging.info('Writing final output files')
    size        = len(predictions)
    hCode       = str(HOST_CODE)
    sCode       = str(SYMB_CODE)
    blastDict   = blastClassification
    hostResults = prefix1 + '_' + fastaName
    symbResults = prefix2 + '_' + fastaName
    if args.outDir:
        outFolder = os.path.abspath(args.outDir)
        if not os.path.isdir(outFolder):
            os.makedirs(outFolder)
        hostHandle  = open(os.path.join(outFolder,hostResults), "w")
        symbHandle  = open(os.path.join(outFolder,symbResults), "w")
    else:    
        hostHandle  = open(hostResults, "w")
        symbHandle  = open(symbResults, "w")
    j           = 0
    p           = 0
    for name, seq in iterFasta(fastaPath):
        name      = (name.split(' ')[0])
        blastCode = blastDict.get(name, 0)
        if predictions[j] == blastCode:
            if predictions[j] == hCode:
                hostHandle.write('>%s\n%s\n' % (name, seq)) 
            elif predictions[j] == sCode:
                symbHandle.write('>%s\n%s\n' % (name, seq))
        if predictions[j] != blastCode and blastCode != 0:
            p += 1
            if blastCode == hCode:
                hostHandle.write('>%s\n%s\n' % (name, seq))
            elif blastCode == sCode:
                symbHandle.write('>%s\n%s\n' % (name, seq))
        if blastCode == 0:
            if predictions[j] == hCode:
                hostHandle.write('>%s\n%s\n' % (name, seq))
            elif predictions[j] == sCode:
                symbHandle.write('>%s\n%s\n' % (name, seq))
        j += 1
        if j > size:
            logging.warning('[WARNING] Found more sequences than prediction.  This may be caused by dupplicated sequence names. Ignore this warning if training set not obtained through BLAST.')
            break
    if args.species1 and args.species2:
        logging.info('Prediction completed. Sequences are now fully classified; do note that classified sequence of species 1 will be written into %s and species 2 into %s' %(hostResults,symbResults))
    else: 
        logging.info('Found %d contradicting results between blast Classification and SVM prediction.' % p)
    hostHandle.close()
    symbHandle.close()

def predictSVM(args, blastClassification, kmerTrain, kmerTest):
    """Final SVM predictions and combination with the blast results"""

    logging.info('Predicting with SVM optimal parameters')
    svmPredict = checkExecutable('svm-predict')
    svmScale   = checkExecutable('svm-scale')
    kmerTrain  = os.path.join(args.tempDir, kmerTrain)
    modelFile  = kmerTrain + '.model'
    rangeFile  = kmerTrain + '.range'
    fastaPath  = args.queries
    fastaName  = os.path.basename(fastaPath)
    kmerScale  = os.path.join(args.tempDir, fastaName + '.scaled')
    kmerPred   = os.path.join(args.tempDir, fastaName + '.pred')
    kmerFile   = os.path.join(args.tempDir, fastaName + '.kmers')
    sLength    = 0
    computerKmers(args, args.queries, fastaName + '.kmers', HOST_CODE, sLength, "w", True)
    #SVM_Scale
    scaleCmd   = [svmScale,
                  '-r',
                  rangeFile,
                  kmerFile,
                  '>',
                  kmerScale]
    scaleCmd   = ' '.join(scaleCmd)
    retCode    = subprocess.call(scaleCmd, shell=True)
    if not retCode == 0:
        logging.error('[ERROR] Please check inputs. svm-scale not executed or exit with error.')
        sys.exit(1)
    #SVM_predict
    predictCmd = [svmPredict,
                  kmerScale,
                  modelFile,
                  kmerPred]
    predictCmd = ' '.join(predictCmd)
    #subprocess.Popen(predictCmd, shell=True)
    retCode    = subprocess.call(predictCmd, shell=True)
    if not retCode == 0:
        logging.error('[ERROR] Please check inputs. svm-predict not executed or exit with error.')
        sys.exit(1)
    #parse_Prediction
    predictions = loadSVMPredictions(kmerPred)
    writeOutput(args, predictions, blastClassification, fastaPath, fastaName, HOST_NAME, SYMB_NAME)

def tempPathCheck(args):
    """Check if the temporary folder exists, else creates it"""

    logging.info('Checking for temporary folder')
    tempFolder = os.path.abspath(args.tempDir)
    if not os.path.isdir(tempFolder):
        os.makedirs(tempFolder)

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
    parser.add_argument('queries',
                        help='The input queries sequences')
    #subparsers = parser.add_subparsers(help='Choose between option_1 or option_2 input format')
    #group = parser.add_mutually_exclusive_group()
    #parser_1 = subparsers.add_parser('option_1', help='Provide raw protein sequences, and perform blast before preparation for SVM')
    #parser_2 = subparsers.add_parser('option_2', help='Provide blast results as an input, directly start the preparation for SVM')

    parser.add_argument('-A',
                        '--species1',
                        type=str,
                        help='The input first species sequences [unclassified] (single species)')
    parser.add_argument('-B',
                        '--species2',
                        type=str,
                        help='The input 2nd species sequences [unclassified] (single species)')                    
    parser.add_argument('-H',
                        '--hostSeq',
                        type=str,
                        help='The input host sequences (single species)')
    parser.add_argument('-S',
                        '--symbSeq',
                        type=str,
                        help='The input symbiont sequences (single species)')
    parser.add_argument('-b',
                        '--blastResults',
                        type=str,
                        help='Blast results obtained')
    parser.add_argument('-T',
                        '--blastType',
                        type=str,
                        default='blastn',
                        choices=('blastx','blastn','tblastx'),
                        help='Type of blast search to be performed')
    parser.add_argument('-p',
                        '--nbThreads',
                        type=int,
                        default='1',
                        help='Number of threads to run the BLAST search and SVM')
    parser.add_argument('-e',
                        '--maxBestEvalue',
                        type=float,
                        default='1e-20',
                        help='Maximum value for the best e-value')
    ### TODO implement the possibility to have less testing than training
    #parser.add_argument('--trainingTestingRatio',
    #                    type=float,
    #                    default='2',
    #                    help='Value used to divide classfied sequences into testing and training sets')
    parser.add_argument('-n',
                        '--numberOfSeq',
                        type=int,
                        default='0',
                        help='Maximum number of training & testing sequences')
    parser.add_argument('-s',
                        '--minSeqSize',
                        type=int,
                        default='0',
                        help='Minimum sequence size for training & testing sequences')    
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
    parser.add_argument('-v',
                        '--verboseMode',
                        action='store_true',
                        help='Turn Verbose mode on?')
    parser.add_argument('-t',
                        '--tempDir',
                        type=str,
                        default='temp',
                        help='Name of temporary directory')
    parser.add_argument('-o',
                        '--outDir',
                        type=str,
                        default='',
                        help='Name of optional output directory')
    parser.add_argument('-X',
                        '--clearTemp',
                        action='store_true',
                        help='Clear all temporary data upon completion?')
    parser.add_argument('-z',
                        '--stopAfter',
                        type=str,
                        choices=('db','runBlast','parseBlast','kmers','SVM'),
                        help='Optional exit upon completion of stage.')
    parser.add_argument('-R',
                        '--restart',
                        action='store_true',
                        help='Continue process from last exit stage.')
    args = parser.parse_args()
    if args.minWordSize > args.maxWordSize:
        logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
        sys.exit(1)
    return args

def main():
    """Banzai !!!"""

    args    = mainArgs()
    options = PsyTransOptions(args)
    programList = [args.blastType,'svm-train','svm-scale','svm-predict']
    logName = options._getSuffix()
    tempDirName = os.path.basename(args.tempDir)
    if not os.path.basename(args.tempDir):
        tempDirName = tempDirName.split('/')[-2]
    logName = logName + '_' + tempDirName + '_psytrans.log'
    logging.basicConfig(level=logging.DEBUG, format=("%(asctime)s - %(funcName)s - %(message)s"), filename=logName, filemode='w')
    
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(funcName)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    
    restart = args.restart
    
    if not args.blastResults:
        if (not(args.species1 and args.species2) and not (args.hostSeq and args.symbSeq)) or (args.species1 and args.symbSeq) or (args.species2 and args.symbSeq)\
        or (args.species1 and args.hostSeq) or (args.species2 and args.hostSeq):
            logging.error('[ERROR] Either provide the host and symbiont sequences OR the species fasta sequences')
            sys.exit(1)
    else:
        if (args.blastResults and (args.hostSeq or args.symbSeq or args.species1 or args.species2)):
            logging.error('[ERROR] Either provide species fasta sequences OR the host & symbiont sequences OR the output(s) of the blast results')
            sys.exit(1)
    
    if args.verboseMode:
        logging.getLogger().setLevel(logging.DEBUG)
    logging.info("Arguments parsed. Starting...")
    tempPathCheck(args)

    for program in programList:
        if not checkExecutable(program):
            logging.warning('%s program not found. This might cause some problem at the later stage. Please ensure this program is installed correctly.' %program)
            
    # Start from the input sequences
    if not (args.species1 and args.species2):
        if not args.blastResults:
            #Step 1
            fastaPath  = options.getFastaDbPath()
            dbPath     = options.getDbPath()
            if not (restart and options.checkPoint("writeDatabase.done")):
                writeDatabase(args, options, fastaPath)
                if checkExecutable('makeblastdb'):
                    makeDB(args, options)
                else:
                    logging.error('[ERROR] makeblastdb not found. Exiting')
                    sys.exit(1)
            if args.stopAfter == 'db':
                logging.info('Stop after "db" requested, exiting now')
                sys.exit(0)
            #Step 2
            if not (restart and options.checkPoint("runBlast.done")):
                if checkExecutable(args.blastType):
                    runBlastThreads(args, options)
                else:
                    logging.error('[ERROR] blastx not found. Exiting')
                    sys.exit(0)
            if args.stopAfter == 'runBlast':
                logging.info('Stop after "runBlast" requested, exiting now')
                sys.exit(0)
        # Start from the user-provied blast results
        elif not os.path.exists(args.blastResults):
            logging.error('[ERROR] Could not find user-provided blast results (%s). Exiting' % args.blastResults)
            sys.exit(1)

        #Step 3
        if not (restart and options.checkPoint("parseBlast.done")):
            querries = parseBlast(args, options)
            trainingClassification, blastClassification = classifyFromBlast(querries, args)
            seqSplit(args, options, trainingClassification, blastClassification)
    else:
        if not (restart and options.checkPoint("parseBlast.done")):
            trainingSplit(args, options, HOST_CODE, SYMB_CODE)   
        
    if args.stopAfter == 'parseBlast':
        logging.info('Stop after "parseBlast" requested, exiting now')
        sys.exit(0)
    #Step 4
    #Kmer preparation
    kmerTrain = options.getTrainPath()
    kmerTest  = options.getTestPath()
    if not (restart and options.checkPoint("kmers.done")):
        prepareTrainingKmers(args, options, kmerTrain, kmerTest)
    if args.stopAfter == 'kmers':
        logging.info('Stop after "kmers" requested, exiting now')
        sys.exit(0)

    #Step 5
    if not (restart and options.checkPoint("svm.done")):
        if checkExecutable('svm-train') and checkExecutable('svm-scale') and checkExecutable('svm-predict'):  ### TODO look for them only once
            doSVMEasy(args, options, kmerTrain, kmerTest)
        else:
            logging.error('[ERROR] Failed to find some of the libsvm commands. Make sure that svm-train, svm-scale and svm-predict are installed.')
            sys.exit(1)
    if args.stopAfter == 'SVM':
        logging.info('Stop after "SVM" requested, exiting now')
        sys.exit(0)

    #Step 6
    blastClassification = loadBlastClassification(options)
    predictSVM(args, blastClassification, kmerTrain, kmerTest)
    logging.info("SVM classification completed successfully.")

    if args.clearTemp:
        shutil.rmtree(args.tempDir)

if __name__ == '__main__':
    main()

# vim:ts=4:sw=4:sts=4:et:ai:

