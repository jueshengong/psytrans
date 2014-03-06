#!/usr/bin/env python

import argparse
import gzip
import logging
import os.path
import subprocess
import sys
import shutil
import os
import numpy
import threading
import traceback
from math import exp, log, log10

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
LIBSVM_DIR    = 'libsvm'

LETTERS = ('A', 'T', 'G', 'C')

####################################################################
####################################################################
### Class to store the various paths used throughout the program ###
####################################################################
####################################################################

class PsyTransOptions:
    """This class consists of attributes to allow database and file paths to be obtained conveniently.
    Attributes are referred from args, and written to write filenames, filepaths, naming conventions."""

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
        """returns the file path to the blast database"""
        if not self.dbPath:
            self.dbPath = os.path.join(self.args.tempDir, DB_NAME)
        return self.dbPath

    def getFastaDbPath(self):
        """returns the fasta file path of the blast database"""
        if not self.fastaDbPath:
            self.fastaDbPath = os.path.join(self.args.tempDir, DB_FASTA)
        return self.fastaDbPath

    def getChunkList(self):
        """To obtain the list of fasta chunk paths"""
        if not self.chunkList:
            blastInput  = self.args.queries
            fastaPrefix = os.path.basename(blastInput)
            for i in xrange(self.args.nbThreads):
                fastaChunkName = '%s_chunk_%06d' % (fastaPrefix, i)
                fastaChunkName = os.path.join(self.args.tempDir, fastaChunkName)
                self.chunkList.append(fastaChunkName)
        return self.chunkList

    def getThreadBlastList(self):
        """To obtain the list of threaded blast output paths"""
        if not self.threadBlastList:
            for i in xrange(self.args.nbThreads):
                blastThreadFile = '%s.%06d' % (BLAST_FILE, i)
                blastThreadPath = os.path.join(self.args.tempDir, blastThreadFile)
                self.threadBlastList.append(blastThreadPath)
        return self.threadBlastList

    def getBlastResultsPath(self):
        """To obtain the path to the Blast results"""
        if self.args.blastResults:
            return self.args.blastResults
        if not self.blastResultsPath:
            self.blastResultsPath = os.path.join(self.args.tempDir, BLAST_FILE)
        return self.blastResultsPath

    def getCheckPointPath(self, dFile):
        """to obtain the path for the chekpoint (.done) file"""
        return os.path.join(self.args.tempDir, dFile)

    def createCheckPoint(self, cpFile):
        """to create the checkpoint file in the temporary directory, for continuation of the script
        when the option -R is being enabled"""
        path = self.getCheckPointPath(cpFile)
        open(path, 'w').close()

    def checkPoint(self, dFile):
        """ to check if a particular checkpoint has been created"""
        path = self.getCheckPointPath(dFile)
        return os.path.exists(path)

    def _getNumberOfSequences(self):
        """ obtain the length description of the kmer file name"""
        if self.args.numberOfSeq == 0:
            length = 'all'
        else:
            length = self.args.numberOfSeq
        return str(length)

    def _getSuffix(self):
        """Creates the suffix of the SVM input files"""
        if not self.suffix:
            suffix = self._getNumberOfSequences()
            self.mink = str(self.args.minWordSize)
            self.maxk = str(self.args.maxWordSize)
            self.suffix = suffix + '_c' + self.mink + '_k' + self.maxk
        return self.suffix

    def getTrainFile(self):
        """writing the filename and type of the kmer training file as an input for the SVM process"""
        if not self.trainFile:
            fName = self._getSuffix()
            self.trainFile = 'Training' + '_' + fName + '.txt'
        return str(self.trainFile)

    def getTestFile(self):
        """writing the filename and type of the kmer testing file as an input for the SVM process"""
        if not self.testFile:
            fName = self._getSuffix()
            self.testFile = 'Testing' + '_' + fName + '.txt'
        return str(self.testFile)

    def getHostTrainPath(self):
        """get training path of host sequences"""
        if not self.hostTrainPath:
            self.hostTrainPath = os.path.join(self.args.tempDir, HOST_TRAINING)
        return self.hostTrainPath

    def getHostTestPath(self):
        """get testing path of host sequences"""
        if not self.hostTestPath:
            self.hostTestPath = os.path.join(self.args.tempDir, HOST_TESTING)
        return self.hostTestPath

    def getSymbTrainPath(self):
        """get training path of host sequences"""
        if not self.symbTrainPath:
            self.symbTrainPath = os.path.join(self.args.tempDir, SYMB_TRAINING)
        return self.symbTrainPath

    def getSymbTestPath(self):
        """get testing path of host sequences"""
        if not self.symbTestPath:
            self.symbTestPath = os.path.join(self.args.tempDir, SYMB_TESTING)
        return self.symbTestPath

    def getBlastSortPath(self):
        """get testing path of host sequences"""
        if not self.blastSortPath:
            self.blastSortPath = os.path.join(self.args.tempDir, BLAST_SORT)
        return self.blastSortPath

    def getSVMOutPath(self):
        """get svm output path"""
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
                yield (name, ''.join(seq))
            name = line
            seq = []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))
    handle.close()

#####################################
#####################################
### Make training set using BLAST ###
#####################################
#####################################

def writeDatabase(args, options, fastaPath):
    """writes multiple sources to a fasta file. This function also calls iterFasta() in the process."""
    logging.info('Creating Database.')
    hostPath      = args.hostSeq
    symbPath      = args.symbSeq
    targetpath = open(fastaPath, 'w')
    #Writing Host Sequences to target database
    i = 0
    for name, seq in iterFasta(hostPath):
        i += 1
        name = '>' + HOST_NAME +'_' + str(i)
        targetpath.write('%s\n%s\n' % (name, seq))
    j = 0
    for name, seq in iterFasta(symbPath):
        j += 1
        name = '>' + SYMB_NAME +'_' + str(j)
        targetpath.write('%s\n%s\n' % (name, seq))
    targetpath.close()
    options.createCheckPoint('writedatabase.done')

def makeDB(args, options):
    """calls the makeblastdb process from the command line via subprocess module.
    The database files by default will be created and stored at the temporary folder."""
    dbPath    = options.getDbPath()
    fastaPath = options.getFastaDbPath()
    logPath   = dbPath + '.log'
    makeDBCmd = ['makeblastdb',
                 '-title',
                 DB_NAME,
                 '-in',
                 fastaPath,
                 '-dbtype prot', ### TODO enable tblastx in the future
                 '-out ',
                 dbPath,
                 '-logfile',
                 logPath]
    makeDBCmd = ' '.join(makeDBCmd)
    submakeDB = subprocess.call(makeDBCmd, shell=True)
    if not submakeDB == 0:
        logging.warning('[ERROR] Please check logfile. Blast Database not created.')
        sys.exit()
    options.createCheckPoint('makeDB.done')

def splitBlastInput(args, options):
    """To split the input .fasta file into chunks for parallel BLAST search"""
    logging.info('Splitting sequences into %d chunks' % args.nbThreads)
    chunkList = options.getChunkList()
    handles   = []
    for i in xrange(args.nbThreads):
        handle = open(chunkList[i], 'w')
        handles.append(handle)
    #writing to each chunk .fasta
    i = 0
    for name, seq in iterFasta(args.queries):
        handles[i % args.nbThreads].write('%s\n%s\n' % (name, seq))
        i += 1
    for i in xrange(args.nbThreads):
        handles[i].close()

def runBlast(args, options, threadId):
    """calls the blastx process from the command line via subprocess module.
    Result obtained from the BLAST search. The output format of the result by
    default is set to '6' (tab-seaparated without comments)."""
    #Define BLAST variables
    logging.info('Performing Blast search with thread %d' % threadId)
    eVal         = '%.2e' % args.maxBestEvalue
    dbPath       = options.getDbPath()
    blastOutput  = options.getThreadBlastList()[threadId]
    blastCmd = ['blastx',
                '-evalue',
                eVal,
                '-query',
                options.getChunkList()[threadId],
                '-db',
                dbPath,
                '-outfmt 6',
                '-out',
                blastOutput]
    blastCmd = ' '.join(blastCmd)
    retCode  = subprocess.call(blastCmd, shell=True)
    if not retCode == 0:
        logging.error('[ERROR] Please check blastlogfile. Blast search not executed or exit with an error.')
        sys.exit(1)

def mergeBlastOutput(args, options):
    blastOut       = options.getBlastResultsPath()
    logging.info('Merging Blast results into: %s ' % blastOut)
    blastOutHandle = open(blastOut, 'w')
    for i in xrange(args.nbThreads):
        threadPath   = options.getThreadBlastList()[i]
        threadHandle = open(threadPath)
        for line in threadHandle:
            blastOutHandle.write(line)
        threadHandle.close()
    blastOutHandle.close()

def runBlastThreads(args, options):
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
    """parses the blast result given or previously obtained, to later be used to prepare training and testing set of
    unambiguous sequences. This function returns an object called [querries] which summarises the blast results."""
    logging.info('Parsing Blast results')
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
    """classifies the entries in the list [querries] into ambiguous and unambiguous sequences.
    the measurement criteria used here is just the eValue obtaiend from the BLAST results at the moment.
    The function returns a dictionary [classification], with a complete categorisation of ambiguous sequences
    as well as unambiguous (host or symbiont) sequences. This function calls parseBlast() in the process."""

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
    logging.info('Found %d likely  host hits' % hostClassified)
    logging.info('Found %d likely symbiont hits' % symbClassified)
    return trainingClassification, blastClassification


def seqSplit(args, options, trainingClassification, blastClassification):
    """writes the relevant unambiguous sequences into 4 fasta files: training.fasta for host sequences,
    testing.fasta for host sequences, training.fasta for symb sequences and testing.fasta for symb sequences.
    The size/ratio of the training to testing set can be modified by the user from the script option
    (--trainingTestingRatio)."""
    logging.info('Splitting training and testing sequences')
    m = 0
    j = 0
    handle    = open(args.queries)
    hostTrain = open(options.getHostTrainPath(), 'w')
    hostTest  = open(options.getHostTestPath(), 'w')
    symbTrain = open(options.getSymbTrainPath(), 'w')
    symbTest  = open(options.getSymbTestPath(), 'w')
    blastSort = open(options.getBlastSortPath(), 'w')
    for name, seq in iterFasta(args.queries):
        identity = identity = (name.split(' ')[0])[1:]
        seqClass = trainingClassification.get(identity, 0)
        if seqClass == HOST_CODE:
            if m % 2 == 0:
                hostTrain.write('>%s\n%s\n' % (identity, seq))
            else :
                hostTest.write('>%s\n%s\n' % (identity, seq))
            m += 1
        elif seqClass == SYMB_CODE:
            if j % 2 == 0:
                symbTrain.write('>%s\n%s\n' % (identity, seq))
            else :
                symbTest.write('>%s\n%s\n' % (identity, seq))
            j += 1
    for blastId in blastClassification:
        blastCode = blastClassification[blastId]
        blastSort.write('%s\t%d\n' % (blastId, blastCode))
    handle.close()
    hostTest.close()
    hostTrain.close()
    symbTest.close()
    symbTrain.close()
    blastSort.close()
    options.createCheckPoint('parseBlast.done')


############################
############################
### Compute Kmer vectors ###
############################
############################

def prepareMaps(k, maxk, kmers):
    """prepares the kmer maps for the specified kmer range."""
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


def computerKmers(args, path, outfile, code, mode, computeAll):
    """computes the kmer counts throughout the kmer range for each sequence, and write the output to a .txt file.
    Each kmer counts will be scaled accordingly with the sequence size. This function calls prepareMaps() in the process."""
    logging.info('Computing kmers for %s' % path)
    label = int(code)
    if computeAll:
        length = 0
    else:
        length = args.numberOfSeq
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
        counts[i] = numpy.zeros(4 ** i, numpy.double)
    # Iterate over sequences
    nSeqs   = 0
    for name, seq in iterFasta(path):
        if length > 0 and nSeqs >= length:
            break
        #seq   = str(seq.seq).upper()
        size   = len(seq)
        n      = 0
        handle.write('%d' % label)
        # For each kmer value
        for i in xrange(kMin, kMax + 1):
            kCounts  = counts[i]
            # For word in the sequence
            for j in xrange(size - i + 1):
                word = seq[j:j + i]
                kMap = maps[i - kMin]
                idx  = kMap[word]
                kCounts[idx] += 1
            freq = kCounts / kCounts.sum()
            for f in freq:
                n += 1
                if f != 0:
                    handle.write(' %d:%.3e' % (n, f))
        handle.write('\n')
        nSeqs += 1
        # Reset counts
        for i in xrange(kMin, kMax + 1):
            counts[i].fill(0)
    # Trace
    logging.info('Processed %d sequences' % nSeqs)
    handle.close()

def prepareTrainingKmers(args, options, kmerTrain, kmerTest):
    logging.info('Preparing kmers for training')
    """computes the kmer counts for each of the training and testing sequences.
    The function outputs two files:
    a training file and a testing file to be used as inputs for the SVM training."""
    logging.info('Computing kmer counts')
    hostTrainpath = options.getHostTrainPath()
    hostTestpath  = options.getHostTestPath()
    symbTrainpath = options.getSymbTrainPath()
    symbTestpath  = options.getSymbTestPath()
    computerKmers(args, hostTrainpath, kmerTrain, HOST_CODE, "w", False)
    computerKmers(args, hostTestpath, kmerTest, HOST_CODE, "w", False)
    computerKmers(args, symbTrainpath, kmerTrain, SYMB_CODE, "a", False)
    computerKmers(args, symbTestpath, kmerTest, SYMB_CODE, "a", False)
    options.createCheckPoint('kmers.done')

##################################################################
##################################################################
### SVM computations, based on svm-easy / svm-grid from libsvm ###
##################################################################
##################################################################

def doSVMEasy(args, options, kmerTrain, kmerTest):
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
    cmdPredict      = [svmPredict,
                       scaledTestFile,
                       modelFile,
                       predictTestFile]
    cmdPredict      = ' '.join(cmdPredict)
    subprocess.call(cmdPredict, shell=True)
    logging.info('Prediction in: %s' % predictTestFile)
    options.createCheckPoint('svm.done')

def calculateSVMGridJobs():
    """Calculates the coordinates of the search space"""

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

# used to notify the worker to stop
class SVMGridWorkerStopToken:
        pass

class SVMGridWorker(threading.Thread):

    def __init__(self, name, jobQueue, resultQueue, dataPath):
        threading.Thread.__init__(self)
        self.name        = name
        self.jobQueue    = jobQueue
        self.resultQueue = resultQueue
        self.dataPath    = dataPath

    def run(self):
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
    # put jobs in queue
    logging.info('Optimising SVM parameters')
    jobs        = calculateSVMGridJobs()
    jobQueue    = Queue.Queue(0)
    resultQueue = Queue.Queue(0)
    for line in jobs:
        for (c, g) in line:
            jobQueue.put((c, g))
    jobQueue._put = jobQueue.queue.appendleft

    # fire local workers
    for i in xrange(args.nbThreads):
        worker = SVMGridWorker('Worker%03d' % i, jobQueue, resultQueue, dataPath)
        worker.start()

    doneJobs   = {}
    svmOutPath = options.getSVMOutPath()
    resultFile = open(svmOutPath, 'w')
    bestRate   = -1
    bestC1     = None
    bestG1     = None

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
    logging.info('Writing final output files')
    size        = len(predictions)
    hCode       = str(HOST_CODE)
    sCode       = str(SYMB_CODE)
    blastDict   = blastClassification
    hostResults = prefix1 + '_' + fastaName
    symbResults = prefix2 + '_' + fastaName
    hostHandle  = open(hostResults, "w")
    symbHandle  = open(symbResults, "w")
    j           = 0
    for name, seq in iterFasta(fastaPath):
        name      = (name.split(' ')[0])[1:]
        blastCode = blastDict.get(name, 0)
        if predictions[j] == blastCode:
            if predictions[j] == hCode:
                hostHandle.write('>%s\n%s\n' % (name, seq)) 
            elif predictions[j] == sCode:
                symbHandle.write('>%s\n%s\n' % (name, seq))
        if predictions[j] != blastCode and blastCode != 0:
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
        if j >= size:
            logging.warning('[WARNING] Found more sequences than prediction.  This may be caused by dupplicated sequence names.')
            break
    hostHandle.close()
    symbHandle.close()

def predictSVM(args, blastClassification, kmerTrain, kmerTest):
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
    computerKmers(args, args.queries, fastaName + '.kmers', HOST_CODE, "w", True)
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
    """Checks if the temporary folder exists, else creates it."""
    logging.info('Checking for temporary folder')
    tempFolder = os.path.abspath(args.tempDir)
    if not os.path.isdir(tempFolder):
        os.makedirs(tempFolder)

def checkExecutable(program):
    """Checks whether the program is installed and executable."""
    # First check in $PATH
    path = os.getenv('PATH')
    for d in path.split(os.path.pathsep):
        exe = os.path.join(d, program)
        if os.path.exists(exe) and os.access(exe, os.X_OK):
            return exe
    # Then check in the subdirectory
    root = os.path.dirname(os.path.abspath(sys.argv[0]))
    exe  = os.path.join(root, LIBSVM_DIR, program)
    if os.path.exists(exe) and os.access(exe, os.X_OK):
        return exe


def mainArgs():
    parser = argparse.ArgumentParser(description='Performs SVM CLassification of Host-Symbiont Sequences')
    parser.add_argument('queries',
                        help='The input queries sequences')
    #subparsers = parser.add_subparsers(help='Choose between option_1 or option_2 input format')
    #group = parser.add_mutually_exclusive_group()
    #parser_1 = subparsers.add_parser('option_1', help='Provide raw protein sequences, and perform blast before preparation for SVM')
    #parser_2 = subparsers.add_parser('option_2', help='Provide blast results as an input, directly start the preparation for SVM')

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
                        default='1000',
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
    parser.add_argument('-v',
                        '--verboseMode',
                        action='store_true',
                        help='Turn Verbose mode on?')
    parser.add_argument('-t',
                        '--tempDir',
                        type=str,
                        default='temp',
                        help='Name of temporary directory')
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
        sys.stderr.write('Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
        sys.exit(1)
    return args

def main():
    logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
    args    = mainArgs()
    options = PsyTransOptions(args)
    restart = args.restart
    if not (args.hostSeq and args.symbSeq and not args.blastResults) and \
        not (args.blastResults and not (args.hostSeq or args.symbSeq)):
        logging.error('[ERROR] Either provide the host and symbiont sequences OR the output(s) of the blast results')
        sys.exit()
    if args.verboseMode:
        logging.getLogger().setLevel(logging.DEBUG)
    logging.info("Arguments parsed. Starting...")
    tempPathCheck(args)

    # Start from the input sequences
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
            sys.exit()
        #Step 2
        if not (restart and options.checkPoint("runBlast.done")):
            if checkExecutable('blastx'):
                runBlastThreads(args, options)
            else:
                logging.error('[ERROR] blastx not found. Exiting')
                sys.exit(1)
        if args.stopAfter == 'runBlast':
            logging.info('Stop after "runBlast" requested, exiting now')
            sys.exit()
    # Start from the user-provied blast results
    elif not os.path.exists(args.blastResults):
        logging.error('[ERROR] Could not find user-provided blast results (%s). Exiting' % args.blastResults)
        sys.exit(1)

    #Step 3
    if not (restart and options.checkPoint("parseBlast.done")):
        querries = parseBlast(args, options)
        trainingClassification, blastClassification = classifyFromBlast(querries, args)
        seqSplit(args, options, trainingClassification, blastClassification)
    if args.stopAfter == 'parseBlast':
        logging.info('Stop after "parseBlast" requested, exiting now')
        return

    #Step 4
    #Kmer preparation
    kmerTrain = options.getTrainFile()
    kmerTest  = options.getTestFile()
    if not (restart and options.checkPoint("kmers.done")):
        prepareTrainingKmers(args, options, kmerTrain, kmerTest)
    if args.stopAfter == 'kmers':
        logging.info('Stop after "kmers" requested, exiting now')
        return

    #Step 5
    if not (restart and options.checkPoint("svm.done")):
        if checkExecutable('svm-train') and checkExecutable('svm-scale') and checkExecutable('svm-predict'):  ### FIXME look for them only once
            doSVMEasy(args, options, kmerTrain, kmerTest)
        else:
            logging.error('[ERROR] Failed to find some of the libsvm commands. Make sure that svm-train, svm-scale and svm-predict are installed.')
            sys.exit(1)
    if args.stopAfter == 'SVM':
        logging.info('Stop after "SVM" requested, exiting now')
        return

    #Step 6

    blastClassification = loadBlastClassification(options)
    predictSVM(args, blastClassification, kmerTrain, kmerTest)
    logging.info("SVM classification completed successfully.")

    if args.clearTemp:
        shutil.rmtree(args.tempDir)

if __name__ == '__main__':
    main()

# vim:ts=4:sw=4:sts=4:et:ai:
