#!/usr/bin/env python

import argparse
import array
import gzip
import logging
import os.path
import sys
import string
import tempfile
import numpy
import random
import math
from pypr.clustering import *
import scipy.ndimage
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.mlab import PCA


########################
########################
### Global constants ###
########################
########################

HOST_CODE  = 1
SYMB_CODE  = 2
NBINS      = 80   ###NOT ADVISABLE TO SET NBIN LOWER THAN 50
CONF_LEVEL = 0.95
THRESHOLD  = 2
PC_DIM     = 2   # Total principal components to retain, advised to be not more than 4
PC_1       = 1   # Order of principal component to compare: 1 = PC1, etc
PC_2       = 2
N_CLASSIFY = 2   # Do Not Change. Only for multiple classification. Not supported atm.

#Parameters to optimize Gaussian Model Fitting
MAX_ITER  = 50
TRIALS    = 10
TRIM      = 0

#Binned Matrix files (reference)
RAW_MATRIX      = 'rawMatrix.txt' # Binned PCA matrix in actual counts
AFTER_ERODE     = 'AfterErode.txt'
BIN_MATRIX      = '01Matrix.txt'
FILTERED_MATRIX = 'filtered.txt'


LETTERS = ('A', 'T', 'G', 'C')

COMPLEMENT_TABLE = string.maketrans('ATGCatgc', 'TACGtacg')


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
        self.inputFile         = None
        self.logName           = None
        self.plotPath          = None
        self.spec1Path         = None
        self.spec2Path         = None
        self.suffix            = None
        self.tempDir           = None
        self.createTempDir()


    def _getNumberOfSequences(self):
        """Return the length part of the kmer file name"""
        if self.args.numberOfSeq == 0:
            #effective set size is capped by min(spec1train, spec2train)
            length = 'max'  
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
    
    def createTempDir(self):
        if not os.path.exists(self.args.tempDir):
            os.mkdir(self.args.tempDir)
        self.tempDir = tempfile.mkdtemp(prefix='autotrain_',
                                        suffix='_temp',
                                        dir=self.args.tempDir)
    def getInputFile(self):
        if not self.inputFile:
            self.inputFile = self.args.queries
        return self.inputFile
                                                
    def getLogPath(self):
        if not self.logName:
            self.logName = 'autotrain_' + self._getSuffix() + '.log'
        return self.logName
        
    def getPlotPath(self):
        if not self.plotPath:
            if self.args.outputPlotName:
                self.plotPath = self.args.outputPlotName
            else:
                end_str       = '_PC' + str(PC_1) + '_v_PC' + str(PC_2)
                self.plotPath = 'autotrain_' + self._getSuffix() + end_str
        return self.plotPath

    def getSpec1Path(self):
        if not self.spec1Path:
            self.spec1Path = str(self.args.spec1Path) + '.fasta'
        return self.spec1Path

    def getSpec2Path(self):
        if not self.spec2Path:
            self.spec2Path = str(self.args.spec2Path) + '.fasta'
        return self.spec2Path
        

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

def seqCount(path):
    """Counts the number of sequences of a fasta file"""
    c = 0
    if path.endswith('.gz'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        if line.startswith(">"): 
            c += 1
    return c

def writeMat(options, filename, matrix):
    pathname = os.path.join(options.tempDir, filename)
    handle = open(pathname, 'w')
    for i in xrange(len(matrix)):
        line = []
        for j in xrange(len(matrix[i])):
            entry = str(int(matrix[i, j]))
            line.append(entry)
        handle.write("%s \n" %line)
    handle.close()

######################################
####STAT DISTRIBUTIONS###############
#####################################


def multiGaussianPdf(x, mu, sigma):
    """Function to map data-point to a multivariate gaussian distribution"""
    size = len(x)
    sigma = numpy.matrix(sigma)
    if size == len(mu) and (size, size) == sigma.shape:
        det = numpy.linalg.det(sigma)
        if det == 0:
            raise NameError("The covariance matrix can't be singular")
        constant = 1.0/ ( math.pow((2*math.pi),float(size)/2) * math.pow(det,1.0/2) )
        newMu = numpy.matrix(x - mu)
        inv = sigma.I        
        expo = math.pow(math.e, -0.5 * (newMu * inv * newMu.T))
        return constant * expo
    else:
        raise NameError("The dimensions of the input don't match")


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



def computeKmerMat(options, full=False):
    """Compute the kmer counts throughout the kmer range for each sequence, and
    write the output to a file.  Each kmer counts will be scaled accordingly
    with the sequence size."""
    path = options.getInputFile()
    logging.info('Computing kmers for %s' % path)
    sCounts = seqCount(path)
    if full:
        length = 0
        size = sCounts
        randList = range(0, sCounts) #indexing the list
    else:
        length = options.args.numberOfSeq
        size = length      
        randList = random.sample(xrange(sCounts), length) #indexing the list used for training PC projection
        randList.sort()
    randDict = dict.fromkeys(randList)
    # Prepare all maps
    kMin    = options.args.minWordSize
    kMax    = options.args.maxWordSize
    maps    = []
    logging.info('Preparing kmer maps')
    for i in xrange(kMin, kMax + 1):
        maps.append(prepareMaps(0, i, ['']))	    
    # Initialise counts
    counts     = {}
    kmerArray = []
    kmerSize = 0
    for i in xrange(kMin, kMax + 1):
        counts[i] = array.array('d', [0 for x in xrange(4 ** i)])
    	kmerRange = 4 ** i
    	kmerSize  = kmerSize + kmerRange
    #Initialise data Matrix	
    kmerMatrix = numpy.zeros(shape=(size,kmerSize))
    # Iterate over sequences
    nSeqs    = 0
    position = 0
    for name, seq in iterFasta(path):
        size = len(seq)
        if length > 0 and nSeqs >= length:
            break
        if not position in randDict:
            position += 1
            continue
        n      = 0
        seqs = (seq, )
        if options.args.bothStrands:
            seqs = (seq, revComp(seq))
        # For each kmer value
        for i in xrange(kMin, kMax + 1):
            kCounts  = counts[i]
            # For word in the sequence
            for j in xrange(size - i + 1):
                word  = seq[j:j + i]
                kMap  = maps[i - kMin]
                idx   = kMap.get(word,None)
                if idx is None:
                    continue
                kCounts[idx]  += 1
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
    ##NOTE: randList is needed later on to act as (randList order - seq order) reference that maps back to the queried sequence
    return kmerMatrix, randList

##########################################
########### RUN PCA ######################
##########################################


def doPCA(dataMatrix):
    """ Running the PCA module from matplotlib.mlab which returns the pca object"""
    logging.info('Computing PCA')
    pcaObject = PCA(dataMatrix)
    logging.info("Successfully completed PCA")
    return pcaObject

def projectMat(pca, matrix, dim):
    """perform the pca projection to the data matrix based on input pca object"""
    projected = numpy.zeros(shape=(len(matrix),dim))
    logging.info('Writing PC projections to data matrix')
    for i in xrange(len(matrix)):
        projected[i] = pca.project(matrix[i,:])[:dim]
    return projected
        	    
def runPCA(trainMatrix, dim):
    """Wrapper function to run PCA, and return both the full projected PCs on the entire data"""
    pcaObject = doPCA(trainMatrix)
    projected = projectMat(pcaObject, trainMatrix, dim)
    return projected
    
def plotDensityPCA(options, pcaMatrix):
    """ Plotting the PCA contour plot"""
    plotPath = options.getPlotPath()
    axisX = PC_1 - 1 #First PC correspond to entry [0]
    axisY = PC_2 - 1
    pcaXCoord = pcaMatrix[:,axisX]
    pcaXCoord = pcaXCoord.real
    pcaYCoord = pcaMatrix[:,axisY]
    pcaYCoord = pcaYCoord.real
    H, edgeX, edgeY = numpy.histogram2d(pcaXCoord, pcaYCoord, bins = NBINS)
    H = numpy.rot90(H)
    H = numpy.flipud(H)
    # mask zeroes
    maskedH = numpy.ma.masked_where(H==0, H)
    #Plot the histogram
    fig2 = plt.figure()
    plt.pcolormesh(edgeX, edgeY, maskedH)
    plt.xlabel('Pricinpal Component %d' % PC_1)
    plt.ylabel('Principal Component %d' % PC_2) 
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    fig2.savefig(plotPath, format=options.args.format)
    axisData = plt.axis()
    return axisData

#####################################################
######CLASSIFICATION VIA POSTERIOR LIKELIHOOD#######
####################################################

def clusterViaPosterior(dp, centre, cov, weight, trim):
    '''Function to cluster a given data point to a Gaussian Component based on Posterior likelihood'''   
    def calcPosteriorProb(dp, mu, sigma, weight):
        postProb = multiGaussianPdf(dp, mu, sigma)*weight
        return postProb
    
    def pDelta(entries, pos):
        deltas = []
        for i in entries:
            d = numpy.abs(i - entries[pos])
            deltas.append(d)
        return max(deltas)
                
    size = len(weight)
    posteriors = []
    for j in xrange(len(weight)):
        mu    = centre[j]
        sigma = cov[j]
        prob = calcPosteriorProb(dp, mu, sigma, weight[j])
        posteriors.append(prob)
    normPosteriors = [v/sum(posteriors) for v in posteriors]
    for i in xrange(len(normPosteriors)):
        if pDelta(normPosteriors, i) >= float(1-trim)/size and normPosteriors[i] == max(normPosteriors):
            label = i + 1
            return label
        else:
            label = -1
    return label
    

def labelClustersViaPost(dataMatrix, centre, cov, weight, trim):
    '''label each data point based on Posterior likelihood of Gaussian Components. Returns a dictionary object containing label.'''
    length = len(dataMatrix)
    clusterDict = {}
    for i in xrange(length):
        dp      = (dataMatrix[i,0], dataMatrix[i,1])
        dpGroup = clusterViaPosterior(dp, centre, cov, weight, trim)
        clusterDict[i] = dpGroup
    return clusterDict


def fitGMM(dataMatrix, maxIter, trials, trim):
    '''Wrapper function for posterior GMM clustering'''
    nClassify = N_CLASSIFY  
    logging.info("Optimizing GMM parameters based on Expectation Maximisation")
    centre, cov, weight, lg = gmm.em(dataMatrix, nClassify, max_iter=maxIter, verbose=False, iter_call=None, delta_stop=9.9999999999999995e-07, init_kw={}, max_tries=trials, diag_add=0.001)
    logging.info("Model obtained. log-likelihood of fitted model: %d" %lg)
    clusterDict = labelClustersViaPost(dataMatrix, centre, cov, weight, trim)
    return clusterDict

     

#########################################
########## BIN FUNCTIONS#################
#########################################

def slotTag(X, Y):
    """A simple function to provide a name-tag for the bin dictionary entries"""
    slotId = str(int(X)) + '_' + str(int(Y))
    return slotId

def get_axis(matrix):
    '''A function to obtain the axis limits for the Bin-matrix'''
    minX = math.floor(min(matrix[:,0])/10)*10
    maxX = math.floor(max(matrix[:,1])/10)*10
    minY = math.floor(min(matrix[:,0])/10)*10
    maxY = math.floor(max(matrix[:,1])/10)*10 
    axisData = [minX, maxX, minY, maxY]
    return axisData

def trimMatrix(matrix, confidence):
    """Function to trim the tail distribution of the matrix data""" 
    size = len(matrix)
    if confidence == 1:
        return matrix
    elif confidence == 0:
        logging.error("confidence level have to be a value between 0 and 1")
        sys.exit(1)
    trim = math.ceil(size*(1-confidence)/2)
    ind = numpy.lexsort((matrix[:,1],matrix[:,0]))
    sortedMat = matrix[ind]
    sortedMat = sortedMat[0+trim-1:-trim]
    #keep track of the unsorted sequence id
    origIndex     = ind[0+trim-1:-trim]
    return sortedMat, origIndex

def calcBins(options, matrix, origIndex, axisData, nbinX, nbinY, margin):
    """Function to convert the projected PCA matrix data into Bin counts"""
    binDict= {}
    minX = float(axisData[0])
    maxX = float(axisData[1])
    minY = float(axisData[2])
    maxY = float(axisData[3])
    stepX = (maxX - minX)/nbinX
    stepY = (maxY - minY)/nbinY
    binMatrix = numpy.zeros(shape=(nbinX,nbinY))
    for i in xrange(len(matrix)):
        xpt = matrix[i, 0]
        ypt = matrix[i, 1]
        slotX = math.ceil((xpt - minX)/stepX) - 1
        slotY = math.ceil((ypt - minY)/stepY) - 1
        if not (0 <= slotX < (nbinX - 1)) or not (0 <= slotY < (nbinY - 1)):
            continue
        binMatrix[slotX, slotY] += 1
        slotId = slotTag(slotX, slotY)
        #append original sequence Id to dictionary for the correct bin slot
        binDict.setdefault(slotId, []).append(origIndex[i])
    writeMat(options, RAW_MATRIX, binMatrix) 
    return binMatrix, binDict

def posteriorErode(untrimmedData, binMatrix, clusterDict, origIndex, axisData, nbinX, nbinY):
    '''Perform primary erosion and cluster using posterior porbability cut-offs'''                  
    matP = untrimmedData[origIndex]
    minX = float(axisData[0])
    maxX = float(axisData[1])
    minY = float(axisData[2])
    maxY = float(axisData[3])
    stepX = (maxX - minX)/nbinX
    stepY = (maxY - minY)/nbinY
    for i in xrange(len(matP)):
        xpt = matP[i, 0]
        ypt = matP[i, 1]
        slotX = math.ceil((xpt - minX)/stepX) - 1
        slotY = math.ceil((ypt - minY)/stepY) - 1
        if not (0 <= slotX < (nbinX - 1)) or not (0 <= slotY < (nbinY - 1)):
            continue
        if binMatrix[slotX, slotY] == 0:
            continue
        if clusterDict.get(origIndex[i], 0) == -1:
            binMatrix[slotX, slotY] = 0
                
    return binMatrix
    
        
def thresholding(matrix, percentile):
    """To obtain the appropriate value to as a threshold, based on percentile of counts, before performing erosion."""
    countlist = []
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix[0])):
            count = matrix[i, j]
            if not count == 0:
                countlist.append(count)
    val = numpy.percentile(countlist, percentile)
    return val

########################################
####NEIGHBOURHOOD FUNCTIONS#############
########################################
    
def neighbourhood(matrix,i,j):
    """This function returns the neighbourhood (1 cell difference) of a specific matrix cell."""
    neighbour = []
    for binX in range(i-1, i+2):
        for binY in range(j-1, j+2):
            if -1<i<=len(matrix)  and -1<j<=len(matrix[0]) and (i != binX or j != binY) \
            and 0<=binX<len(matrix) and 0<=binY<len(matrix[0]):
                entry = (binX,binY)
                neighbour.append(entry)
    return neighbour   

def hasZeroAround(matrix, i, j):
    """This function checks whether there is a [0] valued neighbour on a specific matrix entry [i,j]"""
    isZero = False
    for coord in neighbourhood(matrix, i, j):
        if matrix[coord] == 0:
            isZero = True
            return isZero
        else:
            continue
    return isZero 

def convertNeighbour(origin, i ,j):
    """This function converts the non-zero/unclassified neighbourhood of a matrix entry [i,j] into the value of the matrix mat[i,j]=val itself, and updates the matrix object"""
    matrix = numpy.copy(origin)
    val = matrix[i, j] 
    for coord in neighbourhood(matrix, i, j):
        if matrix[coord] == 0 or matrix[coord] == val:
            continue
        else:
            matrix[coord] = val
    return matrix

def seedLabel(origin, unspecified, label):
    """this function seeds the unspecified region of the matrix, with a predefined label"""
    matrix = numpy.copy(origin)
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix[i])):
            if matrix[i,j] == int(unspecified):
                matrix[i,j] = label
                return matrix
    
def checkRegion(matrix, label):
    """This function checks whether entries with value [label] exist in the matrix"""
    logging.debug("checkin region for matrix on label %d" %label)
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix[i])):
            if matrix[i,j] == int(label):
                return True
            else:    
                continue
    return False
                
def convertRegion(matrix, label):
    """This function converts the corresponding region of a matrix connected to a specific label"""
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix[i])):
            if matrix[i,j] == int(label):
                matrix = convertNeighbour(matrix, i ,j)
            else:    
                continue
    return matrix   

def getRegionSize(matrix, label):
    """This function returns the size of the region of a specific label"""
    count = 0
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix[i])):
            if matrix[i,j] == label:
                count += 1
            else:    
                continue
    return count
                
#####################################
######## EROSION FUNCTIONS###########
#####################################

def manualErode(options, matrix):
    """This function performs a 1-step erosion on the bin-matrix"""
    original = matrix
    matB = numpy.copy(matrix)
    erodeList = []
    for i in xrange(len(original)):
        for j in xrange(len(original[i])):
            if original[i,j] == 0:
                continue    
            if hasZeroAround(original, i, j):
                matB[i,j] = 0
    writeMat(options, AFTER_ERODE, matB)
    return matB

def labelRegions(matrix):
    """The wrapper function to perform labelling and conversion of regions, until all regions have been fully classified."""
    mat = numpy.copy(matrix)
    label = 2
    unspecified = 1
    clusters = 0
    while checkRegion(mat, unspecified):
        mat = seedLabel(mat, unspecified, label)
        mat = convertRegion(mat, label)
        label += 1
        clusters += 1
    labelled = label -1 
    logging.info("found %d regions" %clusters)
    return label, mat
    
def filterRegions(regions, origin):
    """Based on the number of species to be classified, this function returns the largest region to be used for classification if multiple 
    clusters exists, and returns the bin-matrix with only nClassify=2 regions, and a list of the corresponding region labels """
    matrix = numpy.copy(origin)
    nspecies = N_CLASSIFY
    clusters = 0
    regionList  = range(0, regions+1)
    regionCount = [0,0]
    for i in regionList:
        if i<2:
            continue
        area = getRegionSize(matrix, i)
        regionCount.append(area)
    logging.debug("regionCounts [0,1,2,...]")
    logging.debug(regionCount)    
    nspecList = numpy.argsort(regionCount)[-nspecies:]
    if not nspecies == (regions-1):
        for i in regionCount:
            if i >= 9:
                clusters += 1
        if clusters > nspecies:
            logging.warning('Found multiple clustered regions. Please ensure dataset is clean. Results might be inaccurate.') 
    for i in xrange(len(matrix)):
        for j in xrange(len(matrix[i])):
            if matrix[i,j] in nspecList:
                continue
            else:
                matrix[i,j] = 0
    return matrix, nspecList
    
    
def erodeMap(options, matrix, threshold):
    """This is a wrapper function to perform all the erosion, filtering of regions, and classification 
    of clusters. The function require the raw Bin-matrix with original counts and a threshold parameter to be specified.
    The function returns the bin-matrix with only nClassify=2 regions, and a list of the corresponding region labels."""
    ##Filter via threshold##
    for i in xrange(len(matrix)):
        if not 0 < i < (len(matrix) - 2):
            continue
        for j in xrange(len(matrix[i])):
            if not 0 < j < (len(matrix[i]) - 2):
                continue
            if 0 < matrix[i,j] < threshold:
                matrix[i,j] = 0   
    erodeMat = scipy.sign(matrix)
    #Raw Binary Matrix file in txt format (in tmp)
    writeMat(options, BIN_MATRIX, erodeMat)
    #perform binary erosion via iteration
    eroded = manualErode(options, erodeMat)
    regions, eroded = labelRegions(eroded)
    filtered, nspecList = filterRegions(regions, eroded)
    writeMat(options, FILTERED_MATRIX, filtered)                                           
    return filtered, nspecList


def filterTrainSeq(filteredMat, binDict, nspecList):
    """This function retrieve sequences in the final classified bin-matrix regions, and parse them in lists."""
    hCode = int(nspecList[0])
    sCode = int(nspecList[1])
    logging.debug("hCode, is %d sCode is %d" %(hCode, sCode))
    hostTrainList = []
    symbTrainList = []
    for i in xrange(len(filteredMat)):
        if sum(filteredMat[i]) == 0:
            continue
        for j in xrange(len(filteredMat[i])):
            if filteredMat[i,j] == 0:
                continue
            if filteredMat[i,j] == hCode:
                slotId = slotTag(i,j)
                entry  = binDict.get(slotId,[])
                hostTrainList.extend(entry)
            elif filteredMat[i,j] == sCode:
                slotId = slotTag(i,j)
                entry  = binDict.get(slotId,[])
                symbTrainList.extend(entry) 
    return hostTrainList, symbTrainList


def writeTrainFile(options, hTrainList, sTrainList, randList):
    """This function writes the corresponding sequence from the training lists into fasta files"""
    spec1list = []
    spec2list = []
    limit = min(len(hTrainList),len(sTrainList)) #to ensure same training size
    for i in hTrainList:
        origSeq = [randList[i]]
        spec1list.extend(origSeq)
    for i in sTrainList:
        origSeq = [randList[i]]
        spec2list.extend(origSeq)
    spec1list = dict.fromkeys(spec1list)
    spec2list = dict.fromkeys(spec2list)
    path = options.getInputFile()
    handle1 = open(options.getSpec1Path(), 'w')
    handle2 = open(options.getSpec2Path(), 'w')
    k = 0
    hlim = 0
    slim = 0
    for name, seq in iterFasta(path):
        if hlim >= limit and slim >= limit:
            break 
        if k in spec1list and hlim < limit:
            handle1.write('>%s\n%s\n' % (name, seq))
            k += 1
            hlim += 1
        elif k in spec2list and slim < limit:
            handle2.write('>%s\n%s\n' % (name, seq))
            k += 1
            slim += 1
        else:
            k += 1
            continue
    logging.info('Training set of %d sequences have been prepared' %limit)
    handle1.close()
    handle2.close()
    return hTrainList, sTrainList
    
    
def clusterBins(options, matrix, confidence, percentile, clusterDict):
    """Wrapper function to obtain the eroded bin-matrix, and then the filtered bin-matrix."""  
    logging.info('Morphing bin-maps via erosion') 
    axisData = get_axis(matrix)
    sortedMatrix, origIndex = trimMatrix(matrix, confidence)
    binMatrix, binDict = calcBins(options, sortedMatrix, origIndex, axisData, NBINS, NBINS, 1000)
    binMatrix = posteriorErode(matrix, binMatrix, clusterDict, origIndex, axisData, NBINS, NBINS)
    threshold = thresholding(binMatrix, percentile)
    filteredMat, nspecList = erodeMap(options, binMatrix, threshold)   
    return filteredMat, binDict, nspecList


def autoTrainPCA(options, confidence, matrix, randList, percentile):
    """Wrapper function for the entire PCA auto-training process"""
    clusterDict = fitGMM(matrix, 50, 10, 0)
    finalMatrix, binDict, nspecList = clusterBins(options, matrix, confidence, percentile, clusterDict) 
    hTrainList, sTrainList = filterTrainSeq(finalMatrix, binDict, nspecList)
    hTrainList, sTrainList = writeTrainFile(options, hTrainList, sTrainList, randList)
                

def mainArgs():
    """Process command-line arguments"""

    parser = argparse.ArgumentParser(description='Perform auto-training of queried fasta sequences')
    parser.add_argument('queries',
                        help='The query fasta sequences')
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
    parser.add_argument('-p',
                        '--spec1Path',
                        type=str,
                        default='spec1',
                        help='Output Training filename of the 1st species')
    parser.add_argument('-q',
                        '--spec2Path',
                        type=str,
                        default='spec2',
                        help='Output Training filename of the 2nd species')  
    parser.add_argument('-r',
                        '--bothStrands',
                        action='store_true',
                        help='Compute kmers for the forward and reverse strands')
    parser.add_argument('-t',
                        '--tempDir',
                        type=str,
                        default='temp',
                        help='Location (prefix) of the temporary directory')                                  
    parser.add_argument('-o',
                        '--outputPlotName',
                        type=str,
                        default='',
                        help='Name of the output PCA plot')
    parser.add_argument('-f',
                        '--format',
                        type=str,
                        default='png',
                        help='Format of the output PCA plot')
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='Turn Verbose mode on?')
    args = parser.parse_args()
    if args.minWordSize > args.maxWordSize:
        logging.error('[ERROR] Minimum kmer size (-c/--minKmerSize) must be less than Maximum kmer size (-k/--maxKmerSize)\n')
        sys.exit(1)
    return args


def main():
    """Autotrain Banzai!!!"""

    #Set Seed
    random.seed(122)

    #Get options
    args    = mainArgs()
    options = pcaOptions(args)

    #Setting up script logging to .log file and console
    logFormat = "%(asctime)s - %(funcName)s - %(message)s"
    logging.basicConfig(level=logging.INFO,
                        format=logFormat,
                        filename=options.getLogPath(), filemode="w")
    console   = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter(logFormat)
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    
    #Initialize process by computing kmer vectors for the sample
    logging.info("Arguments parsed. Starting auto-train...")
    trainMatrix, sampledList = computeKmerMat(options, full=True)    
    logging.info('kmer Matrix created. preparing PCA') 
    
    #Extract PC projections from a subsample, and project to the whole dataset
    projected = runPCA(trainMatrix, PC_DIM)
    
    #Plots PC1 against PC2
    plotDensityPCA(options, projected)
    logging.info('performing PCA Auto-training') 
    
    #Performs the auto-training: GMM fitting, Posterior likelihood filtering, Erosion, Thresholding
    autoTrainPCA(options, CONF_LEVEL, projected, sampledList, 50)
    logging.info("Auto-Training completed successfully.")

if __name__ == '__main__':
    main()
# vim:ts=4:sw=4:sts=4:et:ai:
