#!/usr/bin/env python

import argparse
import array
import gzip
import logging
import os.path
import sys
import numpy
import random
import math
from pypr.clustering import *
from matplotlib.mlab import PCA
from numpy.testing import assert_array_almost_equal
import statsmodels.api as sm
import scipy.ndimage
#from statsmodels.sandbox.tools import pca
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from matplotlib.mlab import PCA

#from statsmodels.sandbox.tools.cross_val import LeaveOneOut



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

BINARIES_DIR    = 'binaries'

LETTERS = ('A', 'T', 'G', 'C')

#random.seed(123)
random.seed(122)

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

def writeMat(filename, matrix):
    handle = open(filename, 'w')
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
    size = len(x)
    sigma = numpy.matrix(sigma)
    if size == len(mu) and (size, size) == sigma.shape:
        det = numpy.linalg.det(sigma)
        if det == 0:
            raise NameError("The covariance matrix can't be singular")
        constant = 1.0/ ( math.pow((2*math.pi),float(size)/2) * math.pow(det,1.0/2) )
        newMu = numpy.matrix(x - mu)
        inv = sigma.I        
        exp = math.pow(math.e, -0.5 * (newMu * inv * newMu.T))
        return constant * exp
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



def computeKmerMat(args, full=False):
    """Compute the kmer counts throughout the kmer range for each sequence, and
    write the output to a file.  Each kmer counts will be scaled accordingly
    with the sequence size."""
    path = args.queries
    logging.info('Computing kmers for %s' % path)
    sCounts = seqCount(path)
    if full:
        length = 0
        size = sCounts
    else:
        length  = args.numberOfSeq
        size = length
    if length == 0:    
        randList = range(0, sCounts)
    else:
        #Check fasta size and create sorted random sequence list        
        randList = random.sample(xrange(sCounts), length)
        randList.sort()
    randDict = dict.fromkeys(randList)
    # Prepare all maps
    kMin    = args.minWordSize
    kMax    = args.maxWordSize
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
        #handle.write('%d' % label)
        # For each kmer value
        for i in xrange(kMin, kMax + 1):
            kCounts  = counts[i]
            # For word in the sequence
            for j in xrange(size - i + 1):
                word = seq[j:j + i]
                kMap = maps[i - kMin]
                idx  = kMap.get(word,None)
                if idx is None:
                    continue
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
    return kmerMatrix, randList

##########################################
########### RUN PCA ######################
##########################################


def doPCA(dataMatrix):
    """ Running the PCA module from matplotlib.mlab which returns the pca object"""
    logging.info('Computing PCA')
    pcaObject = PCA(dataMatrix)
    #a, b, c, d = pca(dataMatrix, keepdim=2, normalize=True, demean=True)
    logging.info("Successfully completed PCA")
    return pcaObject


def plotDensityPCA(args, pcaMatrix, axis1, axis2, graphCode):
    """ Plotting the PCA contour plot"""
    outName = args.outputFileName
    PCPrefix = 'PC' + str(axis1) + 'and PC' + str(axis2) + '_' + str(graphCode)
    outName = outName + '_PCAdensity_' + PCPrefix + '_' + str(NBINS)
    axisX = axis1 - 1
    axisY = axis2 - 1
    pcaXCoord = pcaMatrix[:,axisX]
    pcaXCoord = pcaXCoord.real
    pcaYCoord = pcaMatrix[:,axisY]
    pcaYCoord = pcaYCoord.real
    #fig2.set_title('Density plot of main PCA components')
    H, edgeX, edgeY = numpy.histogram2d(pcaXCoord, pcaYCoord, bins = NBINS)
    H = numpy.rot90(H)
    H = numpy.flipud(H)
    # mask zeroes
    maskedH = numpy.ma.masked_where(H==0, H)
    #Plot the histogram
    fig2 = plt.figure()
    plt.pcolormesh(edgeX, edgeY, maskedH)
    plt.xlabel('Pricinpal Component %d' % axis1)
    plt.ylabel('Principal Component %d' % axis2) 
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    fig2.savefig(outName, format='png')
    axisData = plt.axis()
    #show()
    return axisData

def obtainReducedY(dim,pcaObject):
    """obtain the PC of the original training data"""
    pcaMatrix = pcaObject.Y[:,:dim]
    return pcaMatrix

def projectMat(pca, matrix, dim):
    """perform the pca projection to the data matrix based on input pca object"""
    projected = numpy.zeros(shape=(len(matrix),dim))
    logging.info('Writing PC projections to data matrix')
    for i in xrange(len(matrix)):
        projected[i] = pca.project(matrix[i,:])[:dim]
    return projected
        
	    
def runPCA(args, trainMatrix, dim):
    """Wrapper function to run PCA, and return both the full projected PCs on the entire data"""
    pcaObject = doPCA(trainMatrix)
    #trainedMatrix = obtainReducedY(dim,pcaObject)  [FOR DEBUGGING ONLY]
    projected = projectMat(pcaObject, trainMatrix, dim)
    #xreduced, factors, evals, evecs = doPCA(kmerMatrix.T)
    return projected
    

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
    nClassify = 2
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
        print "confidence level have to be a value between 0 and 1"
        sys.exit(1)
    trim = math.ceil(size*(1-confidence)/2)
    ind = numpy.lexsort((matrix[:,1],matrix[:,0]))
    sortedMat = matrix[ind]
    sortedMat = sortedMat[0+trim-1:-trim]
    #keep track of the unsorted sequence id
    origIndex     = ind[0+trim-1:-trim]
    return sortedMat, origIndex


def getBinSlots(args, threshold):
    """Function to determine the required number of bins (per axis)"""
    seqSize = args.numberOfSeq
    nBins   = math.floor((seqSize/ threshold)**0.5)
    return nBins


def calcBins(matrix, origIndex, axisData, nbinX, nbinY, margin):
    """Function to convert the projected PCA matrix data into Bin counts"""
    binDict= {}
    #minX = math.floor(min(matrix[:,0])*margin)/margin
    #minY = math.floor(min(matrix[:,1])*margin)/margin
    #maxX = math.ceil(max(matrix[:,0])*margin)/margin
    #maxY = math.ceil(max(matrix[:,1])*margin)/margin
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
    writeMat('rawMatrix.txt', binMatrix) 
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
        #else:
            #binMatrix[slotX, slotY] = clusterDict.get(origIndex[i],0)
                
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
    """This function converts the corresponding region of a matrix to connected to a specific label"""
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
########MORPH FUNCTIONS##############
#####################################

def manualErode(matrix):
    """This function performs a 1-step erosion on the bin-matrix"""
    original = matrix
    matB = numpy.copy(matrix)
    erodeList = []
    #FOR DEBUGGING ONLY
    writeMat('originalMatrix.txt', matB)
    for i in xrange(len(original)):
        for j in xrange(len(original[i])):
            if original[i,j] == 0:
                continue    
            if hasZeroAround(original, i, j):
                matB[i,j] = 0
    writeMat('AfterErode.txt', matB)
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
    
def filterRegions(args, regions, origin):
    """Based on the number of species to be classified, this function returns the largest region to be used for classification if multiple 
    clusters exists, and returns the bin-matrix with only nClassify=2 regions, and a list of the corresponding region labels """
    matrix = numpy.copy(origin)
    nspecies = 2
    clusters = 0
    regionList  = range(0, regions+1)
    regionCount = [0,0]
    for i in regionList:
        if i<2:
            continue
        area = getRegionSize(matrix, i)
        regionCount.append(area)
    print "regionCounts [0,1,2,...]"
    print regionCount    
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
    
    
def erodeMap(args, matrix, threshold):
    """This is a wrapper function to perform all the erosion, filtering of regions, and classification 
    of clusters. The function require the raw Bin-matrix with original counts and a threshold parameter to be specified.
    The function returns the bin-matrix with only nClassify=2 regions, and a list of the corresponding region labels."""
    iterations = 1
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
    #TEST
    writeMat('01Matrix.txt', erodeMat)
    #perform binary erosion via iteration
    eroded = manualErode(erodeMat)
    test = getRegionSize(eroded, 1) 
    regions, eroded = labelRegions(eroded)
    filtered, nspecList = filterRegions(args, regions, eroded)
    writeMat('filtered.txt', filtered) 
    #LIST OF THE RELEVANT REGION TO BE USED TO CLASSIFY SPECIES ## FOR DEBUG ONLY 
    print "nspecList given as"
    print nspecList                                            
    return filtered, nspecList


####MAYBE NOT...
####TODO can try to implement multiple classification (n>2) ###

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
                #print "hCode found at i : %d, j : %d" %(i,j)
                slotId = slotTag(i,j)
                entry  = binDict.get(slotId,[])
                hostTrainList.extend(entry)
            elif filteredMat[i,j] == sCode:
                #print "sCode found at i : %d, j : %d" %(i,j)
                slotId = slotTag(i,j)
                entry  = binDict.get(slotId,[])
                symbTrainList.extend(entry) 
    return hostTrainList, symbTrainList


def writeTrainFile(args, hTrainList, sTrainList, spec1Path, spec2Path, randList):
    """This function writes the corresponding sequence from the training lists into fasta files"""
    spec1list = []
    spec2list = []
    limit = min(len(hTrainList),len(sTrainList))
    for i in hTrainList:
        origSeq = [randList[i]]
        spec1list.extend(origSeq)
    for i in sTrainList:
        origSeq = [randList[i]]
        spec2list.extend(origSeq)
    spec1list = dict.fromkeys(spec1list)
    spec2list = dict.fromkeys(spec2list)
    path = args.queries
    handle1 = open(spec1Path+'.fasta', 'w')
    handle2 = open(spec2Path+'.fasta', 'w')
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
    
    
def replot(args, dataMatrix, hTrainList, sTrainList):    
    xlist = []
    ylist = []
    i = 0
    for h in hTrainList:
        xp = dataMatrix[h,0]
        yp = dataMatrix[h,1]
        xlist.append(xp)
        ylist.append(yp)
    for s in sTrainList:
        xp = dataMatrix[s,0]
        yp = dataMatrix[s,1]
        xlist.append(xp)
        ylist.append(yp)
    mat2 = numpy.vstack([xlist, ylist])
    plotDensityPCA(args, mat2.T, 1, 2, 3)    
    

def clusterBins(args, matrix, confidence, percentile, clusterDict):
    """Wrapper function to obtain the eroded bin-matrix, and then the filtered bin-matrix."""  
    logging.info('Morphing bin-maps via erosion') 
    axisData = get_axis(matrix)
    sortedMatrix, origIndex = trimMatrix(matrix, confidence)
    #nBins = getBinSlots(args, threshold)
    binMatrix, binDict = calcBins(sortedMatrix, origIndex, axisData, NBINS, NBINS, 1000)
    binMatrix = posteriorErode(matrix, binMatrix, clusterDict, origIndex, axisData, NBINS, NBINS)
    threshold = thresholding(binMatrix, percentile)
    filteredMat, nspecList = erodeMap(args, binMatrix, threshold)   
    return filteredMat, binDict, nspecList


def autoTrainPCA(args, confidence, matrix, randList, percentile):
    """Wrapper function for the entire PCA auto-training process"""
    clusterDict = fitGMM(matrix, 50, 10, 0)
    finalMatrix, binDict, nspecList = clusterBins(args, matrix, confidence, percentile, clusterDict) 
    hTrainList, sTrainList = filterTrainSeq(finalMatrix, binDict, nspecList)
    spec1Path = args.spec1Path
    spec2Path = args.spec2Path
    hTrainList, sTrainList = writeTrainFile(args, hTrainList, sTrainList, spec1Path, spec2Path, randList)
    ####FIXME For DEBUG only####
    replot(args, matrix, hTrainList, sTrainList)
                
### WAS PREVIOUSLY USED WHEN TRY TO RUN THINGS IN R#####                
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
                        help='The query fasta sequences')
    #subparsers = parser.add_subparsers(help='Choose between option_1 or option_2 input format')
    #group = parser.add_mutually_exclusive_group()
    #parser_1 = subparsers.add_parser('option_1', help='Provide raw protein sequences, and perform blast before preparation for SVM')
    #parser_2 = subparsers.add_parser('option_2', help='Provide blast results as an input, directly start the preparation for SVM')
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
    trainMatrix, sampledList = computeKmerMat(args, full=True)
    
    logging.info('kmer Matrix created. preparing PCA') 
    projected = runPCA(args, trainMatrix, 2)
    plotDensityPCA(args, projected, 1, 2, 1)
    logging.info('performing PCA Auto-training') 
    autoTrainPCA(args, CONF_LEVEL, projected, sampledList, 80)
    
    logging.info("Auto-Training completed successfully.")

if __name__ == '__main__':
    main()
# vim:ts=4:sw=4:sts=4:et:ai:
