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

if sys.version_info[0] < 3:
	from Queue import Queue
else:
	from queue import Queue


HOST_NAME  = "host"
SYMB_NAME  = "symb"
DB_NAME    = "HostSymbDB"
DB_FASTA   = DB_NAME + '.fasta'
BLAST_FILE = HOST_NAME + SYMB_NAME + '_test_OUTPUT.txt'
HOST_CODE  = 1
SYMB_CODE  = 2

HOST_TRAINING = HOST_NAME + '_training.fasta'
HOST_TESTING  = HOST_NAME + '_testing.fasta'
SYMB_TRAINING = SYMB_NAME + '_training.fasta'
SYMB_TESTING  = SYMB_NAME + '_testing.fasta'
LIBSVM_DIR    = 'libsvm'

LETTERS = ('A', 'T', 'G', 'C')

###TODO  Read about classes first#####
class PsyTransOptions:
""" this class consists of attributes to allow database and file paths to be obtained conveniently. 
	Attributes are referred from args, and written to write filenames, filepaths, naming conventions."""
	def __init__(self, args):
		self.args              = args 
		self.dbPath            = None
		self.fastaDbPath       = None
		self.blastResultsPath  = None
		self.blastInternalPath = None
		self.suffix            = None
		self.inputFile         = None
		self.trainFile         = None
		self.testFile          = None

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
		
	def isBlastGiven(self):
	"""check if any blast results are being provided by the user"""
		self.blastPath = self.args.blastResults
		if self.blastPath:
			return True
		else:
			return False	
		    
	def getBlastResultsPath(self): 
	"""to obtain the path to the user's Blast Result"""
		cwd = os.getcwd()
		if not self.blastResultsPath:
			self.blastResultsPath = os.path.join(cwd, str(self.args.blastResults))
		return self.blastResultsPath
	
	def getBlastInternalPath(self):
	"""to obtain the path to the Blast Result internally obtained by the script"""
		if not self.blastInternalPath:
			self.blastInternalPath = os.path.join(self.args.tempDir, BLAST_FILE)
		return self.blastInternalPath
		
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
		
	def getLength(self):
	""" obtain the length description of the kmer file name"""
		if self.args.numberOfSeq == 0:
			length = 'all'
		else:
			length = self.args.numberOfSeq
		return str(length)
		
	def getSuffix(self):
	""" writing of a naming convention."""
		if not self.suffix:
			suffix = self.getLength()
			self.mink = str(self.args.minWordSize)
			self.maxk = str(self.args.maxWordSize)
			self.suffix = suffix + '_c' + self.mink + '_k' + self.maxk
		return self.suffix
		
	def getTrainFile(self):
	"""writing the filename and type of the kmer training file as an input for the SVM process"""
		if not self.trainFile:
			fName = self.getSuffix()
			self.trainFile = 'Training' + '_' + fName + '.txt'
		return str(self.trainFile)
		
	def getTestFile(self):
	"""writing the filename and type of the kmer testing file as an input for the SVM process"""
		if not self.testFile:
			fName = self.getSuffix()
			self.testFile = 'Testing' + '_' + fName + '.txt'
		return str(self.testFile)
		
		
#######################################



def iterFasta(path):
	"""Iterates over the sequences of a fasta file"""
	logging.debug("Loading fasta files from %s" % path)
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
 
#Functions for 02    
def writeDatabase(args, options, fastaPath):
	"""writes multiple sources to a fasta file. This function also calls iterFasta() in the process."""
	logging.debug('Creating Database.')
	hostPath      = args.hostSeq
	symbPath      = args.symbSeq
	targetpath = open(fastaPath, 'w')
	#Writing Host Sequences to target database
	i = 0
	for name, seq in iterFasta(hostPath):
		i += 1
		name = '>' + HOST_NAME +'_' + str(i)
		targetpath.write(name + "\n")
		targetpath.write(seq + "\n")
	j = 0
	for name, seq in iterFasta(symbPath):
		j += 1
		name = '>' + SYMB_NAME +'_' + str(j)
		targetpath.write(name + "\n")
		targetpath.write(seq + "\n")
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
		logging.warning('Error detected. Please check logfile. Blast Database not created.')
		sys.exit()
	options.createCheckPoint('makeDB.done')
	    

#Function for 03 BLAST SEARCH
def runBlast(args, options):
	"""calls the blastx process from the command line via subprocess module. Result obtained
	from the BLAST search will be compressed into a zipped file. The output format of the result by default 
	is set to '7' (tab-seaparated with comments)."""
    #Define BLAST variables
	logging.info('Performing Blast Search')
	eVal         = '%.2e' % args.maxBestEvalue
	dbPath       = options.getDbPath()
	blastOutput  = options.getBlastInternalPath()
	blastCmd = ['blastx',
	            '-evalue',
	            eVal,
	            '-query',
	            args.queries,
                '-db', 
                dbPath,
	            '-outfmt 7',
	            '-out',
	            blastOutput]
	blastCmd = ' '.join(blastCmd)
	retCode  = subprocess.call(blastCmd, shell=True)
	if not retCode == 0:
		logging.warning('Error detected. Please check blastlogfile. Blast Search not executed or exit with error.')
		sys.exit(1)
	options.checkPoint('runBlast')
	logging.debug('Blast finished')


#Functions for 04 ParseBlast up to Preparing sets   
def parseBlast(args, options):
	"""parses the blast result given or previously obtained, to later be used to prepare training and testing set of 
	unambiguous sequences. This function returns an object called [querries] which summarises the blast results."""
	if not args.blastResults:
		path = options.getBlastInternalPath()
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
		if not qName in querries:
			querries[qName] = []
			n += 1
		hit = (hName, evalue)
		querries[qName].append(hit)
	logging.info('Parsed %d blast records' % n)
	logging.info('Found %d queries hits' % len(querries)) 
	handle.close()
	return querries

def classify(querries, args):
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

	classification = {}
	for qName in querries:
		hits = querries[qName]
		hits.sort(sortHits)
		hasCoral        = False
		hasZoox         = False
		coralBestEvalue = -1
		zooxBestEvalue  = -1
		for hName, evalue in hits:
			if hName.startswith('coral'):
				if not hasCoral:
					hasCoral        = True
					coralBestEvalue = evalue
			else :
				if not hasZoox:
					hasZoox        = True
					zooxBestEvalue = evalue
        # WARNING Could be improved
		if hasCoral and not hasZoox and coralBestEvalue <= args.maxBestEvalue:
			classification[qName] = HOST_CODE
		elif hasZoox and not hasCoral and zooxBestEvalue <= args.maxBestEvalue:
			classification[qName] = SYMB_CODE
	return classification

def seqSplit(args, options, classification):
	"""writes the relevant unambiguous sequences into 4 fasta files: training.fasta for host sequences,
	testing.fasta for host sequences, training.fasta for symb sequences and testing.fasta for symb sequences.
	The size/ratio of the training to testing set can be modified by the user from the script option 
	(--trainingTestingRatio)."""
	m = 0
	j = 0
	d = args.trainingTestingRatio 
	handle = open(args.queries)
	hostTrain = open(os.path.join(args.tempDir, HOST_TRAINING), 'w')
	hostTest  = open(os.path.join(args.tempDir, HOST_TESTING), 'w')
	symbTrain = open(os.path.join(args.tempDir, SYMB_TRAINING), 'w')
	symbTest  = open(os.path.join(args.tempDir, SYMB_TESTING), 'w')
	for name, seq in iterFasta(args.queries):
		identity = identity = (name.split(' ')[0])[1:]
		seqClass = classification.get(identity, 0)
		if seqClass == HOST_CODE:
			if m % d == 0:
				hostTrain.write(">" + identity + "\n")
				hostTrain.write(seq + "\n")
			else :
				hostTest.write(">" + identity + "\n")
				hostTest.write(seq + "\n")
			m += 1
		elif seqClass == SYMB_CODE:
			if j % d == 0:
				symbTrain.write(">" + identity + "\n")
				symbTrain.write(seq + "\n")
			else :
				symbTest.write(">" + identity + "\n")
				symbTest.write(seq + "\n")
			j += 1

	handle.close()
	hostTest.close()
	hostTrain.close()
	symbTest.close()
	symbTrain.close()
	logging.info('Found %d unambiguous hits' % len(classification))
	logging.info('Found %d host only hits' % m)
	logging.info('Found %d symbion only hits' % j)
	options.createCheckPoint('parseBlast.done')  


#Functions for 05 Kmer Computations
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


def computerKmers(path, outfile, code, args, mode):
	"""computes the kmer counts throughout the kmer range for each sequence, and write the output to a .txt file.
	Each kmer counts will be scaled accordingly with the sequence size. This function calls prepareMaps() in the process."""
	label = int(code)
    # Prepare all maps
	kMin    = args.minWordSize
	kMax    = args.maxWordSize
	maps    = []
	for i in xrange(kMin, kMax + 1):
		maps.append(prepareMaps(0, i, ['']))
    # Initialise output
	out     = str(outfile)
	outPath = os.path.join(args.tempDir,out)
	handle  = open(outPath, mode)
    # Initialise counts
	counts  = {}
	for i in xrange(kMin, kMax + 1):
		counts[i] = numpy.zeros(4 ** i, numpy.double)
    # Iterate over sequences
	nSeqs   = 0
	for name, seq in iterFasta(path):
		if args.numberOfSeq > 0 and nSeqs >= args.numberOfSeq:
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

def iterKmers(args, options, kmerTrain, kmerTest):
	"""computes the kmer count map for each of the training and testing sequences. The function outputs two files:
	a training file and a testing file to be used as inputs for the SVM training. This function calls computerKmers() in the process"""
	hostTrainpath = os.path.join(args.tempDir, HOST_TRAINING)
	hostTestpath  = os.path.join(args.tempDir, HOST_TESTING)
	symbTrainpath = os.path.join(args.tempDir, SYMB_TRAINING)
	symbTestpath  = os.path.join(args.tempDir, SYMB_TESTING)	
	mink = str(args.minWordSize)
	maxk = str(args.maxWordSize)
	logging.debug('Obtaining kmer counts')
	logging.debug('Creating kmer maps')
	computerKmers(hostTrainpath, kmerTrain, HOST_CODE, args, "w")
	computerKmers(hostTestpath, kmerTest, HOST_CODE, args, "w")
	computerKmers(symbTrainpath, kmerTrain, SYMB_CODE, args, "a")
	computerKmers(symbTestpath, kmerTest, SYMB_CODE, args, "a")
	options.createCheckPoint('kmers.done')

#################################################################
################################################################
#Calling chlin\'s SVM script'###################################
#################################################################
#################################################################

def callSVM(args, options, kmerTrain, kmerTest):
	"""calls the easy.py script through the command line using the subprocess module. 
	That script in turn calls grid.py which performs the training for the SVM."""
	easyPath  = os.path.join("libsvm", "easy.py")
	trainPath = os.path.join(args.tempDir, kmerTrain)
	testPath  = os.path.join(args.tempDir, kmerTest)
	easyCmd = ['python ', easyPath, ' ', trainPath, ' ', testPath] 
	easyCmd = ''.join(easyCmd)
	logging.debug('Calling easy.py python script')
	subprocess.Popen(easyCmd, shell=True)
	options.createCheckPoint('svm.done')


#Prediction SVM
def loadPredictions(path):
	handle      = open(path)
	content     = handle.read()
	predictions = content.strip().split('\n')
	handle.close()
	return predictions


def writeOutput(predictions, fastaPath, prefix1, prefix2):
	size = len(predictions)
	out1 = open(prefix1 + fastaPath, "w")
	out2 = open(prefix2 + fastaPath, "w")
	j = 0
	if j < xrange(size):
		assert predictions[j] in '12'
		for name, seq in iterFasta(fastaPath):
			if predictions[i] == '1':
				out1.write(name + "\n")
				out1.write(seq + "\n")
				j += 1
			else:
				out2.write(name + "\n")
				out2.write(seq + "\n")
				j += 1
	out1.close()
	out2.close()


def predictSVM(args, kmerTrain, kmerTest):
	scaledFile      = kmerTrain + '.scale'
	modelFile       = kmerTrain + '.model'
	rangeFile       = kmerTrain + '.range' 
	scaledTestFile  = kmerTest  + '.scale'
	predictTestFile = kmerTest  + '.predict'
	logging.info('Predicting and Classifying Sequences')
	fastaPath = args.queries
	kmerScale = '.'.join(fastaPath.split('.')[:-1]) + '.scaled'
	kmerScale = os.path.join(args.tempDir, kmerScale)
	kmerPred  = '.'.join(fastaPath.split('.')[:-1]) + '.pred'
	kmerPred = os.path.join(args.tempDir, kmerPred)
	kmerFile  = '.'.join(fastaPath.split('.')[:-1]) + '.kmers'
	kmerOut = os.path.join(args.tempDir, kmerFile)
	###FIXME Need to work on writing the output of SVM scripts to Temp, Currently every output from this script
	## goes to Temp, but the grid.py & easy.py writes scale data files to current DIR. CHECK LATER ##### 
	#This function by default writes output file to temp DIR
	computerKmers(args.queries, kmerFile, HOST_CODE, args, "w")
	#SVM_Scale
	scaleCmd   = [ svmscale_exe, ' -r ', rangeFile, ' ', kmerOut,' > ', kmerScale]
	scaleCmd   = ''.join(scaleCmd)
	subprocess.Popen(scaleCmd, shell=True)
	#SVM_predict
	predictCmd = [ svmpredict_exe, ' -r ', kmerScale, ' ', modelFile, ' ', kmerPred]
	predictCmd = ''.join(predictCmd)
	subprocess.Popen(predictCmd, shell=True)
	#parse_Prediction
	predictions = loadPredictions(kmerPred)
	writeOutput(predictions, fastaPath, HOST_NAME, SYMB_NAME)
	logging.info('Classification of Sequences is now completed.')
	### TODO Implement Threading?


def tempPathCheck(args):
	"""to check if the temporary folder exists, else it will be created."""
	logging.info('Checking for temporary file folder')
	tempFolder = os.path.abspath(args.tempDir) 
	if not os.path.isdir(tempFolder):
		os.makedirs(tempFolder)

def checkExecutable(program):  
	"""to check whether the required program(.exe) has been installed and can be executable."""
	path = os.getenv('PATH')
	root = os.path.dirname(os.path.abspath(sys.argv[0]))
	libPath = os.path.join(os.getcwd(), LIBSVM_DIR)
	for d in path.split(os.path.pathsep):  
		exe = os.path.join(d, program) 
		if os.path.exists(exe) and os.access(exe, os.X_OK):
			return exe
	for dir_, _, files in os.walk(libPath):
		dir_ = os.path.join(dir_, program)
		if os.path.exists(dir_) and os.access(line, os.X_OK):
			return dir_


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
                        #nargs='+',
                        help='Blast results obtained')
	parser.add_argument('-e',
                        '--maxBestEvalue',
                        type=float,
                        default='1e-20',
                        help='Maximum value for the best e-value')
	parser.add_argument('--trainingTestingRatio',
                        type=float,
                        default='2',
                        help='Value used to divide classfied sequences into testing and training sets')
    ### FIXME check how this plays with the above option
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
                        type = bool,
                        default=True,
                        help='Turn Verbose mode on?')
	parser.add_argument('-t',
                        '--tempDir',
                        type=str,
                        default='temp',
                        help='Name of temporary directory')
	parser.add_argument('-X',
                        '--clearTemp',
                        type=bool,
                        default=False,
                        help='Clear all temporary data upon completion?')
	parser.add_argument('-z',
                        '--stopAfter',
                        type=str,
                        choices=('db','runBlast','parseBlast','kmers','SVM'),
                        help='Optional exit upon completion of stage.')
	parser.add_argument('-R',
                        '--restart',
                        type=bool,
                        default=False,
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
		logging.error('Either provide the host and symbiont sequences OR the output(s) of the blast results')
		sys.exit()
	if args.verboseMode:
		logging.getLogger().setLevel(logging.DEBUG)
	logging.info("Arguments parsed. Starting...")
	tempPathCheck(args)

	blastPath = options.getBlastResultsPath() ### FIXME dont double guess the user
	if not (options.isBlastGiven() and os.path.exists(blastPath)):
		if not options.isBlastGiven():
    	#Step 1
			fastaPath  = options.getFastaDbPath()
			dbPath     = options.getDbPath()
			if not (restart and options.checkPoint("writeDatabase.done")):
				writeDatabase(args, options, fastaPath)
				if checkExecutable('makeblastdb'):
					makeDB(args, options)
				else:
					logging.warning('makeblastdb not found. Exiting')
					sys.exit()
			if args.stopAfter == 'db':
				logging.info('Stop after "db" requested, exiting now')
				sys.exit()

    		#Step 2
			if not (restart and options.checkPoint("runBlast.done")):
				if checkExecutable('blastx'):
					runBlast(args, options)
				else:
					logging.warning('blastx not detected in PATH. Please check that the program is being installed correctly. Exiting')
					sys.exit()
			if args.stopAfter == 'runBlast':
				logging.info('Stop after "runBlast" requested, exiting now')
				sys.exit()
				
		else:
			logging.warning('Path to blast result does not exist. Exiting.')
			sys.exit()

    #Step 3
	if not (restart and options.checkPoint("parseBlast.done")):
		querries       = parseBlast(args, options)
		classification = classify(querries, args)
		seqSplit(args, options, classification)
	if args.stopAfter == 'parseBlast':
		logging.info('Stop after "parseBlast" requested, exiting now')
		return

    #Step 4
	#Kmer preparation
	kmerTrain = options.getTrainFile()
	kmerTest  = options.getTestFile()
	if not (restart and options.checkPoint("kmers.done")):
		iterKmers(args, options, kmerTrain, kmerTest)
	if args.stopAfter == 'kmers':
		logging.info('Stop after "kmers" requested, exiting now')
		return

    #Step 5
	if not (restart and options.checkPoint("svm.done")):
		if checkExecutable('svm-train') and checkExecutable('svm-scale') and checkExecutable('svm-predict'):
			callSVM(args, options, kmerTrain, kmerTest)
		else:
			logging.warning('SVM package not complete. please check that svm-train, svm-scale and svm-predict has been correctly installed')
			sys.exit()
	if args.stopAfter == 'SVM':
		logging.info('Stop after "SVM" requested, exiting now')
		return

    #Step 6
	predictSVM(args, kmerTrain, kmerTest)
	logging.info("SVM classification completed successfully.")
	
	if args.clearTemp:
		shutil.rmtree(args.tempDir)

if __name__ == '__main__':  
	main()
