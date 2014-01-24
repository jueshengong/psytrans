#!/usr/bin/python

import argparse
import gzip
import logging
import sys

HOST_NAME = "Coral"
SYMB_NAME = "Zoox"
HOST_CODE = 0
SYMB_CODE = 1


def renameSeq(args):
    logging.info("Relabelling sequences")
    hostName = HOST_NAME
    symbName = SYMB_NAME
    hostPath = args.hostSeq
    symbPath = args.symbSeq
    return hostName, symbName, hostPath, symbPath


def iterFasta(path):
    name, seq = None, []
    logging.debug("Loading fasta files from %s" % path)
    if path.endswith('.gz'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        i = 1
	    line = line.strip()
            if not line:
                continue
        if line[0] == ">":
            if name:
                yield (name, seq)
            name = line.append('>' + species + '_' + str(i))
            seq = ""
            i += 1
        else:
            seq += line[:-1]
         
        # now a sequence is loaded with name and seq
        yield (name, seq)
    handle.close()
 
    
def writeDatabase(host, symb):    
    logging.debug("Constructing database")
    targetpath = open('target_db.fasta', 'w')
    #Writing Host Sequences to target database
    hPath = host
    i=0
    for name, seq in iterFasta(hPath):
        i += 1
        name = '>'+ hostName +'_'+str(i)
        targetpath.write(name + "\n")
        targetpath.write(seq + "\n") 
    #Writing Symbiont Sequences to target database
    sPath = symb
    j=0
    for name, seq in iterFasta(sPath):
        j += 1
        name = '>'+ symbName +'_'+str(j)
        targetpath.write(name + "\n")
        targetpath.write(seq + "\n") 
    targetpath.close()   
    logging.debug('Database has been created. Preparing for Blast..')
    
    
def runBlast():
    logging.info('Performing Blast Search')
    logging.debug('Blast started')
    logging.debug('Blast finished')
    
def parseBlastAndPrepareSets(args):
    if args.blastResults != '-':
       logging.debug('Blast Results obtained. Converting..')
    logging.info('Parsing Blast Results')
    logging.info('Prepare Training Sets')

def prepareInput():
    logging.info('Computing kmer tables')
    logging.debug('Obtaining kmer counts')
    logging.debug('Creating kmer maps')

#chlin\'s SVM script runs here'
def svmRun():
    logging.info('Performing SVM')

def predictClass():
    logging.info('Predicting and Classifying Sequences')

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
                        help= 'The input host sequences (single species)')
    parser.add_argument('-S',
                        '--symbSeq',
                        type=str,
                        help = 'The input symbiont sequences (single species)')
    
    parser.add_argument('-b',
                        '--blastResults',
                        type=str,
                        nargs='+',
                        help='Blast results obtained')
    parser.add_argument('-e',
                        '--maxBestEvalue',
                        type=float,
                        default='1e-20',
                        help='Maximum value for the best e-value')
    parser.add_argument('-d',
                        '--divisionOfSets',
                        type=float,
                        default='2',
                        help='Value used to divide classfied sequences into testing and training sets')
    parser.add_argument('-n',
                        '--numberOfSeq',
                        type=int,
                        default = '1000',
                        help='Maximum number of training & testing sequences')
    parser.add_argument('-c',
                        '--minWordSize',
                        type = int,
                        default = '1',
                        help = 'Minimum value of DNA word length')
    parser.add_argument('-k',
                        '--maxWordSize',
                        type = float,
                        default = '4',
                        help = 'Maxmimum value of DNA word length')
    parser.add_argument('-v',
                        '--verboseMode',
                        type = bool,
                        default = True,
                        help = 'Turn Verbose mode on?')
    parser.add_argument('-X',
                        '--clearTemp',
                        type = bool,
                        default = False,
                        help = 'Clear all temporary data upon completion?')
    args = parser.parse_args()
    return args

def main():
    logging.basicConfig(level=logging.INFO, format=("%(asctime)s - %(funcName)s - %(message)s"))
    args = mainArgs()
    
    if not (args.hostSeq and args.symbSeq and not args.blastResults) and \
       not (args.blastResults and not (args.hostSeq or args.symbSeq)):
       logging.error('Either provide the host and symbion sequences OR the output(s) of the blast results')
       sys.exit()
    if args.verboseMode:
       logging.basicConfig(level=logging.DEBUG, format=("%(asctime)s - %(message)s"))
    logging.info("Arguments parsed. Starting...")

    #Step 1
    renameSeq(args)
    #Step 2
    writeDatabase(hostPath, symbPath)
    #Step 3
    runBlast()
    #Step 4
    parseBlastAndPrepareSets(args)
    #Step 5
    prepareInput()
    #Step 6
    svmRun()
    #Step 7
    predictClass()
    logging.warning("SVM classification completed successfully.")

if __name__ == '__main__':
    main()
