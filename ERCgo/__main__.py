import cli
import in_out
import gaf
import shared_go
import hog_comp_ids
import os

#######################
# Start of ERCgo main #
#######################
print('#######################################')
print('Starting ERCnet Gene Ontology Analysis!')
print('#######################################')

print('=======================================================================================================')
print('Conducting preliminary steps for ERCgo analysis: collecting, generating, and organizing files and data')
print('=======================================================================================================')

print('---------------------------------------------------------------------------------------------------')
print('Parsing cli arguments, fetching required files from input path, and resolve output file path')
print('---------------------------------------------------------------------------------------------------')

##Parse user input from the command line
args = cli.runParser()

##Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

##Fetch ERCnet files and GAF from parser input directory
print('Searching for ERCnet output files...')
#ERCnet ouput edge file
ercNetEdgeFilePath = in_out.findEdgeFile(argsDict['input'])
#ERCnet output vertices file
ercNetVerticesFilePath = in_out.findVerticesFile(argsDict['input'])
#GAF file
gafFilePath = in_out.findGafFile(argsDict['input'])

print('\n')

##Check for existing output directory, create a directory if it does not exist
print('Checking and resolving output cli argument if needed...')
masterOutPath = in_out.checkOutputDirectory(argsDict['output'])

print('\n\n')

###############
# Process GAF #
###############

print('---------------------------------------------------------------------------------------------------')
print('Processing GAF')
print('---------------------------------------------------------------------------------------------------')

##Create dictionary of genes and associated GO terms by reading and processing GAF
geneGoDict = gaf.processGaf(gafFilePath, masterOutPath)
print('\n')

#################################
# Create HOG_COMP ID dictionary #
#################################

print('---------------------------------------------------------------------------------------------------')
print('Create HOG-COMP ID dictionary')
print('---------------------------------------------------------------------------------------------------')

print('> Generating dictionary of {HOG_ID:Comprehensive_ID} to convert edge file HOG IDs to comprehensive IDs')
##Create dictionary of hog id equivalent comp ids for edge file conversion
hogCompDict = hog_comp_ids.generateHogCompDict(ercNetVerticesFilePath)

print('\n')

##################
# ERCgo Analysis #
##################

print('=======================================================================================================')
print('Prepartion complete. Beginning ERCgo analysis...')
print('=======================================================================================================')

##ERCnet data
print('---------------------------------------------------------------------------------------------------')
print('Starting analysis of ERCnet data')
print('---------------------------------------------------------------------------------------------------')

edgeFileName = 'ERCnet_Network'
print('Now processing: ' + edgeFileName)

#Get GO terms for each gene pair and create table
geneGoDF, writePath, frequencies = shared_go.generateSharedGOTable(masterOutPath, ercNetEdgeFilePath, hogCompDict, geneGoDict, edgeFileName)

#Analyze GO term sets returned for each gene pair, find intersection, and calculate score
shared_go.analyzeSharedGo(geneGoDF, masterOutPath, writePath, frequencies, edgeFileName)
print('\n')

print('#####################')
print('Go analysis complete! ')
print('#####################')
print('\n')
