import cli
import in_out
import gaf
import shared_go
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

print('\n')

##Create new directory with the job_name where ERCgo output will be written
#Check for user defined output, otherwise get ERCgo path
#Check for existing directory at the output/job_name path, delete if it exists to start fresh
print('Checking and resolving output cli argument if needed...')
if argsDict['output'] is None:
  outputDir = os.getcwd()
else:
  outputDir = argsDict['output']
masterOutPath = in_out.checkOutputDirectory(outputDir + '/' + argsDict['job_name'])

#Make directory to write conversion files for later reference
intermediateFilesPath = masterOutPath + '/Intermediate_Files_' + argsDict['job_name']
os.makedirs(intermediateFilesPath)

print('\n\n')

###############
# Process GAF #
###############

print('---------------------------------------------------------------------------------------------------')
print('Processing GAF')
print('---------------------------------------------------------------------------------------------------')

##Create dictionary of genes and associated GO terms by reading and processing GAF
#GAF file
gafFilePath = in_out.findGafFile(argsDict['input'])
geneGoDict = gaf.processGaf(gafFilePath, masterOutPath)
print('\n')

########################################
# Prep ERCnet hits output for analysis #
########################################

if argsDict['analysis'] == 'hits':
  in_out.formatErcNetDataHits(argsDict, intermediateFilesPath)

elif argsDict['analysis'] == 'full':
  print('Full')
elif argsDict['analysis'] == 'both':
  print('Both')

print('#####################')
print('Go analysis complete! ')
print('#####################')
print('\n')
