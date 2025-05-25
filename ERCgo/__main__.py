import cli
import in_out
import gaf
import shared_go
import os
import population

def analysisPipeline(genePairsDF, genePairsFilePath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict):
  goTermsFreq = population.calculatePopulationFrequencies(genePairsDF, geneGoDict, intermediateFilesPath, argsDict['job_name'])
  genePairsGoDF, genePairsGoPath = shared_go.collectGoTerms(genePairsFilePath, geneGoDict, intermediateFilesPath, argsDict)
  shared_go.analyzeSharedGo(genePairsGoDF, masterOutPath, genePairsGoPath, goTermsFreq, argsDict)

#######################
# Start of ERCgo main #
#######################
print('#######################################')
print('Starting ERCnet Gene Ontology Analysis!')
print('#######################################')

##Parse user input from the command line
args = cli.runParser()

##Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

print('=======================================================================================================')
print('Conducting preliminary steps for ERCgo analysis: collecting files and checking output directory')
print('=======================================================================================================')
print('---------------------------------------------------------------------------------------------------')
print('INPUT')
print('---------------------------------------------------------------------------------------------------')

##Check input
print('> Verify correct input files are present for specified analysis...')
gafFilePath = in_out.findGafFile(argsDict['input'])
interactomeFilePath = in_out.findInteractomeFile(argsDict['input'])
print('> DONE')

print('---------------------------------------------------------------------------------------------------')
print('OUTPUT')
print('---------------------------------------------------------------------------------------------------')
##Create new directory with the job_name where ERCgo output will be written
#Check for user defined output, otherwise get ERCgo path
#Check for existing directory at the output/job_name path, delete if it exists to start fresh
print('> Check cli output argument and resolve if needed...')
if argsDict['output'] is None:
  cwd = os.getcwd()
  outputDir = cwd + '/OUTPUT'
else:
  outputDir = argsDict['output']

masterOutPath = in_out.checkOutputDirectory(outputDir + '/' + argsDict['job_name'] + '_OUT')

#Make directory to write conversion files for later reference
intermediateFilesPath = masterOutPath + '/Intermediate_Files_' + argsDict['job_name']
os.makedirs(intermediateFilesPath)
print('> DONE')

###############
# Process GAF #
###############

print('---------------------------------------------------------------------------------------------------')
print('Processing GAF')
print('---------------------------------------------------------------------------------------------------')

##Create dictionary of genes and associated GO terms by reading and processing GAF
#GAF file
geneGoDict = gaf.processGaf(gafFilePath, intermediateFilesPath)
print('> DONE')

########################################
# Preprocess input data #
########################################
print('=======================================================================================================')
print('Preprocess input data for analysis')
print('=======================================================================================================')

genePairsDF, genePairsPath = in_out.formatInteractomeData(argsDict, intermediateFilesPath, interactomeFilePath)

########################################
# Analysis #
########################################
print('=======================================================================================================')
print('GO term analysis')
print('=======================================================================================================')

analysisPipeline(genePairsDF, genePairsPath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict)

print('#####################')
print('Go analysis complete!')
print('#####################')
print('\n')
