import cli
import in_out
import gaf
import shared_go
import os
import population

def analysisPipeline(genePairsDF, genePairsFilePath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict):
  goTermsFreq = population.calculatePopulationFrequencies(genePairsDF, geneGoDict, masterOutPath, argsDict['job_name'])
  genePairsGoDF, genePairsGoPath = shared_go.collectGoTerms(genePairsFilePath, geneGoDict, intermediateFilesPath, argsDict['job_name'])
  shared_go.analyzeSharedGo(genePairsGoDF, masterOutPath, genePairsGoPath, goTermsFreq, argsDict['job_name'])

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
print('Conducting preliminary steps for ERCgo analysis: collecting, generating, and formatting files and data')
print('=======================================================================================================')
print('---------------------------------------------------------------------------------------------------')
print('INPUT')
print('---------------------------------------------------------------------------------------------------')

##Check input
print('> Verify correct input files are present for specified analysis...')
gafFilePath = in_out.findGafFile(argsDict['input'])
if argsDict['analysis'] == 'hits':
  edgeFilePath = in_out.findEdgeFile(argsDict['input'])
  verticesFilePath = in_out.findVerticesFile(argsDict['input'])

elif argsDict['analysis'] == 'full':
  ercResultsFilePath = in_out.findErcResultsFile(argsDict['input'])

elif argsDict['analysis'] == 'both':
  edgeFilePath = in_out.findEdgeFile(argsDict['input'])
  verticesFilePath = in_out.findVerticesFile(argsDict['input'])
  ercResultsFilePath = in_out.findErcResultsFile(argsDict['input'])
print('> DONE')

print('---------------------------------------------------------------------------------------------------')
print('OUTPUT')
print('---------------------------------------------------------------------------------------------------')
##Create new directory with the job_name where ERCgo output will be written
#Check for user defined output, otherwise get ERCgo path
#Check for existing directory at the output/job_name path, delete if it exists to start fresh
print('> Check cli output argument and resolve if needed...')
if argsDict['output'] is None:
  outputDir = os.getcwd()
else:
  outputDir = argsDict['output']

masterOutPath = in_out.checkOutputDirectory(outputDir + '/' + argsDict['job_name'])

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
geneGoDict = gaf.processGaf(gafFilePath, masterOutPath)
print('> DONE')

########################################
# Format input data #
########################################
print('=======================================================================================================')
print('Format input data into [COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval] table for analysis')
print('=======================================================================================================')

if argsDict['analysis'] == 'hits':
  genePairsDF, genePairsDropNaPath = in_out.formatErcNetDataHits(argsDict, intermediateFilesPath, edgeFilePath, verticesFilePath)

elif argsDict['analysis'] == 'full':
  genePairsDF, genePairsDropNaPath = in_out.formatFullResults(argsDict, intermediateFilesPath, ercResultsFilePath)

elif argsDict['analysis'] == 'both':
  hitGenePairsDF, hitGenePairsDropNaPath = in_out.formatErcNetDataHits(argsDict, intermediateFilesPath, edgeFilePath, verticesFilePath)
  fullGenePairsDF, fullGenePairsDropNaPath = in_out.formatFullResults(argsDict, intermediateFilesPath, ercResultsFilePath)

########################################
# Analysis #
########################################
print('=======================================================================================================')
print('GO term analysis')
print('=======================================================================================================')

if argsDict['analysis'] == 'hits':
  analysisPipeline(genePairsDF, genePairsDropNaPath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict)

elif argsDict['analysis'] == 'full':
  analysisPipeline(genePairsDF, genePairsDropNaPath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict)

elif argsDict['analysis'] == 'both':
  analysisPipeline(hitGenePairsDF, hitGenePairsDropNaPath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict)
  analysisPipeline(fullGenePairsDF, fullGenePairsDropNaPath, geneGoDict, masterOutPath, intermediateFilesPath, argsDict)

print('#####################')
print('Go analysis complete! ')
print('#####################')
print('\n')
