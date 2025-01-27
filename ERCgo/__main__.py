import cli
import hog_comp_ids
import goatools_GAF
import shared_go
import randomize
import os

def collectFiles(path, endStr = '.tsv'):
    files = []
    for file in os.listdir(path):
        if file.endswith(endStr):
            files.append(file)
    files.sort()
    return files

def getName(file):
    return file[14:20]

def edgeFileSharedGO(edgeFilePath, verticesFilePath, masterOutPath, edgeFileName):

  intermediateFilesPath = masterOutPath + '/Intermediate_Files_' + edgeFileName
  os.makedirs(intermediateFilesPath)

  print('> Generate dictionary of {HOG_ID:Comprehensive_ID} to convert edge file HOG IDs to comprehensive IDs')
  hogCompDict = hog_comp_ids.generateLookupDict(verticesFilePath)

  print('> Generate table of [HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]')
  hogCompPath = intermediateFilesPath + '/HOG_COMP_TABLE_FULL_' + edgeFileName + '.tsv'
  geneLookupDF_FULL = hog_comp_ids.generateHogCompTable(edgeFilePath, hogCompDict, hogCompPath)

  print('> Drop rows that contain None values, ie there was no COMP_ID match for a given HOG_ID')
  hogCompNaPath = intermediateFilesPath + '/HOG_COMP_TABLE_DROP_NA_' + edgeFileName + '.tsv'
  compPairsNaPath = intermediateFilesPath + '/COMP_PAIRS_DROP_NA_' + edgeFileName + '.tsv'
  geneLookupDF_DROP = hog_comp_ids.dropNaRows(geneLookupDF_FULL, hogCompNaPath, compPairsNaPath)
  print('------------------------------------')

  print('> Utilize GOATOOLS package to extract GO terms for COMP_IDs')
  goAssocDF = goatools_GAF.generateGeneGoAssocDF(argsDict['gaf'])

  print('> Write [COMP_ID, GO_Terms] table to tsv: ID_GO_TERMS_TABLE.tsv')
  compGoPath = intermediateFilesPath + '/ID_GO_TERMS_TABLE_' + edgeFileName + '.tsv'
  goAssocDF.to_csv(compGoPath, sep='\t', index=False)
  print('------------------------------------')

  print('> Generate dictionary of {Gene_ID:GO_Terms} to find GO terms associated with genes in ERCnet identified gene pairs')
  assoc_dict = shared_go.generateAssocDict(goAssocDF)
  print('------------------------------------')

  print('> Collect GO terms associated with ERCnet pairs and generate table')
  compPairsWritePath = intermediateFilesPath + '/COMP_GO_TABLE_' + edgeFileName + '.tsv'
  genePairWithGoDF = shared_go.genePairGO(hogCompNaPath, assoc_dict, compPairsWritePath)
  print('------------------------------------')

  print('> Compare GO terms in gene pairs and determine intersection, if any')
  sharedGoWritePath = masterOutPath + '/SHARED_GO_TABLE_' + edgeFileName + '.tsv'
  shared_go.compareGoTerms(compPairsWritePath, genePairWithGoDF, sharedGoWritePath)


print('Starting ERCnet Gene Ontology Analysis!')
print('========================================================================')

#Parse user input from the command line
args = cli.runParser()
#Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

print('Searching for ERCnet output files...')
ercNetEdgeFilePath = hog_comp_ids.findEdgeFile(argsDict['input'])
ercNetVerticesFilePath = hog_comp_ids.findVerticesFile(argsDict['input'])
print('------------------------------------')

if not os.path.exists(argsDict['output']):
    print('A directory does not exist at this path.')
    print('Making directory for writing all output files...')
    os.makedirs(argsDict['output'])
    print('Successful')
else:
    print('Existing output directory found.')

print('========================================================================')

print('Generating randomized versions of ERCnet edge file for statistical analysis')
masterOut = argsDict['output']

###Generate and analyze randomized verions of edge file for stats analysis
randDirPath = randomize.randomDirectory(masterOut)
randomize.generateRandomizedFiles(ercNetEdgeFilePath, randDirPath, 3)

print('========================================================================')

randEdgeFilesList = collectFiles(randDirPath)

for edgeFile in randEdgeFilesList:
  print('Now processing: ' + getName(edgeFile))
  edgeFileSharedGO(randDirPath + '/' + edgeFile, ercNetVerticesFilePath, masterOut, getName(edgeFile))
  print('========================================================================')

edgeFileName = 'ERCnet_Network'
print('Now processing: ' + edgeFileName)
edgeFileSharedGO(ercNetEdgeFilePath, ercNetVerticesFilePath, masterOut, edgeFileName)
print('========================================================================')

print('Go analysis complete!')
