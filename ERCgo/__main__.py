import cli
import hog_comp_ids
import goatools_GAF
import shared_go
import randomize

print('Starting ERCnet Gene Ontology Analysis!')
print('====================================')

#Parse user input from the command line
args = cli.runParser()
#Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

print('Searching for ERCnet output files...')
edgeFilePath = hog_comp_ids.findEdgeFile(argsDict['input'])
verticesFilePath = hog_comp_ids.findVerticesFile(argsDict['input'])
print('------------------------------------')

###Generate and analyze randomized verions of edge file for stats analysis
randDirPath = randomize.randomDirectory(argsDict['output'])
randomize.generateRandomizedFiles(edgeFilePath, randDirPath, 3)

print('------------------------------------')

print('Generate dictionary of {HOG_ID:Comprehensive_ID} to convert edge file HOG IDs to comprehensive IDs')
hogCompDict = hog_comp_ids.generateLookupDict(verticesFilePath)

print('Generate table of [HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]')
geneLookupDF_FULL = hog_comp_ids.generateHogCompTable(edgeFilePath, hogCompDict, argsDict['output'])

print('Drop rows that contain None values, ie there was no COMP_ID match for a given HOG_ID')
geneLookupDF_DROP = hog_comp_ids.dropNaRows(geneLookupDF_FULL, argsDict['output'])
print('------------------------------------')

print('Utilize GOATOOLS package to extract GO terms for COMP_IDs')
goAssocDF = goatools_GAF.generateGeneGoAssocDF(argsDict['gaf'])
print('Write [COMP_ID, GO_Terms] table to tsv: ID_GO_TERMS_TABLE.tsv')
goAssocDF.to_csv(argsDict['output'] + '/ID_GO_TERMS_TABLE.tsv', sep='\t', index=False)
print('------------------------------------')

print('Generate dictionary of {Gene_ID:GO_Terms} to find GO terms associated with genes in ERCnet identified gene pairs')
assoc_dict = shared_go.generateAssocDict(goAssocDF)
print('------------------------------------')

print('Collect GO terms associated with ERCnet pairs and generate table')
genePairWithGoDF = shared_go.genePairGO(argsDict['output'] + '/COMP_PAIRS_DROP_NA.tsv', assoc_dict, argsDict['output'])
print('------------------------------------')

print('Compare GO terms in gene pairs and determine intersection, if any')
shared_go.compareGoTerms(argsDict['output'] + '/COMP_GO_TABLE.tsv', genePairWithGoDF, argsDict['output'])
