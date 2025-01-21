import cli
import hog_comp_ids
import goatools_GAF
import shared_go

#Parse user input from the command line
args = cli.runParser()
#Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

edgeFilePath = hog_comp_ids.findEdgeFile(argsDict['input'])
verticesFilePath = hog_comp_ids.findVerticesFile(argsDict['input'])

hogCompDict = hog_comp_ids.generateLookupDict(verticesFilePath)

geneLookupDF_FULL = hog_comp_ids.generateHogCompTable(edgeFilePath, hogCompDict, argsDict['output'])

geneLookupDF_DROP = hog_comp_ids.dropNaRows(geneLookupDF_FULL, argsDict['output'])

goAssocDF = goatools_GAF.generateGeneGoAssocDF(argsDict['gaf'])

goAssocDF.to_csv(argsDict['output'] + '/ID_GO_TERMS_TABLE.tsv', sep='\t', index=False)

assoc_dict = shared_go.generateAssocDict(goAssocDF)

shared_go.genePairGO(argsDict['output'] + '/COMP_PAIRS_DROP_NA.tsv', assoc_dict, argsDict['output'])