import cli
import hog_comp_ids
import goatools_GAF

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
