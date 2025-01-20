import cli
import hog_comp_ids
import pandas

#Parse user input from the command line
args = cli.runParser()
#Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

#print(argsDict)
edgeFilePath = hog_comp_ids.findEdgeFile(argsDict['input'])
verticesFilePath = hog_comp_ids.findVerticesFile(argsDict['input'])

hogCompDict = hog_comp_ids.generate_lookup_dict(verticesFilePath)

geneLookupDF_FULL = hog_comp_ids.generateHogCompTable(edgeFilePath, hogCompDict)

geneLookupDF_DROP = hog_comp_ids.dropNaRows(geneLookupDF_FULL)
