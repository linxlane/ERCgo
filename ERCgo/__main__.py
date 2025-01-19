import cli
import file_system
import prep

#Parse user input from the command line
args = cli.runParser()
#Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

#print(argsDict)
edgeFilePath = file_system.findEdgeFile(argsDict['input'])
verticesFilePath = file_system.findVerticesFile(argsDict['input'])

hogCompDict = prep.generate_lookup_dict(verticesFilePath)

with open(edgeFilePath, 'r') as edgeFile:
  #Skip first line with column titles
  next(edgeFile)
  for line in edgeFile:
    hogPair = line.strip().split('\t')
    compIDA = prep.lookup(hogPair[0], hogCompDict)
    compIDB = prep.lookup(hogPair[1], hogCompDict)
    print(compIDA + ' --> ' + compIDB)