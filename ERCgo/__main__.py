import cli
import file_system

#Parse user input from the command line
args = cli.runParser()
#Generate dictionary of flags parsed in cli.py
argsDict = vars(args)

print(argsDict)
edgeFilePath = file_system.findEdgeFile(argsDict['input'])
verticesFilePath = file_system.findVerticesFile(argsDict['input'])
print(edgeFilePath)


