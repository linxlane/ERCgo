import os
import pandas
import random
import glob

def randomDirectory(outDirectory):
  randomDirPath = outDirectory + '/randomEdgeFiles'
      
  if not os.path.exists(randomDirPath):
    os.makedirs(randomDirPath)
    print('Successful', end='')
  else:
    print('Directory for randomized edge files found. That\'s weird...')
  
  return randomDirPath

def randomizeBColUnmatchedA(colA, colB):
  rand_GeneB_HOG_col = []
  
  for hogA in colA:
    randHogB = random.choice(colB)
    while randHogB == hogA:
      randHogB = random.choice(colB)
    rand_GeneB_HOG_col.append(randHogB)

  return rand_GeneB_HOG_col

def generateRandomizedFiles(edgeFile, randomDirectory, reps):
    #Read edge file into dataframe and extract gene A and gene B columns as lists
    hogPairsDF = pandas.read_csv(edgeFile, sep='\t')
    hogAColList = hogPairsDF['GeneA_HOG'].tolist()
    hogBColList = hogPairsDF['GeneB_HOG'].tolist()

    #Create new empty dataframe for to store same A column and new randomized B column
    hogARandHogB_DF = pandas.DataFrame(data=hogAColList)

    #Loop genereates X versions of ERCnet edge file where X is the value of -r
    for i in range(1, reps + 1):
      print('> Generating random edge file: ' + str(i) + '...', flush=True, end='')
      randColBList = randomizeBColUnmatchedA(hogAColList, hogBColList)
      hogARandHogB_DF['Rand_GeneB_HOG'] = randColBList
      newFilePath = randomDirectory + '/' + 'network_edges_' + 'rand_' + str(i) + '.tsv'
      hogARandHogB_DF.to_csv(newFilePath, sep='\t', index=False)
      print('Complete!')

def collectRandomizedEdgeFiles(randDirPath):
  #Collect randomized edge files for analysis
  randEdgeFilesList = glob.glob(randDirPath + '/*.tsv')
  randEdgeFilesList.sort()
  return randEdgeFilesList
