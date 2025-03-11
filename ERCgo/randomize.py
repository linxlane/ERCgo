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

def randomizeABColsUnmatched(colA, colB):
  rand_GeneA_HOG_col = []
  rand_GeneB_HOG_col = []
  
  num_rows = len(colA)

  for i in range(num_rows):
    #Pick a gene from col A and col B randomly
    randHogA = random.choice(colA)
    randHogB = random.choice(colB)
    
    #Make sure that randomly selected gene A and B aren't the same
    while randHogB == randHogA:
      randHogA = random.choice(colA)
      randHogB = random.choice(colB)
    
    #Add to column list
    rand_GeneA_HOG_col.append(randHogA)
    rand_GeneB_HOG_col.append(randHogB)

  return rand_GeneA_HOG_col, rand_GeneB_HOG_col

def generateRandomizedFiles(edgeFile, randomDirectory, reps, method):
    #Read edge file into dataframe and extract gene A and gene B columns as lists
    hogPairsDF = pandas.read_csv(edgeFile, sep='\t')
    hogAColList = hogPairsDF['GeneA_HOG'].tolist()
    hogBColList = hogPairsDF['GeneB_HOG'].tolist()

    #Create new empty dataframe for to store columns
    randomizedDF = pandas.DataFrame()

    ##Loop genereates X versions of ERCnet edge file where X is the value of -r
    #Keep A column the same and randomize genes in B column
    if method == 'B':
      print('Randomizing B column only')
      for i in range(1, reps + 1):
        print('> Generating random edge file: ' + str(i) + '...', flush=True, end='')
        randColBList = randomizeBColUnmatchedA(hogAColList, hogBColList)
        randomizedDF['GeneA_HOG'] = hogAColList
        randomizedDF['Rand_GeneB_HOG'] = randColBList
        newFilePath = randomDirectory + '/' + 'network_edges_' + 'rand_B_' + str(i) + '.tsv'
        randomizedDF.to_csv(newFilePath, sep='\t', index=False)
        print('Complete!')

    #Randomize both A and B columns
    if method == 'AB':
      print('Randomizing A and B columns')
      for i in range(1, reps + 1):
        print('> Generating random edge file: ' + str(i) + '...', flush=True, end='')
        randColAList, randColBList = randomizeABColsUnmatched(hogAColList, hogBColList)
        randomizedDF['Rand_GeneA_HOG'] = randColAList
        randomizedDF['Rand_GeneB_HOG'] = randColBList
        newFilePath = randomDirectory + '/' + 'network_edges_' + 'rand_AB_' + str(i) + '.tsv'
        randomizedDF.to_csv(newFilePath, sep='\t', index=False)
        print('Complete!')


def collectRandomizedEdgeFiles(randDirPath):
  #Collect randomized edge files for analysis
  randEdgeFilesList = glob.glob(randDirPath + '/*.tsv')
  randEdgeFilesList.sort()
  return randEdgeFilesList
