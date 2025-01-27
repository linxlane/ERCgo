import os
import pandas

def randomDirectory(outDirectory):
  randomDirPath = outDirectory + '/randomEdgeFiles'
      
  if not os.path.exists(randomDirPath):
    print('Making directory for writing randomized edge files...')
    os.makedirs(randomDirPath)
    print('Successful')
  else:
    print('Directory for randomized edge files found.')
  
  return randomDirPath

def generateRandomizedFiles(edgeFile, randomDirectory, reps):
    GeneA_HOG_col = pandas.read_csv(edgeFile, sep='\t', usecols=['GeneA_HOG'])
    GeneB_HOG_col = pandas.read_csv(edgeFile, sep='\t', usecols=['GeneB_HOG'])

    for i in range(1, reps + 1):
      print("Randomizing ERCnet edge network data and writing in directory")
      rand_GeneB_HOG_col = GeneB_HOG_col.sample(frac=1).reset_index(drop=True)
      hogARandHogB = pandas.concat([GeneA_HOG_col, rand_GeneB_HOG_col], axis=1)
      newFilePath = randomDirectory + '/' + 'network_edges_' + 'rand_' + str(i) + '.tsv'
      hogARandHogB.to_csv(newFilePath, sep='\t', index=False)
      print('Finished randomization of edge file ' + str(i))
