import os
import pandas

def randomDirectory(outDirectory):
  randomDirPath = outDirectory + '/randomSets'
      
  if not os.path.exists(randomDirPath):
    print('Making directory for writing random replicates...')
    os.makedirs(randomDirPath)
    print('Successful')
  else:
    print('Directory for random replicates files found.')
  
  return randomDirPath

def generateRandomizedFiles(edgeFile, randomDirectory, reps):
    GeneA_HOG_col = pandas.read_csv(edgeFile, sep='\t', usecols=['GeneA_HOG'])
    print(GeneA_HOG_col)
    GeneB_HOG_col = pandas.read_csv(edgeFile, sep='\t', usecols=['GeneB_HOG'])
    print(GeneB_HOG_col)

    for i in range(1, reps + 1):
      print("Randomizing ERCnet edge network data and writing in directory")
      rand_GeneB_HOG_col = GeneB_HOG_col.sample(frac=1).reset_index(drop=True)
      print(rand_GeneB_HOG_col)    
      hogARandHogB = pandas.concat([GeneA_HOG_col, rand_GeneB_HOG_col], axis=1)
      newFilePath = randomDirectory + '/' + 'network_edges_' + 'rand_' + str(i) + '.tsv'
      hogARandHogB.to_csv(newFilePath, sep='\t', index=False)
      print('Finished randomization of randomized set ' + str(i))