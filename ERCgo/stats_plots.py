import argparse
import glob
import os
import shutil
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy.stats import mannwhitneyu, rankdata

def checkOutputDirectory(outPath):
  if not os.path.exists(outPath):
    print('> A directory does not exist at this path.')
    print('> Making directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print('Successful', end='\n')
  else:
    print('> Existing output directory found.')
    print('> Deleting directory to start fresh.')
    shutil.rmtree(os.path.abspath(outPath))
    print('> Making new directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print('Successful', end='\n')
    
  return outPath

def wilcoxonTest(sharedGoList):
  ercData = pandas.read_csv(sharedGoList[0], sep='\t', usecols=['Overlap_Score'])
  ercNetDist = list(ercData['Overlap_Score'])
  #print(ercDataList)

  randData = pandas.read_csv(sharedGoList[1], sep='\t', usecols=['Overlap_Score'])
  randDist = list(randData['Overlap_Score'])
  #print(randDataList)

  # Mann-Whitney U Test

  print("Sample sizes:", len(ercNetDist), len(randDist))
  print("Number of zeros in ercNetDist:", ercNetDist.count(0.0))
  print("Number of zeros in randDist:", randDist.count(0.0))
  print('\n')

  mean_ERCnet = np.mean(ercNetDist)
  mean_rand = np.mean(randDist)

  print(f"Mean of ercNetDist: {mean_ERCnet}")
  print(f"Mean of randDist: {mean_rand}")
  print('\n')

  # Compute ranks (Mann-Whitney works with ranks)
  combined_data = np.concatenate((ercNetDist, randDist))
  ranks = rankdata(combined_data)
  ranks1 = ranks[:len(ercNetDist)]
  ranks2 = ranks[len(ercNetDist):]

  mean_rank1 = np.mean(ranks1)
  mean_rank2 = np.mean(ranks2)

  print(f"Mean Rank of ercNetDist: {mean_rank1}")
  print(f"Mean Rank of randDist: {mean_rank2}")
  print('\n')

  if mean_rank1 > mean_rank2:
      print("ercNetDist has larger values on average.")
  else:
      print("randDist has larger values on average.")

  print("\nFirst 10 ranks of ercNetDist:", ranks1[:10])
  print("First 10 ranks of randDist:", ranks2[:10])

  # Mann-Whitney U Test
  u_stat, p_mw = mannwhitneyu(ercNetDist, randDist, alternative='greater')

  # Compute expected U under null hypothesis
  n1, n2 = len(ercNetDist), len(randDist)
  expected_u = (n1 * n2) / 2  # Expected value if distributions are the same
  std_u = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)  # Standard deviation under null

  print(f"\nMann-Whitney U test statistic: {u_stat}")
  print(f"Expected U under H0: {expected_u}")
  print(f"Standard deviation of U: {std_u}")
  print(f"Z-score: {(u_stat - expected_u) / std_u}")
  print(f"P-value: {p_mw}")

##Start of main stat and plotting script
print('################################################')
print('Starting Statistical Analysis of ERCgo Analysis!')
print('################################################')

##Parse CLI arguments
print('---------------------------------------------------------------------------------------------------')
print('Parsing cli arguments')
print('---------------------------------------------------------------------------------------------------')
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, metavar='dir_path',
  help='''Path to ERCnet output files and GAF file which will be used in the GO analysis''')

parser.add_argument('-o', '--output', required=True, metavar='dir_path',
  help='''Path to output folder where ERCnet GO analysis information will be written''')

args = parser.parse_args()
argsDict = vars(args)

outputDirectory = checkOutputDirectory(argsDict['output'])

##Collect analysis tables from ERCgo output directory
print('---------------------------------------------------------------------------------------------------')
print('Get GO analysis files from input ERCgo run')
print('---------------------------------------------------------------------------------------------------')
sharedGOAnalysisFilesList = glob.glob(argsDict['input'] + '/GO_ANALYSIS_*.tsv')

if len(sharedGOAnalysisFilesList) > 0:
  sharedGOAnalysisFilesList.sort()
  print('Analysis files found!')
  print('\n')

  print('---------------------------------------------------------------------------------------------------')
  print('Perform Wilcoxon Statistical Test')
  print('---------------------------------------------------------------------------------------------------')
  print('\n')
  
  wilcoxonTest(sharedGOAnalysisFilesList)
  
  print('\n')
  print('###########################################')
  print('Statistical analysis and plotting complete!')
  print('###########################################')
  print('\n')

else:
  print('No analysis files found. Please try again')  
  print('\n')
