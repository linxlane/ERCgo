import argparse
import glob
import os
import shutil
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from scipy.stats import mannwhitneyu, rankdata, pearsonr

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

def trendline(ercData):
  # Calculate Pearson correlation and p-value
  pearson_corr, pearson_pval = pearsonr(ercData['Overlap_Score'], ercData['P_R2'])
  print('pearson_corr: ' + str(pearson_corr))
  print('pearson_pval: ' + str(pearson_pval))

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

args = parser.parse_args()
argsDict = vars(args)

##Collect analysis tables from ERCgo output directory
print('---------------------------------------------------------------------------------------------------')
print('Get GO analysis files from input ERCgo run')
print('---------------------------------------------------------------------------------------------------')
sharedGOAnalysisFilesList = glob.glob(argsDict['input'] + '/GO_ANALYSIS_*.tsv')

if len(sharedGOAnalysisFilesList) > 0:
  sharedGOAnalysisFilesList.sort()
  print(sharedGOAnalysisFilesList)
  print('Analysis files found!')
  print('\n')


  df = pandas.read_csv(sharedGOAnalysisFilesList[0], sep='\t')
  trendline(df)
  
  print('\n')
  print('###########################################')
  print('Statistical analysis and plotting complete!')
  print('###########################################')
  print('\n')

else:
  print('No analysis files found. Please try again')  
  print('\n')
