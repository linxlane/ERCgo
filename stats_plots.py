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

def pearsonCorrelation(ercData):
  # Calculate Pearson correlation and p-value
  pearson_corr, pearson_pval = pearsonr(ercData['Overlap_Score'], ercData['P_Pval'])
  print('pearson_corr: ' + str(pearson_corr))
  print('pearson_pval: ' + str(pearson_pval))

def negLog10(col):
  values = -np.log10(col)
  return values

def scatterPlot(ercData):
  #print(ercData.head(50))
  #print('------------------------------------')

  smallValue = ercData['Overlap_Score'][ercData['Overlap_Score'] != 0].min()
  replaceZerosDF = ercData.replace(to_replace=0, value=smallValue)
  #print(replaceZerosDF.head(50))

  replaceZerosDF['negLog10'] = negLog10(replaceZerosDF['P_Pval'])

  alphaAlphaMask = replaceZerosDF['Color'] == 'Alpha-Alpha'
  betaBetaMask = replaceZerosDF['Color'] == 'Beta-Beta'
  rpnRpnMask = replaceZerosDF['Color'] == 'RPN-RPN'
  rptRptMask = replaceZerosDF['Color'] == 'RPT-RPT'
  otherMask = replaceZerosDF['Color'].str.contains('Other')
  noInterestMask = replaceZerosDF['Color'] == 'Not of interest'

  alphaAlphaDF = replaceZerosDF[alphaAlphaMask].reset_index()
  betaBetaDF = replaceZerosDF[betaBetaMask].reset_index()
  rpnRpnDF = replaceZerosDF[rpnRpnMask].reset_index()
  rptRptDF = replaceZerosDF[rptRptMask].reset_index()
  otherDF = replaceZerosDF[otherMask].reset_index()



  allPoints = sns.scatterplot(data=replaceZerosDF[noInterestMask], x='Overlap_Score', y='negLog10', marker='X', color='green', zorder=1)
  rprRprPlot = sns.scatterplot(data=replaceZerosDF[rptRptMask], x='Overlap_Score', y='negLog10', marker='o', color='yellow',  zorder=2)
  rpnRpnPlot = sns.scatterplot(data=replaceZerosDF[rpnRpnMask], x='Overlap_Score', y='negLog10', marker='o', color='purple',  zorder=2)
  betaBetaPlot = sns.scatterplot(data=replaceZerosDF[betaBetaMask], x='Overlap_Score', y='negLog10', marker='o', color='orange',  zorder=2)
  alphaAlphaPlot = sns.scatterplot(data=replaceZerosDF[alphaAlphaMask], x='Overlap_Score', y='negLog10', marker='o', color='blue', zorder=2)
  otherPlot = sns.scatterplot(data=replaceZerosDF[otherMask], x='Overlap_Score', y='negLog10', marker='o', color='red', zorder=2)

  '''
  for point in range(len(rpnRpnDF)):
    label = getAGI(rpnRpnDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(rpnRpnDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=rpnRpnDF['Overlap_Score'][point], y=rpnRpnDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')

  for point in range(len(rptRptDF)):
    label = getAGI(rptRptDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(rptRptDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=rptRptDF['Overlap_Score'][point], y=rptRptDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')


  for point in range(len(betaBetaDF)):
    label = getAGI(betaBetaDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(betaBetaDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=betaBetaDF['Overlap_Score'][point], y=betaBetaDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='orange')
  
  for point in range(len(alphaAlphaDF)):
    label = getAGI(alphaAlphaDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(alphaAlphaDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=alphaAlphaDF['Overlap_Score'][point], y=alphaAlphaDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  
  for point in range(len(otherDF)):
    label = getAGI(otherDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(otherDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=otherDF['Overlap_Score'][point], y=otherDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  '''

  plt.xscale('log')
  
  plt.ylabel('-log10(P_Pval)')

  plt.title('BXB')

  plt.show()

##Start of main stat and plotting script
print('################################################')
print('Starting Statistical Analysis of ERCgo Analysis!')
print('################################################')
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, metavar='file_path',
  help='''Path to ERCgo GO_ANALYSIS file''')

args = parser.parse_args()
argsDict = vars(args)

goAnalysisFilePath = argsDict['input']

print('---------------------------------------------------------------------------------------------------')
print('Read GO analysis file into a dataframe')
print('---------------------------------------------------------------------------------------------------')
try: 
  goAnalysisDf = pandas.read_csv(goAnalysisFilePath, sep='\t')
  print('Successful!')
except:
  print('There was a problem reading the provided Go analysis file. Please check your input and try again.')

print('---------------------------------------------------------------------------------------------------')
print('Correlation Statistics')
print('---------------------------------------------------------------------------------------------------')
pearsonCorrelation(goAnalysisDf)

print('---------------------------------------------------------------------------------------------------')
print('Scatterplot')
print('---------------------------------------------------------------------------------------------------')
scatterPlot(goAnalysisDf)

print('\n')
print('###########################################')
print('Statistical analysis and plotting complete!')
print('###########################################')
print('\n')
