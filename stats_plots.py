import argparse
import sys
import glob
import os
import shutil
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, rankdata, pearsonr, spearmanr, linregress


def checkOutputDirectory(outPath):
  if not os.path.exists(outPath):
    print('> A directory does not exist at this path.')
    print('> Making directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print('Successful', end='\n')
  else:
    print('> Existing output directory found.')
    print('> Successful', end='\n')

  return outPath


def pearsonCorrelation(goAnalysisDf):
  # Calculate Pearson correlation and p-value
  pearson_corr, pearson_pval = pearsonr(goAnalysisDf['Overlap_Score'], goAnalysisDf['Confidence_Value'])
  print('pearson_corr: ' + str(pearson_corr))
  print('pearson_pval: ' + str(pearson_pval))
  


def spearmanCorrelation(goAnalysisDf):
  # Calculate Pearson correlation and p-value
  spearman_corr, spearman_pval = spearmanr(goAnalysisDf['Overlap_Score'], goAnalysisDf['Confidence_Value'])
  print('spearman_corr: ' + str(spearman_corr))
  print('spearman_pval: ' + str(spearman_pval))


def negLog10(col):
  values = -np.log10(col)
  return values


def scatterPlot(goAnalysisDf, argsDict):
  # print(goAnalysisDf.head(50))
  # print('------------------------------------')

  smallValue = goAnalysisDf['Overlap_Score'][goAnalysisDf['Overlap_Score'] != 0].min()
  replaceZerosDF = goAnalysisDf.replace(to_replace=0, value=smallValue)
  # print(replaceZerosDF.head(50))

  # replaceZerosDF['negLog10'] = negLog10(replaceZerosDF['P_Pval'])

  clpClpMask = replaceZerosDF['Color'] == 'Clp-Clp interaction'
  clpInterestMask = replaceZerosDF['Color'] == 'Clp-Interest interaction'
  interestClpMask = replaceZerosDF['Color'] == 'Interest-Clp interaction'
  interestInterestMask = replaceZerosDF['Color'] == 'Interest-Interest'
  noInterestMask = replaceZerosDF['Color'] == 'Not of interest'

  clpClpDF = replaceZerosDF[clpClpMask].reset_index()
  clpInterestDF = replaceZerosDF[clpInterestMask].reset_index()
  interestClpDF = replaceZerosDF[interestClpMask].reset_index()
  interestInterestDF = replaceZerosDF[interestInterestMask].reset_index()

  allPoints = sns.scatterplot(
      data=replaceZerosDF[noInterestMask],
      x='Overlap_Score',
      y='Confidence_Value',
      marker='X',
      color='green',
      zorder=1,
  )
  rprRprPlot = sns.scatterplot(
      data=replaceZerosDF[interestInterestMask],
      x='Overlap_Score',
      y='Confidence_Value',
      marker='o',
      color='yellow',
      zorder=2,
  )
  rpnRpnPlot = sns.scatterplot(
      data=replaceZerosDF[interestClpMask],
      x='Overlap_Score',
      y='Confidence_Value',
      marker='o',
      color='blue',
      zorder=2,
  )
  betaBetaPlot = sns.scatterplot(
      data=replaceZerosDF[clpInterestMask],
      x='Overlap_Score',
      y='Confidence_Value',
      marker='o',
      color='blue',
      zorder=2,
  )
  alphaAlphaPlot = sns.scatterplot(
      data=replaceZerosDF[clpClpMask],
      x='Overlap_Score',
      y='Confidence_Value',
      marker='o',
      color='red',
      zorder=2,
  )

  '''
  for point in range(len(rpnRpnDF)):
    label = getAGI(rpnRpnDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(rpnRpnDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=rpnRpnDF['Overlap_Score'][point], y=rpnRpnDF['Confidence_Value'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')

  for point in range(len(rptRptDF)):
    label = getAGI(rptRptDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(rptRptDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=rptRptDF['Overlap_Score'][point], y=rptRptDF['Confidence_Value'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')


  for point in range(len(betaBetaDF)):
    label = getAGI(betaBetaDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(betaBetaDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=betaBetaDF['Overlap_Score'][point], y=betaBetaDF['Confidence_Value'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='orange')
  
  for point in range(len(alphaAlphaDF)):
    label = getAGI(alphaAlphaDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(alphaAlphaDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=alphaAlphaDF['Overlap_Score'][point], y=alphaAlphaDF['Confidence_Value'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  
  for point in range(len(otherDF)):
    label = getAGI(otherDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(otherDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=otherDF['Overlap_Score'][point], y=otherDF['Confidence_Value'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  '''

  #plt.xscale('log')
  #plt.yscale('log')
  plt.xlabel('log(Overlap_Score)')
  plt.ylabel('log(Confidence_Value)')
  plt.title('Interactome Hits')
  #plt.show()
  plt.savefig(argsDict['output'] + '/Interactome_scatterplot.pdf', format='pdf')


def filterHits(goAnalysisDf):
    hitsDf = goAnalysisDf[goAnalysisDf['Confidence_Value'] >= 1]
    nonHitsFilter = goAnalysisDf[goAnalysisDf['Confidence_Value'] == 0]
    return hitsDf, nonHitsFilter


def plotPropKde(nonHitsProps, hitsProp):
    plt.figure()
    sns.kdeplot(nonHitsProps)
    plt.axvline(x=hitsProp, color='red', linestyle='--')
    plt.title('Proportion KDE')
    plt.xlabel('len(nonHitsSample > 0)/len(nonHitsSample)')
    plt.savefig(argsDict['output'] + '/permutation_prop_KDE.pdf', format='pdf')


def plotMeanKde(nonHitsMeans, hitsMean):
    plt.figure()
    sns.kdeplot(nonHitsMeans)
    plt.axvline(x=hitsMean, color='red', linestyle='--')
    plt.title('Mean KDE')
    plt.xlabel('1000 means of non-hit samples')
    plt.savefig(argsDict['output'] + '/permutation_mean_KDE.pdf', format='pdf')


def mannwhitney(hits, nonhits):
    # Perform the one-sided Mann-Whitney U test (sample1 > sample2)
    stat, p = mannwhitneyu(hits, nonhits, alternative='greater')

    print(f'Statistic: {stat}')
    print(f'P-value: {p}')


def permutationTest(hits, nonHits):
    hitsMean = hits.mean()
    hitsProp = len(hits[hits > 0]) / len(hits)

    nonHitMeansList = []
    nonHitPropList = []
    hitsLength = len(hits)
    for i in range(1000):
        nonHitsSample = nonHits.sample(hitsLength)
        nonHitsMean = nonHitsSample.mean()
        nonHitsPropSamp = len(nonHitsSample[nonHitsSample > 0]) / len(nonHitsSample)
        nonHitMeansList.append(nonHitsMean)
        nonHitPropList.append(nonHitsPropSamp)

    nonHitsPropFull = len(nonHits[nonHits > 0]) / len(nonHits)

    print('> Writing to file...')
    with open(argsDict['output'] + '/permutation_stats.txt', "a") as f:
      f.write('Hits Proportion: ' + str(hitsProp) + '\n')
      f.write('Hits Mean: ' + str(hitsMean) + '\n')
      f.write('Non-hits Proportion: ' + str(nonHitsPropFull) + '\n')
      f.close()
    print('> Successful!')

    return hitsMean, hitsProp, nonHitMeansList, nonHitPropList


def plotFullData(fullErcData):
    agg = ds.Canvas().points(fullErcData, 'Overlap_Score', 'Confidence_Value')
    ds.tf.set_background(ds.tf.shade(agg, cmap=cc.fire), 'black')


##Start of main stat and plotting script
print('################################################')
print('Starting Statistical Analysis of ERCgo Analysis!')
print('################################################')
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, metavar='file_path', help='''Path to ERCgo GO_ANALYSIS file''')
parser.add_argument('-m', '--mode', required=True, choices=['hits', 'full'], help='''Specifies what stats and plots to genereate depending on the input analysis type.''')
parser.add_argument('-o', '--output', required=True, metavar='dir_path', help='''Path to output directory. If this directory does not exist at runtime, stats_plots will create it.''')


args = parser.parse_args()
argsDict = vars(args)

print('---------------------------------------------------------------------------------------------------')
print('Check output directory')
print('---------------------------------------------------------------------------------------------------')
checkOutputDirectory(argsDict['output'])

print('---------------------------------------------------------------------------------------------------')
print('Read GO analysis file into a dataframe')
print('---------------------------------------------------------------------------------------------------')
goAnalysisFilePath = argsDict['input']
try:
    goAnalysisDf = pandas.read_csv(goAnalysisFilePath, sep='\t', engine='python')
    print('> Successful!')
except:
    sys.exit('There was a problem reading the provided Go analysis file. Please check your input and try again. Terminating script.')

print('---------------------------------------------------------------------------------------------------')
print('Correlation Statistics')
print('---------------------------------------------------------------------------------------------------')
print('> Generating correlation statistics...')

slope = linregress(goAnalysisDf['Overlap_Score'], goAnalysisDf['Confidence_Value']).slope
pearson_corr, pearson_pval = pearsonr(goAnalysisDf['Overlap_Score'], goAnalysisDf['Confidence_Value'])
spearman_corr, spearman_pval = spearmanr(goAnalysisDf['Overlap_Score'], goAnalysisDf['Confidence_Value'])

print('> Writing to file...')
with open(argsDict['output'] + '/correlation_stats.txt', "w") as f:
  f.write('All statistics calculated using scipy stats:\n')
  f.write('linregress_slope: ' + str(slope) + '\n')
  f.write('pearson_corr: ' + str(pearson_corr) + '\n')
  f.write('pearson_pval: ' + str(pearson_pval) + '\n')
  f.write('spearman_corr: ' + str(spearman_corr) + '\n')
  f.write('spearman_pval: ' + str(spearman_pval) + '\n')
print('> Successful!')

print('---------------------------------------------------------------------------------------------------')
print('Scatterplot')
print('---------------------------------------------------------------------------------------------------')
if argsDict['mode'] == 'hits':
  print('> Generating scatterplot...')
  scatterPlot(goAnalysisDf, argsDict)
  print('> Successful!')
else:
  print('> Full mode activated. Skip scatterplot.')

print('---------------------------------------------------------------------------------------------------')
print('Permutation Test and KDE')
print('---------------------------------------------------------------------------------------------------')
if argsDict['mode'] == 'full':
  hits, nonHits = filterHits(goAnalysisDf)
  # print(type(hits))
  # print(type(nonHits))
  print('> Calculating basic stats for hits and non-hits')
  print('> Writing to file...')
  with open(argsDict['output'] + '/permutation_stats.txt', "w") as f:
    f.write('Hits Rows: ' + str(len(hits)) + '\n')
    f.write('Non-hit Rows: ' + str(len(nonHits)) + '\n')
    f.write('Hit Max: ' + str(hits['Overlap_Score'].max()) + '\n')
    f.write('Non-hit Max: ' + str(nonHits['Overlap_Score'].max()) + '\n')
    f.close()
  print('> Successful!')

  hitsMean, hitsProp, nonHitsMeansList, nonHitsPropList = permutationTest(hits['Overlap_Score'], nonHits['Overlap_Score'])

  plotPropKde(nonHitsPropList, hitsProp)
  plotMeanKde(nonHitsMeansList, hitsMean)
else:
   print('> Hits mode activated. Skip permutation test.')

print()
print('###########################################')
print('Statistical analysis and plotting complete!')
print('###########################################')
