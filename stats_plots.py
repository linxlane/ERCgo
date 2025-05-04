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

def getAGI(compID, agiDict):
  try:
    agi = agiDict[compID]
  except KeyError:
    agi = 'Other'

  return agi

def trendline(ercData):
  # Calculate Pearson correlation and p-value
  pearson_corr, pearson_pval = pearsonr(ercData['Overlap_Score'], ercData['P_R2'])
  print('pearson_corr: ' + str(pearson_corr))
  print('pearson_pval: ' + str(pearson_pval))

def negLog10(col):
  values = -np.log10(col)
  return values

def scatterPlot(ercData, agiDict):
  #print(ercData.head(50))
  #print('------------------------------------')

  smallValue = 0.00001
  replaceZerosDF = ercData.replace(to_replace=0, value=smallValue)
  #print(replaceZerosDF.head(50))

  replaceZerosDF['negLog10'] = negLog10(replaceZerosDF['P_R2'])

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
  
  plt.ylabel('-log10(P_R2)')

  plt.title('BXB')

  plt.show()

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

#parser.add_argument('-o', '--output', required=True, metavar='dir_path',
#  help='''Path to output folder where ERCnet GO analysis information will be written''')

args = parser.parse_args()
argsDict = vars(args)

#outputDirectory = checkOutputDirectory(argsDict['output'])

agiDict = {'AT2G05840' : 'Alpha1_AT2G05840',
            'AT5G35590' : 'Alpha1_AT5G35590',
            'AT1G16470' : 'Alpha2_AT1G16470',
            'AT1G79210' : 'Alpha2_AT1G79210',
            'AT3G22110' : 'Alpha3_AT3G22110',
            'AT4G15165' : 'Alpha3_AT4G15165',
            'AT3G51260' : 'Alpha4_AT3G51260',
            'AT5G66140' : 'Alpha4_AT5G66140',
            'AT1G53850' : 'Alpha5_AT1G53850',
            'AT3G14290' : 'Alpha5_AT3G14290',
            'AT1G47250' : 'Alpha_B26_AT1G47250',
            'AT5G42790' : 'Alpha6_AT5G42790',
            'AT2G27020' : 'Alpha7_AT2G27020',
            'AT4G31300' : 'Beta1_AT4G31300',
            'AT3G27430' : 'Beta2_AT3G27430',
            'AT5G40580' : 'Beta2_AT5G40580',
            'AT1G21720' : 'Beta3_AT1G21720',
            'AT1G77440' : 'Beta3_AT1G77440',
            'AT3G22630' : 'Beta4_AT3G22630',
            'AT4G14800' : 'Beta4_AT4G14800',
            'AT1G13060' : 'Beta5_AT1G13060',
            'AT3G26340' : 'Beta5_AT3G26340',
            'AT3G60820' : 'Beta6_AT3G60820',
            'AT1G56450' : 'Beta7_AT1G56450',
            'AT4G38630' : 'RPN10_AT4G38630',
            'AT5G23540' : 'RPN11_AT5G23540',
            'AT1G64520' : 'RPN12_AT1G64520',
            'AT5G42040' : 'RPN12_AT5G42040',
            'AT2G20580' : 'RPN1_AT2G20580',
            'AT4G28470' : 'RPN1_AT4G28470',
            'AT1G04810' : 'RPN2_AT1G04810',
            'AT2G32730' : 'RPN2_AT2G32730',
            'AT1G20200' : 'RPN3_AT1G20200',
            'AT1G75990' : 'RPN3_AT1G75990',
            'AT5G09900' : 'RPN5_AT5G09900',
            'AT5G64760' : 'RPN5_AT5G64760',
            'AT1G29150' : 'RPN6_AT1G29150',
            'AT4G24820' : 'RPN7_AT4G24820',
            'AT3G11270' : 'RPN8_AT3G11270',
            'AT5G05780' : 'RPN8_AT5G05780',
            'AT4G19006' : 'RPN9_AT4G19006',
            'AT5G45620' : 'RPN9_AT5G45620',
            'AT1G53750' : 'RPT1_AT1G53750',
            'AT1G53780' : 'RPT1_AT1G53780',
            'AT2G20140' : 'RPT2_AT2G20140',
            'AT4G29040' : 'RPT2_AT4G29040',
            'AT5G58290' : 'RPT3_AT5G58290',
            'AT1G45000' : 'RPT4_AT1G45000',
            'AT5G43010' : 'RPT4_AT5G43010',
            'AT1G09100' : 'RPT5_AT1G09100',
            'AT3G05530' : 'RPT5_AT3G05530',
            'AT5G19990' : 'RPT6_AT5G19990',
            'AT5G20000' : 'RPT6_AT5G20000'}

##Collect analysis tables from ERCgo output directory
print('---------------------------------------------------------------------------------------------------')
print('Get GO analysis files from input ERCgo run')
print('---------------------------------------------------------------------------------------------------')
sharedGOAnalysisFilesList = glob.glob(argsDict['input'] + '/GO_ANALYSIS_*.tsv')

if len(sharedGOAnalysisFilesList) > 0:
  sharedGOAnalysisFilesList.sort()
  print('Analysis files found!')
  print('\n')
  ercData = pandas.read_csv(sharedGOAnalysisFilesList[0], sep='\t')

  trendline(ercData)

  print('---------------------------------------------------------------------------------------------------')
  print('GO Score vs P_R2 Scatter Plot')
  print('---------------------------------------------------------------------------------------------------')
  print('\n')
  scatterPlot(ercData, agiDict)

  print('---------------------------------------------------------------------------------------------------')
  print('Perform Wilcoxon Statistical Test')
  print('---------------------------------------------------------------------------------------------------')
  print('\n')
  
  #wilcoxonTest(sharedGOAnalysisFilesList)
  
  print('\n')
  print('###########################################')
  print('Statistical analysis and plotting complete!')
  print('###########################################')
  print('\n')

else:
  print('No analysis files found. Please try again')  
  print('\n')
