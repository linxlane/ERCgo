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
  ercData = pandas.read_csv(sharedGoList[0], sep='\t', usecols=['overlap_log'])
  ercNetDist = list(ercData['overlap_log'])
  #print(ercDataList)

  randData = pandas.read_csv(sharedGoList[1], sep='\t', usecols=['overlap_log'])
  randDist = list(randData['overlap_log'])
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

def getAGI(compID):
  agiDict = {'AT1G49970' : 'CLPR1',
              'AT1G12410'	: 'CLPR2',
              'AT1G09130'	: 'CLPR3',
              'AT4G17040'	: 'CLPR4',
              'AT1G66670'	: 'CLPP3',
              'AT5G45390'	: 'CLPP4',
              'AT1G02560'	: 'CLPP5',
              'AT1G11750'	: 'CLPP6',
              'AT5G51070'	: 'CLPD',
              'AT1G68660'	: 'CLPS',
              'AT5G50920'	: 'CLPC1',
              'AT3G48870'	: 'CLPC2',
              'AT4G25370'	: 'CLPT1',
              'AT4G12060'	: 'CLPT2',
              'ATCG00670'	: 'CLPP1',
              'AT3G14240' : 'SBT1.5',
              'AT1G19370' : 'mp', #membrane protein;no name 
              'AT3G47470' : 'LHCA4', 
              'AT3G23620' : 'ARPF2', 
              'AT4G11960' : 'PGRL1B',
              'AT1G08130' : 'LIG1', 
              'AT4G12800' : 'PSAL;RALFL24', 
              'AT5G01590' : 'TIC56', 
              'AT4G38590' : 'BGAL14', 
              'AT1G09850' : 'XBCP3', 
              'AT1G52220' : 'CURT1C', 
              'AT2G40360' : 'ATPEIP1;ATPEP1;BOP1', 
              'AT3G09050' : '8a7o', #8-amino-7-oxononanoate synthase 
              'AT2G47990' : 'EDA13;EDA19;SWA1', 
              'AT1G26090' : 'P-loop', #P-loop containing nucleoside triphosphate hydrolases superfamily protein
              'AT4G23940' : 'ARC1;FTSHI1',
              'AT3G60830' : 'ARP7;ATARP7', 
              'AT1G26460' : 'TPR', #Tetratricopeptide repeat (TPR)-like superfamily protein 
              'AT1G06950' : 'ATTIC110;TIC110', 
              'AT2G04270' : 'RNASE E', # RNASE E;RNASE E/G-LIKE;RNE;RNEE/G
              'AT1G10510' : 'emb2004', 
              'AT5G53080' : 'WTG1',
              'AT3G19800' : 'DUF177B', 
              'AT3G12380' : 'ARP5;ATARP5', 
              'AT1G36320' : 'CDB1L'
              }
  
  agi = agiDict[compID]
  return agi

def negLog10(col):
  values = -np.log10(col)
  return values

def logZero(col):
  values = -np.log10(col, out=np.zeros_like(col, dtype=np.float64), where=(col!=0))
  return values

def scatterPlot(sharedGoTable):
  ercData = pandas.read_csv(sharedGoTable, sep='\t')
  ercData['overlap_log'] = logZero(ercData['Overlap_Score'])
  ercData['negLog10'] = negLog10(ercData['P_Pval'])
  print(ercData.head(50))
  
  clpClpMask = ercData['Color'] == 'Clp-Clp interaction'
  clpInterestMask = ercData['Color'] == 'Clp-interest interaction'
  interestInterestMask = ercData['Color'] == 'Interest-Interest interaction'
  noIntAbovePoint005 = ercData['overlap_log'] > 0.005
  noInterestMask = ercData['Color'] == 'Not of interest'

  clpClpDF = ercData[clpClpMask].reset_index()
  clpIntDF = ercData[clpInterestMask].reset_index()
  intIntDF = ercData[interestInterestMask].reset_index()
  noIntAbovePoint005DF = ercData[noIntAbovePoint005].reset_index()


  allPoints = sns.scatterplot(data=ercData[noInterestMask], x='overlap_log', y='negLog10', marker='X', color='green', zorder=1)
  interestInterest = sns.scatterplot(data=ercData[interestInterestMask], x='overlap_log', y='negLog10', marker='o', color='purple',  zorder=2)
  clpInterest = sns.scatterplot(data=ercData[clpInterestMask], x='overlap_log', y='negLog10', marker='o', color='orange',  zorder=2)
  clpClp = sns.scatterplot(data=ercData[clpClpMask], x='overlap_log', y='negLog10', marker='o', color='blue', zorder=2)

  '''
  for point in range(len(noIntAbovePoint005DF)):
    label = getAGI(noIntAbovePoint005DF['COMP_GENE_A'][point]) + '-' + getAGI(noIntAbovePoint005DF['COMP_GENE_B'][point])
    plt.text(x=noIntAbovePoint005DF['overlap_log'][point], y=noIntAbovePoint005DF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='green')
  '''
  
  for point in range(len(intIntDF)):
    label = getAGI(intIntDF['COMP_GENE_A'][point]) + '-' + getAGI(intIntDF['COMP_GENE_B'][point])
    plt.text(x=intIntDF['overlap_log'][point], y=intIntDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')

  for point in range(len(clpIntDF)):
    label = getAGI(clpIntDF['COMP_GENE_A'][point]) + '-' + getAGI(clpIntDF['COMP_GENE_B'][point])
    plt.text(x=clpIntDF['overlap_log'][point], y=clpIntDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='orange')
  
  for point in range(len(clpClpDF)):
    label = getAGI(clpClpDF['COMP_GENE_A'][point]) + '-' + getAGI(clpClpDF['COMP_GENE_B'][point])
    plt.text(x=clpClpDF['overlap_log'][point], y=clpClpDF['negLog10'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')

  '''
  # Calculate Pearson correlation and p-value
  pearson_corr, pearson_pval = pearsonr(ercData['overlap_log'], ercData['negLog10'])
  print(pearson_corr)
  print(pearson_pval)
  '''
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
  help='''Path to ERCgo output''')

#parser.add_argument('-o', '--output', required=True, metavar='dir_path',
#  help='''Path to output folder where ERCgo stats analysis and plots will be written''')

args = parser.parse_args()
argsDict = vars(args)

#outputDirectory = checkOutputDirectory(argsDict['output'])

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
  print('GO Score vs negLog10 Scatter Plot')
  print('---------------------------------------------------------------------------------------------------')
  print('\n')
  scatterPlot(argsDict['input'] + '/GO_ANALYSIS_ERCnet_Network.tsv')
  
  #print('---------------------------------------------------------------------------------------------------')
  #print('Perform Wilcoxon Statistical Test')
  #print('---------------------------------------------------------------------------------------------------')
  #print('\n')
  
  #wilcoxonTest(sharedGOAnalysisFilesList)
  
  print('\n')
  print('###########################################')
  print('Statistical analysis and plotting complete!')
  print('###########################################')
  print('\n')

else:
  print('No analysis files found. Please try again')  
  print('\n')
