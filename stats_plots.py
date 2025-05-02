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

def getAGI(compID, agiDict):
  agi = agiDict[compID]
  return agi

def scatterPlot(sharedGoTable, agiDict):
  ercData = pandas.read_csv(sharedGoTable, sep='\t')

  alphaAlphaMask = ercData['Color'] == 'Alpha-Alpha'
  betaBetaMask = ercData['Color'] == 'Beta-Beta'
  rpnRpnMask = ercData['Color'] == 'RPN-RPN'
  rptRptMask = ercData['Color'] == 'RPT-RPT'
  otherMask = ercData['Color'].str.contains('Other')
  noInterestMask = ercData['Color'] == 'Not of interest'

  alphaAlphaDF = ercData[alphaAlphaMask].reset_index()
  betaBetaDF = ercData[betaBetaMask].reset_index()
  rpnRpnDF = ercData[rpnRpnMask].reset_index()
  rptRptDF = ercData[rptRptMask].reset_index()
  otherDF = ercData[otherMask].reset_index()



  allPoints = sns.scatterplot(data=ercData[noInterestMask], x='Overlap_Score', y='P_R2', marker='X', color='green', zorder=1)
  rprRprPlot = sns.scatterplot(data=ercData[rptRptMask], x='Overlap_Score', y='P_R2', marker='o', color='yellow',  zorder=2)
  rpnRpnPlot = sns.scatterplot(data=ercData[rpnRpnMask], x='Overlap_Score', y='P_R2', marker='o', color='purple',  zorder=2)
  betaBetaPlot = sns.scatterplot(data=ercData[betaBetaMask], x='Overlap_Score', y='P_R2', marker='o', color='orange',  zorder=2)
  alphaAlphaPlot = sns.scatterplot(data=ercData[alphaAlphaMask], x='Overlap_Score', y='P_R2', marker='o', color='blue', zorder=2)
  otherPlot = sns.scatterplot(data=ercData[otherMask], x='Overlap_Score', y='P_R2', marker='o', color='red', zorder=2)

  '''
  for point in range(len(rpnRpnDF)):
    label = getAGI(rpnRpnDF['COMP_GENE_A'][point]) + '-' + getAGI(rpnRpnDF['COMP_GENE_B'][point])
    plt.text(x=rpnRpnDF['Overlap_Score'][point], y=rpnRpnDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')

  for point in range(len(rptRptDF)):
    label = getAGI(rptRptDF['COMP_GENE_A'][point]) + '-' + getAGI(rptRptDF['COMP_GENE_B'][point])
    plt.text(x=rptRptDF['Overlap_Score'][point], y=rptRptDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')


  for point in range(len(betaBetaDF)):
    label = getAGI(betaBetaDF['COMP_GENE_A'][point]) + '-' + getAGI(betaBetaDF['COMP_GENE_B'][point])
    plt.text(x=betaBetaDF['Overlap_Score'][point], y=betaBetaDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='orange')
  
  for point in range(len(alphaAlphaDF)):
    label = getAGI(alphaAlphaDF['COMP_GENE_A'][point]) + '-' + getAGI(alphaAlphaDF['COMP_GENE_B'][point])
    plt.text(x=alphaAlphaDF['Overlap_Score'][point], y=alphaAlphaDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  '''

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

parser.add_argument('-o', '--output', required=True, metavar='dir_path',
  help='''Path to output folder where ERCnet GO analysis information will be written''')

args = parser.parse_args()
argsDict = vars(args)

outputDirectory = checkOutputDirectory(argsDict['output'])

agiDict = {'ARAT_AT2G05840.1_PACid_19639196' : 'Alpha1_AT2G05840',
            'ARAT_AT5G35590.1_PACid_19666663' : 'Alpha1_AT5G35590',
            'ARAT_AT1G16470.1_PACid_19651668' : 'Alpha2_AT1G16470',
            'ARAT_AT1G79210.1_PACid_19650105' : 'Alpha2_AT1G79210',
            'ARAT_AT3G22110.1_PACid_19663555' : 'Alpha3_AT3G22110',
            'ARAT_AT4G15165.1_PACid_19648727' : 'Alpha3_AT4G15165',
            'ARAT_AT3G51260.1_PACid_19664689' : 'Alpha4_AT3G51260',
            'ARAT_AT5G66140.1_PACid_19672619' : 'Alpha4_AT5G66140',
            'ARAT_AT1G53850.1_PACid_19652641' : 'Alpha5_AT1G53850',
            'ARAT_AT3G14290.1_PACid_19658710' : 'Alpha5_AT3G14290',
            'ARAT_AT1G47250.1_PACid_19658242' : 'Alpha_B26_AT1G47250',
            'ARAT_AT5G42790.1_PACid_19670251' : 'Alpha6_AT5G42790',
            'ARAT_AT2G27020.1_PACid_19640509' : 'Alpha7_AT2G27020',
            'ARAT_AT4G31300.2_PACid_19646330' : 'Beta1_AT4G31300',
            'ARAT_AT3G27430.2_PACid_19660822' : 'Beta2_AT3G27430',
            'ARAT_AT5G40580.1_PACid_19671607' : 'Beta2_AT5G40580',
            'ARAT_AT1G21720.1_PACid_19649188' : 'Beta3_AT1G21720',
            'ARAT_AT1G77440.1_PACid_19656196' : 'Beta3_AT1G77440',
            'ARAT_AT3G22630.1_PACid_19662898' : 'Beta4_AT3G22630',
            'ARAT_AT4G14800.2_PACid_19644251' : 'Beta4_AT4G14800',
            'ARAT_AT1G13060.2_PACid_19652499' : 'Beta5_AT1G13060',
            'ARAT_AT3G26340.1_PACid_19658790' : 'Beta5_AT3G26340',
            'ARAT_AT3G60820.1_PACid_19663843' : 'Beta6_AT3G60820',
            'ARAT_AT1G56450.1_PACid_19658067' : 'Beta7_AT1G56450',
            'ARAT_AT4G38630.1_PACid_19646282' : 'RPN10_AT4G38630',
            'ARAT_AT5G23540.1_PACid_19666107' : 'RPN11_AT5G23540',
            'ARAT_AT1G64520.1_PACid_19658060' : 'RPN12_AT1G64520',
            'ARAT_AT5G42040.1_PACid_19665799' : 'RPN12_AT5G42040',
            'ARAT_AT2G20580.1_PACid_19641906' : 'RPN1_AT2G20580',
            'ARAT_AT4G28470.1_PACid_19645981' : 'RPN1_AT4G28470',
            'ARAT_AT1G04810.1_PACid_19651161' : 'RPN2_AT1G04810',
            'ARAT_AT2G32730.1_PACid_19639179' : 'RPN2_AT2G32730',
            'ARAT_AT1G20200.1_PACid_19651422' : 'RPN3_AT1G20200',
            'ARAT_AT1G75990.1_PACid_19652198' : 'RPN3_AT1G75990',
            'ARAT_AT5G09900.3_PACid_19665897' : 'RPN5_AT5G09900',
            'ARAT_AT5G64760.1_PACid_19668432' : 'RPN5_AT5G64760',
            'ARAT_AT1G29150.1_PACid_19650043' : 'RPN6_AT1G29150',
            'ARAT_AT4G24820.1_PACid_19647825' : 'RPN7_AT4G24820',
            'ARAT_AT3G11270.1_PACid_19658597' : 'RPN8_AT3G11270',
            'ARAT_AT5G05780.1_PACid_19670513' : 'RPN8_AT5G05780',
            'ARAT_AT4G19006.1_PACid_19643829' : 'RPN9_AT4G19006',
            'ARAT_AT5G45620.1_PACid_19668970' : 'RPN9_AT5G45620',
            'ARAT_AT1G53750.1_PACid_19657550' : 'RPT1_AT1G53750',
            'ARAT_AT1G53780.2_PACid_19652140' : 'RPT1_AT1G53780',
            'ARAT_AT2G20140.1_PACid_19639123' : 'RPT2_AT2G20140',
            'ARAT_AT4G29040.1_PACid_19646587' : 'RPT2_AT4G29040',
            'ARAT_AT5G58290.1_PACid_19666361' : 'RPT3_AT5G58290',
            'ARAT_AT1G45000.1_PACid_19650223' : 'RPT4_AT1G45000',
            'ARAT_AT5G43010.1_PACid_19670324' : 'RPT4_AT5G43010',
            'ARAT_AT1G09100.1_PACid_19653012' : 'RPT5_AT1G09100',
            'ARAT_AT3G05530.1_PACid_19663123' : 'RPT5_AT3G05530',
            'ARAT_AT5G19990.1_PACid_19667687' : 'RPT6_AT5G19990',
            'ARAT_AT5G20000.1_PACid_19669291' : 'RPT6_AT5G20000'}

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
  print('GO Score vs P_R2 Scatter Plot')
  print('---------------------------------------------------------------------------------------------------')
  print('\n')
  scatterPlot(argsDict['input'] + '/GO_ANALYSIS_ERCnet_Network.tsv', agiDict)

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
