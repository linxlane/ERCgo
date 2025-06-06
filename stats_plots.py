import argparse
import glob
import os
import shutil
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, rankdata, pearsonr, spearmanr, linregress
import matplotlib as mpl

def checkOutputDirectory(outPath):
    if not os.path.exists(outPath):
        print("> A directory does not exist at this path.")
        print("> Making directory for writing all output files...", flush=True, end="")
        os.makedirs(outPath)
        print("Successful", end="\n")
    else:
        print("> Existing output directory found.")
        print("> Deleting directory to start fresh.")
        shutil.rmtree(os.path.abspath(outPath))
        print(
            "> Making new directory for writing all output files...", flush=True, end=""
        )
        os.makedirs(outPath)
        print("Successful", end="\n")

    return outPath


def pearsonCorrelation(ercData):
    # Calculate Pearson correlation and p-value
    pearson_corr, pearson_pval = pearsonr(ercData["Overlap_Score"], ercData["P_R2"])
    print("pearson_corr: " + str(pearson_corr))
    print("pearson_pval: " + str(pearson_pval))


def spearmanCorrelation(ercData):
    # Calculate Pearson correlation and p-value
    spearman_corr, spearman_pval = spearmanr(ercData["Overlap_Score"], ercData["P_R2"])
    print("spearman_corr: " + str(spearman_corr))
    print("spearman_pval: " + str(spearman_pval))


def negLog10(col):
    values = -np.log10(col)
    return values

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

def scatterPlot(ercData):
    # print(ercData.head(50))
    # print('------------------------------------')

    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(12, 8))

    smallValue = ercData["Overlap_Score"][ercData["Overlap_Score"] != 0].min()
    replaceZerosDF = ercData.replace(to_replace=0, value=smallValue)
    # print(replaceZerosDF.head(50))

    # replaceZerosDF['negLog10'] = negLog10(replaceZerosDF['P_Pval'])

    clpClpMask = replaceZerosDF["Color"] == "Clp-Clp interaction"
    clpInterestMask = replaceZerosDF["Color"] == "Clp-Interest interaction"
    interestClpMask = replaceZerosDF["Color"] == "Interest-Clp interaction"
    interestInterestMask = replaceZerosDF["Color"] == "Interest-Interest"
    noIntAbovePoint005 = replaceZerosDF['Overlap_Score'] > 0.01
    noInterestMask = replaceZerosDF["Color"] == "Not of interest"

    clpClpDF = replaceZerosDF[clpClpMask].reset_index()
    clpInterestDF = replaceZerosDF[clpInterestMask].reset_index()
    interestClpDF = replaceZerosDF[interestClpMask].reset_index()
    interestInterestDF = replaceZerosDF[interestInterestMask].reset_index()
    noIntAbovePoint005DF = replaceZerosDF[noIntAbovePoint005].reset_index()
    

    allPoints = sns.scatterplot(
        data=replaceZerosDF[noInterestMask],
        x="Overlap_Score",
        y="P_R2",
        marker="X",
        color="#6B6B6B",
        zorder=1,
    )
    rprRprPlot = sns.scatterplot(
        data=replaceZerosDF[interestInterestMask],
        x="Overlap_Score",
        y="P_R2",
        marker="o",
        color="yellow",
        zorder=2,
    )
    rpnRpnPlot = sns.scatterplot(
        data=replaceZerosDF[interestClpMask],
        x="Overlap_Score",
        y="P_R2",
        marker="o",
        color="#005AB5",
        zorder=2,
    )
    betaBetaPlot = sns.scatterplot(
        data=replaceZerosDF[clpInterestMask],
        x="Overlap_Score",
        y="P_R2",
        marker="o",
        color="#005AB5",
        zorder=2,
    )
    alphaAlphaPlot = sns.scatterplot(
        data=replaceZerosDF[clpClpMask],
        x="Overlap_Score",
        y="P_R2",
        marker="o",
        color="#DC3220",
        zorder=2,
    )
  
    for point in range(len(noIntAbovePoint005DF)):
      label = noIntAbovePoint005DF['COMP_GENE_A'][point] + '-' + noIntAbovePoint005DF['COMP_GENE_B'][point]
      plt.text(x=noIntAbovePoint005DF['Overlap_Score'][point], y=noIntAbovePoint005DF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='#6B6B6B')
    
    for point in range(len(interestInterestDF)):
      label = getAGI(interestInterestDF['COMP_GENE_A'][point]) + '-' + getAGI(interestInterestDF['COMP_GENE_B'][point])
      plt.text(x=interestInterestDF['Overlap_Score'][point], y=interestInterestDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='yellow')

    for point in range(len(interestClpDF)):
      label = getAGI(interestClpDF['COMP_GENE_A'][point]) + '-' + getAGI(interestClpDF['COMP_GENE_B'][point])
      plt.text(x=interestClpDF['Overlap_Score'][point], y=interestClpDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='#005AB5')

    for point in range(len(clpInterestDF)):
      label = getAGI(clpInterestDF['COMP_GENE_A'][point]) + '-' + getAGI(clpInterestDF['COMP_GENE_B'][point])
      plt.text(x=clpInterestDF['Overlap_Score'][point], y=clpInterestDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='#005AB5')

    for point in range(len(clpClpDF)):
      label = getAGI(clpClpDF['COMP_GENE_A'][point]) + '-' + getAGI(clpClpDF['COMP_GENE_B'][point])
      plt.text(x=clpClpDF['Overlap_Score'][point], y=clpClpDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='#DC3220')
    

    plt.xscale("log")
    plt.xlabel('log(Overlap_Score)')
    plt.ylabel("P_R2")
    plt.title("Clp R2T Hits")
    plt.savefig('/Users/linlane/Desktop/clp_paper_GO_figure_red_blue.pdf', format = 'pdf', transparent = True) 
    plt.close()
    #plt.show()


def filterHits(ercData):
    hitsDf = ercData[ercData["P_Pval"] <= 0.0001]
    hitsDf = hitsDf[hitsDf["S_Pval"] <= 0.0001]
    hitsDf = hitsDf[hitsDf["P_R2"] >= 0.5]
    hitsDf = hitsDf[hitsDf["S_R2"] >= 0.5]
    hitsDf = hitsDf[hitsDf["Slope"] > 0]
    #nonHitsDf = ercData[ercData["P_Pval"] > 0.0001]
    #nonHitsDf = nonHitsDf[nonHitsDf["P_R2"] < 0.4]
    print('Length of hits: ' + str(len(hitsDf)))

    nonHitsFilter = ercData[(ercData['P_Pval'] > 0.0001) | (ercData['P_R2'] < 0.5)| (ercData['S_R2'] < 0.5) | (ercData['S_Pval'] > 0.0001)]
    print('Length of non-hits: ' + str(len(hitsDf)))

    print(nonHitsFilter.head())
    print('--------------------------------------')
    return hitsDf, nonHitsFilter


def plotPropKde(nonHitsProps, hitsProp):
    plt.figure()
    sns.kdeplot(nonHitsProps)
    plt.axvline(x=hitsProp, color='red', linestyle='--')
    plt.title('Proportion KDE')
    plt.xlabel('len(nonHitsSample > 0)/len(nonHitsSample)')
    plt.savefig('permutation_prop_KDE.pdf', format='pdf')

def plotMeanKde(nonHitsMeans, hitsMean):
    plt.figure()
    sns.kdeplot(nonHitsMeans)
    plt.axvline(x=hitsMean, color='red', linestyle='--')
    plt.title('Mean KDE')
    plt.xlabel('1000 means of non-hit samples')
    plt.savefig('permutation_mean_KDE.pdf', format='pdf')

def mannwhitney(hits, nonhits):
    # Perform the one-sided Mann-Whitney U test (sample1 > sample2)
    stat, p = mannwhitneyu(hits, nonhits, alternative="greater")

    print(f"Statistic: {stat}")
    print(f"P-value: {p}")


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

    print("hits division: " + str(hitsProp))
    print("hits mean: " + str(hitsMean))
    print(nonHitMeansList)
    print(nonHitPropList)
    nonHitsPropFull = len(nonHits[nonHits > 0]) / len(nonHits)
    print("Full division: " + str(nonHitsPropFull))

    return hitsMean, hitsProp, nonHitMeansList, nonHitPropList


def plotFullData(fullErcData):
    agg = ds.Canvas().points(fullErcData, "Overlap_Score", "P_R2")
    ds.tf.set_background(ds.tf.shade(agg, cmap=cc.fire), "black")


##Start of main stat and plotting script
print("################################################")
print("Starting Statistical Analysis of ERCgo Analysis!")
print("################################################")
parser = argparse.ArgumentParser()

parser.add_argument(
    "-i",
    "--input",
    required=True,
    metavar="file_path",
    help="""Path to ERCgo GO_ANALYSIS file""",
)

args = parser.parse_args()
argsDict = vars(args)

goAnalysisFilePath = argsDict["input"]

print(
    "---------------------------------------------------------------------------------------------------"
)
print("Read GO analysis file into a dataframe")
print(
    "---------------------------------------------------------------------------------------------------"
)
try:
    goAnalysisDf = pandas.read_csv(goAnalysisFilePath, sep="\t")
    print(goAnalysisDf.head())
    print("Successful!")
except:
    print(
        "There was a problem reading the provided Go analysis file. Please check your input and try again."
    )

print(
    "---------------------------------------------------------------------------------------------------"
)
print("Correlation Statistics")
print(
    "---------------------------------------------------------------------------------------------------"
)
slope = linregress(goAnalysisDf["Overlap_Score"], goAnalysisDf["P_R2"]).slope
print(f"linregress slope: {slope}")

pearsonCorrelation(goAnalysisDf)
spearmanCorrelation(goAnalysisDf)

print(
    "---------------------------------------------------------------------------------------------------"
)
print("Scatterplot")
print(
    "---------------------------------------------------------------------------------------------------"
)
#scatterPlot(goAnalysisDf)
#print("Skip")

# print('---------------------------------------------------------------------------------------------------')
# print('Datashader')
# print('---------------------------------------------------------------------------------------------------')
# plotFullData(goAnalysisDf)

print(
    "---------------------------------------------------------------------------------------------------"
)
print("Permutation Test and KDE")
print(
    "---------------------------------------------------------------------------------------------------"
)

hits, nonHits = filterHits(goAnalysisDf)
'''
# print(type(hits))
# print(type(nonHits))
print("Hit Rows")
print(len(hits))
print("Non-Hit Rows")
print(len(nonHits))
print("Hit Max")
# hitMaxLoc = hits.loc[hits['Overlap_Score'] == 0.1361111111111111]
# print(hitMaxLoc)
print(hits["Overlap_Score"].max())
print("Non-Hit Max")
print(nonHits["Overlap_Score"].max())
print("Hit Value Counts")
print(hits["Overlap_Score"].value_counts())
print("Non-Hit Value Counts")
print(nonHits["Overlap_Score"].value_counts())

hitsMean, hitsProp, nonHitsMeansList, nonHitsPropList = permutationTest(hits["Overlap_Score"], nonHits["Overlap_Score"])

plotPropKde(nonHitsPropList, hitsProp)
plotMeanKde(nonHitsMeansList, hitsMean)


# print(hits.head())
# print(nonHits.head())
# print('Mannwhitneyu test')
# mannwhitney(hits['Overlap_Score'], nonHits['Overlap_Score'])
# kde(hits, nonHits)
'''
print("\n")
print("###########################################")
print("Statistical analysis and plotting complete!")
print("###########################################")
print("\n")
