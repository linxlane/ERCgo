import argparse
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


def scatterPlot(ercData):
    # print(ercData.head(50))
    # print('------------------------------------')

    smallValue = ercData["Overlap_Score"][ercData["Overlap_Score"] != 0].min()
    replaceZerosDF = ercData.replace(to_replace=0, value=smallValue)
    # print(replaceZerosDF.head(50))

    # replaceZerosDF['negLog10'] = negLog10(replaceZerosDF['P_Pval'])

    clpClpMask = replaceZerosDF["Color"] == "Clp-Clp interaction"
    clpInterestMask = replaceZerosDF["Color"] == "Clp-Interest interaction"
    interestClpMask = replaceZerosDF["Color"] == "Interest-Clp interaction"
    interestInterestMask = replaceZerosDF["Color"] == "Interest-Interest"
    noInterestMask = replaceZerosDF["Color"] == "Not of interest"

    clpClpDF = replaceZerosDF[clpClpMask].reset_index()
    clpInterestDF = replaceZerosDF[clpInterestMask].reset_index()
    interestClpDF = replaceZerosDF[interestClpMask].reset_index()
    interestInterestDF = replaceZerosDF[interestInterestMask].reset_index()
    

    allPoints = sns.scatterplot(
        data=replaceZerosDF[noInterestMask],
        x="Overlap_Score",
        y="P_R2",
        marker="X",
        color="green",
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
        color="blue",
        zorder=2,
    )
    betaBetaPlot = sns.scatterplot(
        data=replaceZerosDF[clpInterestMask],
        x="Overlap_Score",
        y="P_R2",
        marker="o",
        color="blue",
        zorder=2,
    )
    alphaAlphaPlot = sns.scatterplot(
        data=replaceZerosDF[clpClpMask],
        x="Overlap_Score",
        y="P_R2",
        marker="o",
        color="red",
        zorder=2,
    )

    """
  for point in range(len(rpnRpnDF)):
    label = getAGI(rpnRpnDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(rpnRpnDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=rpnRpnDF['Overlap_Score'][point], y=rpnRpnDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')

  for point in range(len(rptRptDF)):
    label = getAGI(rptRptDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(rptRptDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=rptRptDF['Overlap_Score'][point], y=rptRptDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='purple')


  for point in range(len(betaBetaDF)):
    label = getAGI(betaBetaDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(betaBetaDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=betaBetaDF['Overlap_Score'][point], y=betaBetaDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='orange')
  
  for point in range(len(alphaAlphaDF)):
    label = getAGI(alphaAlphaDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(alphaAlphaDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=alphaAlphaDF['Overlap_Score'][point], y=alphaAlphaDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  
  for point in range(len(otherDF)):
    label = getAGI(otherDF['COMP_GENE_A'][point], agiDict) + '-' + getAGI(otherDF['COMP_GENE_B'][point], agiDict)
    plt.text(x=otherDF['Overlap_Score'][point], y=otherDF['P_R2'][point], s=label, horizontalalignment='center', verticalalignment='bottom', color='blue')
  """

    plt.xscale("log")
    plt.xlabel('log(Overlap_Score)')
    plt.ylabel("P_R2")
    plt.title("Clp R2T Hits")
    plt.show()


def filterHits(ercData):
    hitsDf = ercData[ercData["P_Pval"] <= 0.0001]
    hitsDf = hitsDf[hitsDf["P_R2"] >=0.4]
    hitsDf = hitsDf[hitsDf["Slope"] > 0]
    #nonHitsDf = ercData[ercData["P_Pval"] > 0.0001]
    #nonHitsDf = nonHitsDf[nonHitsDf["P_R2"] < 0.4]

    nonHitsFilter = ercData[(ercData['P_Pval'] > 0.0001) | (ercData['P_R2'] < 0.4)]
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
print("Skip")

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

print("\n")
print("###########################################")
print("Statistical analysis and plotting complete!")
print("###########################################")
print("\n")
