import glob
import sys
import os
import shutil
import hog_comp_ids
import pandas
import re

def findEdgeFile(directory):
  search_result = glob.glob(directory + '/*edges*')
  if len(search_result) == 0:
    sys.exit('  > ERCnet edge file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print(' > ERCnet edge file found.')
    return search_result[0]

  else:
    sys.exit('  > Too many edge files found. Terminating ERCgo.')


def findVerticesFile(directory):
  search_result = glob.glob(directory + '/*vertices*')
  if len(search_result) == 0:
    sys.exit('  > ERCnet vertices file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print(' > ERCnet vertices file found.')
    return search_result[0]

  else:
    sys.exit('  > Too many vertices files found. Terminating ERCgo.')


def findErcResultsFile(directory):
  search_result = glob.glob(directory + '/ERC_results_*')
  if len(search_result) == 0:
    sys.exit('  > ERCnet results file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print(' > ERCnet results file found.')
    return search_result[0]

  else:
    sys.exit('  > Too many results files found. Terminating ERCgo.')


def findGafFile(directory):
  search_result = glob.glob(directory + '/*.gaf')
  if len(search_result) == 0:
    sys.exit('  > GAF file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print(' > GAF file found.')
    return search_result[0]

  else:
    sys.exit('  > Too many GAF files found. Terminating ERCgo.')


def checkOutputDirectory(outPath):
  if not os.path.exists(outPath):
    print(' > A directory does not exist at this path.')
    print(' > Making directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print(' Successful')
  else:
    print(' > Existing output directory found.')
    print(' > Deleting directory to start fresh.')
    shutil.rmtree(os.path.abspath(outPath))
    print(' > Making new directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print(' Successful')
    
  return outPath


def formatErcNetDataHits(argsDict, intermediateFilesPath, edgeFilePath, verticesFilePath):
  print('1. Generating dictionary of {HOG_ID:Comprehensive_ID} to convert edge file HOG IDs to comprehensive IDs', flush=True)
  ##Create dictionary of hog id equivalent comp ids for edge file conversion
  hogCompDict = hog_comp_ids.generateHogCompDict(verticesFilePath)
  print(' > DONE', flush=True)


  ##Convert HOG and comp ids
  print('2. Convert gene pairs for analysis', flush=True)
  print('  > Generate table of [HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]', flush=True)
  hogCompPath = intermediateFilesPath + '/[HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]_FULL_' + argsDict['job_name'] + '.tsv'
  #Create and write tsv that contains equivalent HOG-COMP id values
  geneLookupDF_FULL = hog_comp_ids.generateHogCompTable(edgeFilePath, hogCompDict, hogCompPath)
  print(' > DONE', flush=True)

  print('3. Append R2 and pvalues from edge file', flush=True)
  print('    > Generate table of [HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B, P_R2, P_Pval, S_R2, S_Pval]', flush=True)
  hogCompNetStatsPath = intermediateFilesPath + '/[HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B, P_R2, P_Pval, S_R2, S_Pval]_FULL_' + argsDict['job_name'] + '.tsv'
  hogCompNetStatsDF_Full = hog_comp_ids.appendStats(geneLookupDF_FULL, edgeFilePath, hogCompNetStatsPath)
  print(' > DONE', flush=True)

  print('4. Drop rows that contain None values, ie there was no COMP_ID match for a given HOG_ID', flush=True)
  hogCompNaPath = intermediateFilesPath + '/[HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B]_DROP_NA_' + argsDict['job_name'] + '.tsv'
  genePairsDropNaPath = intermediateFilesPath + '/[COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval]_DROP_NA_' + argsDict['job_name'] + '.tsv'
  hogCompNetStatsDF_DROP = hog_comp_ids.dropNaRows(hogCompNetStatsDF_Full, hogCompNaPath, genePairsDropNaPath)
  print(' > DONE', flush=True)

  return hogCompNetStatsDF_DROP, genePairsDropNaPath


def formatFullResults(argsDict, intermediateFilesPath, ercFilePath):
  genePairsDropNaPath = intermediateFilesPath + '/[COMP_GENE_A_FULL, COMP_GENE_B_FULL]_DROP_NA_' + argsDict['job_name'] + '.tsv'
  genePairsStatsDropNaPath = intermediateFilesPath + '/[COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval]_DROP_NA_' + argsDict['job_name'] + '.tsv'

  print('1. Get gene pair columns from results file', flush=True)
  df = pandas.read_csv(ercFilePath, sep='\t', usecols=['GeneA_ID', 'GeneB_ID'])
  df = df.rename(columns={'GeneA_ID': 'COMP_GENE_A_FULL', 'GeneB_ID': 'COMP_GENE_B_FULL'})
  print(' > DONE', flush=True)

  print('2. Drop rows that contain none values', flush=True)
  dropNaDF = df.dropna()
  dropNaDF.to_csv(genePairsDropNaPath, sep='\t', index=True)
  print(' > DONE', flush=True)

  print('3. Format gene ids', flush=True)
  indexList= []
  colA = []
  colB = []
  with open(genePairsDropNaPath, 'r') as resultsFile:
    #Skip first line with column titles
    next(resultsFile)
    for line in resultsFile:
      rowData = line.strip().split('\t')
      index = int(rowData[0])
      geneAFormat = hog_comp_ids.formatCompID(rowData[1])
      geneBFormat = hog_comp_ids.formatCompID(rowData[2])
      indexList.append(index)
      colA.append(geneAFormat)
      colB.append(geneBFormat)

  formattedIds = pandas.DataFrame()
  formattedIds['index_col'] = indexList
  formattedIds['COMP_GENE_A'] = colA
  formattedIds['COMP_GENE_B'] = colB
  print(' > DONE', flush=True)

  print('4. Get statistic columns of interest from results file', flush=True)
  stats = pandas.read_csv(ercFilePath, sep='\t', usecols=['P_R2', 'P_Pval', 'S_R2', 'S_Pval', 'Slope'])
  print(' > DONE', flush=True)

  print('5. Merge formatted id columns and statistic columns into one table', flush=True)
  mergeDF = formattedIds.merge(stats, left_on='index_col', right_index=True, how='inner').reset_index()
  mergeDF.drop('index', axis=1, inplace=True)
  mergeDF.drop('index_col', axis=1, inplace=True)
  print(' > DONE', flush=True)

  print('6. Write [COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval, Slope]_DROP_NA table', flush=True)
  mergeDF.to_csv(genePairsStatsDropNaPath, sep='\t', index=False)
  print(' > DONE', flush=True)

  return mergeDF, genePairsStatsDropNaPath
