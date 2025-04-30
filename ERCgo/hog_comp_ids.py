import pandas
import re

def generateHogCompDict(path):
  df = pandas.read_csv(path, sep='\t', usecols=['HOG_ID', 'Comprehensive_ID'])
  hog_dict = df.set_index('HOG_ID')['Comprehensive_ID'].to_dict()
  return hog_dict

def lookup(hogGene, hogCompDict):
  compGene = hogCompDict.get(hogGene)
  compGene = formatCompID(compGene)
  return compGene


def formatCompID(input):
  if input.startswith('HOG'):
    output = None
  
  else:
    output = re.search(r'(?<=_)(AT)\w+', input).group() 

  return output


def generateHogCompTable(edgeFilePath, hogCompDict, outputPath):
  geneLookupList = []

  with open(edgeFilePath, 'r') as edgeFile:
    #Skip first line with column titles
    next(edgeFile)
    for line in edgeFile:
      convertList = []
      hogPair = line.strip().split('\t')
      convertList.append(hogPair[0])
      convertList.append(hogPair[1])
      compIDA = lookup(hogPair[0], hogCompDict)
      compIDB = lookup(hogPair[1], hogCompDict)
      convertList.append(compIDA)
      convertList.append(compIDB)
      geneLookupList.append(convertList)

  hogCompDF_FULL = pandas.DataFrame(geneLookupList, columns=['HOG_GENE_A', 'HOG_GENE_B', 'COMP_GENE_A', 'COMP_GENE_B'])
  
  print('  > Write [HOG_GENE_A, HOG_GENE_B, COMP_GENE_A, COMP_GENE_B] table to tsv')
  hogCompDF_FULL.to_csv(outputPath, sep='\t', index=False, na_rep='N/A')
  return hogCompDF_FULL

def appendStats(hogCompDF, edgeFilePath, outputPath):
  edgeFileDF = pandas.read_csv(edgeFilePath, sep='\t', usecols=['P_R2', 'P_Pval', 'S_R2', 'S_Pval'])
  df = pandas.concat([hogCompDF, edgeFileDF], axis=1)
  print('    > Write [HOG_ID_A, HOG_ID_B, COMP_ID_A, COMP_ID_B, P_R2, P_Pval, S_R2, S_Pval] table to tsv')
  df.to_csv(outputPath, sep='\t', index=False, na_rep='N/A')
  return df

def dropNaRows(hogCompDF_FULL, hogCompNaPath, compPairsNaPath):
  hogCompDF_DropNa = hogCompDF_FULL.dropna()
  print('   > ' + str(len(hogCompDF_FULL) - len(hogCompDF_DropNa)) + ' gene pairs dropped due to no matching comprehensive id for one or both HOG id(s).')
  print('  > Write [HOG_GENE_A, HOG_GENE_B, COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval] table without NONE values to tsv')
  hogCompDF_DropNa.to_csv(hogCompNaPath, sep='\t', index=False)
  print('  > Write [COMP_GENE_A, COMP_GENE_B, P_R2, P_Pval, S_R2, S_Pval] table without NONE values to tsv')
  hogCompDF_DropNa.to_csv(compPairsNaPath, sep='\t', index=False, columns=['COMP_GENE_A', 'COMP_GENE_B', 'P_R2', 'P_Pval', 'S_R2', 'S_Pval'])
  return hogCompDF_DropNa
