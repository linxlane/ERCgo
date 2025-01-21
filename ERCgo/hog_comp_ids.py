import glob
import pandas

def findEdgeFile(directory):
  search_result = glob.glob(directory + '/*edges*')
  if len(search_result) == 0:
    print('ERCnet edge file not found.')

  elif len(search_result) == 1:
    print('ERCnet edge file found.')
    return search_result[0]

  else:
    print('Too many edge files found.')

def findVerticesFile(directory):
  search_result = glob.glob(directory + '/*vertices*')
  if len(search_result) == 0:
    print('ERCnet vertices file not found.')

  elif len(search_result) == 1:
    print('ERCnet vertices file found.')
    return search_result[0]

  else:
    print('Too many vertices files found.')

def generateLookupDict(path):
  df = pandas.read_csv(path, sep='\t', usecols=['HOG_ID', 'Comprehensive_ID'])
  hog_dict = df.set_index('HOG_ID')['Comprehensive_ID'].to_dict()
  return hog_dict

def lookup(hogGene, hogCompDict):
  compGene = hogCompDict.get(hogGene)
  compGene = formatCompID(compGene)
  return compGene

def formatCompID(compGene):
  if compGene.startswith('HOG'):
    compGene = None
  
  else:
    compGene = compGene[6:]
    if len(compGene) > 10:
      compGene = compGene[:9]

  return compGene

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

  geneLookupDF_FULL = pandas.DataFrame(geneLookupList, columns=['HOG_GENE_A', 'HOG_GENE_B', 'COMP_GENE_A', 'COMP_GENE_B'])
  geneLookupDF_FULL.to_csv(outputPath + '/HOG_COMP_TABLE_FULL.tsv', sep='\t', index=False, na_rep='N/A')
  return geneLookupDF_FULL

def dropNaRows(geneLookupDF_FULL, outputPath):
  geneLookupDF_DROP = geneLookupDF_FULL.dropna()
  print(str(len(geneLookupDF_FULL) - len(geneLookupDF_DROP)) + ' gene pairs dropped due to no matching comprehensive id for one or both HOG id(s).')
  geneLookupDF_DROP.to_csv(outputPath + '/HOG_COMP_TABLE_DROP_NA.tsv', sep='\t', index=False)
  return geneLookupDF_DROP