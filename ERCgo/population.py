import pandas
from collections import Counter

def populationUniqueGenes(compIDPairsPath):
  genePairsDf = pandas.read_csv(compIDPairsPath, sep='\t')
  uniqueGenes = set(list(genePairsDf['COMP_GENE_A']) + list(genePairsDf['COMP_GENE_B']))
  return uniqueGenes

def countGoTermsFrequency(uniqueIds, assoc_dict):
  #Collect GO terms for each unique gene in the population
  goTerms = []
  keyErrorCounter = 0
  for gene in uniqueIds:
    if gene in assoc_dict.keys():
      goTerms.extend(assoc_dict[gene])
    else:
      keyErrorCounter += 1
  print(f'   > Number of genes in edge file not found in GAF: {keyErrorCounter}')

  #Create counter collection for GO term frequency
  populationCounts = Counter(goTerms)
  return populationCounts

def writePopulationFrequencies(populationWritePath, frequencyCounts):
  with open(populationWritePath, 'w') as file:
    for key, value in frequencyCounts.items():
      file.write(f'{key}\t{value}\n')
