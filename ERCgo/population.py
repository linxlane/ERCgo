import pandas
from collections import Counter

def populationUniqueGenes(genePairsDf):
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


def calculatePopulationFrequencies(genePairs, geneGoDict, outputPath, jobName):
  ##Summarize GO term population for all GO terms associated with genes in the current edge file
  print('2. Gather GO term population of all genes in edge file')
  #Get unique genes from column A and B, returns set of gene comp ids
  print('  > Get unique genes from column A and B')
  uniqueIds = populationUniqueGenes(genePairs)

  #Use counter container to get frequencies of each go term in the population
  print('  > Create counter for GO term frequencies')
  goTermFreq = countGoTermsFrequency(uniqueIds, geneGoDict)
  print('  > Write frequency counter to file')
  populationWritePath = outputPath + '/Population_Frequencies_' + jobName + '.tsv'
  writePopulationFrequencies(populationWritePath, goTermFreq)

  return goTermFreq
