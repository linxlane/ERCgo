from goatools.anno.gaf_reader import GafReader
import pandas

def generateGeneGoAssocDF(gafPath):
  ogaf = GafReader(gafPath)

  ns2assc = ogaf.get_ns2assc()

  geneAndGoList = []

  for namespace, associations in ns2assc.items():
    for protein_id, go_ids in sorted(associations.items()):
      geneGoPair = []
      geneGoPair.append(protein_id)
      geneGoPair.append(sorted(go_ids))
      geneAndGoList.append(geneGoPair)

  geneAndGoListDf = pandas.DataFrame(geneAndGoList, columns=['Gene_ID', 'GO_Terms'])
  return geneAndGoListDf
