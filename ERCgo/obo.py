from goatools.obo_parser import GODag
import sys

def loadGodAg(oboFilePath):
  try:
    go_dag = GODag(oboFilePath)
  except:
    sys.exit('  > There was an error processing the obo file. Please check this file and try again.')

  return go_dag

def getNames(goIds, goDag):
  goNames = []
  for go_id in goIds:
        if go_id in goDag:
            goNames.append(goDag[go_id].name)
        else:
            goNames.append("GO ID not found")
  return goNames