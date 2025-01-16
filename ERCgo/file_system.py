import os
import glob

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