import glob
import sys
import os
import shutil

def findEdgeFile(directory):
  search_result = glob.glob(directory + '/*edges*')
  if len(search_result) == 0:
    sys.exit('> ERCnet edge file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print('> ERCnet edge file found.')
    return search_result[0]

  else:
    sys.exit('> Too many edge files found. Terminating ERCgo.')


def findVerticesFile(directory):
  search_result = glob.glob(directory + '/*vertices*')
  if len(search_result) == 0:
    sys.exit('> ERCnet vertices file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print('> ERCnet vertices file found.')
    return search_result[0]

  else:
    sys.exit('> Too many vertices files found. Terminating ERCgo.')


def findGafFile(directory):
  search_result = glob.glob(directory + '/*.gaf')
  if len(search_result) == 0:
    sys.exit('> GAF file not found. Terminating ERCgo.')

  elif len(search_result) == 1:
    print('> GAF file found.')
    return search_result[0]

  else:
    sys.exit('> Too many GAF files found. Terminating ERCgo.')


def checkOutputDirectory(outPath):
  if not os.path.exists(outPath):
    print('> A directory does not exist at this path.')
    print('> Making directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print('Successful', end='')
  else:
    print('> Existing output directory found.')
    print('> Deleting directory to start fresh.')
    shutil.rmtree(os.path.abspath(outPath))
    print('> Making new directory for writing all output files...', flush=True, end='')
    os.makedirs(outPath)
    print('Successful', end='')
    
  return outPath
