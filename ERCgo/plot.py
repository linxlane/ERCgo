import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

def collectAndPrep(sharedGOFilesList, writePath):
  # Create an empty list to store the dataframes
  df_list = []

  # Loop through the files and read each one into a dataframe
  for file in sharedGOFilesList:
    df = pandas.read_csv(file, sep='\t')
    df_list.append(df)

  # Concatenate all the dataframes into a single dataframe
  final_df = pandas.concat(df_list, ignore_index=True)

  final_df.to_csv(writePath, sep='\t', index=False,)
  
  return final_df

def dropZeros(df):
  return df[df.Number_of_Shared_GO != 0]

def dropOnes(df):
  return df[df.Number_of_Shared_GO != 1]

def seabornKDE(sharedGoData, writePath):
  sns.kdeplot(data=sharedGoData, x='Number_of_Shared_GO', hue='label')
  mpl.rcParams['pdf.fonttype'] = 42
  plt.savefig(writePath, format = 'pdf', transparent = True)
  plt.show()
