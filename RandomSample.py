import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import pathlib
import os

parser = argparse.ArgumentParser()
parser.add_argument('--lower', type=str, help= "Lower redshift limit")
parser.add_argument('--upper', type=str, help= "Higher redshift limit")
args = parser.parse_args()

if args.lower and args.upper:
    lower = args.lower
    upper = args.upper
        
    print("Extracting random sample from redshift bin")
    print("Lower:", args.lower)
    print("Upper:", args.upper)
    
else:
    raise ValueError("Your input is not valid. Please provide the lower and upper limits of the redshift bin you want to test.")
    

    #Create a folder to save the random samples
    
output = "Random_Samples_Speagle/"
pathlib.Path(output).mkdir(parents=True, exist_ok = True)


df=pd.read_csv(lower+'/'+lower+'-'+upper+'.csv')
ids=df['SubhaloID'].to_numpy()
sample=np.random.choice(ids,size=100)
_,where,_=np.intersect1d(ids,sample,return_indices=True)
df[df.index.isin(where)].to_csv(os.path.join(output,'{}'.format(lower)+'sample.csv'))


