import matplotlib.image as mpimg
import argparse
import pathlib
import os
import requests
import h5py
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument('--lower', type=str, help= "Lower redshift limit")
parser.add_argument('--upper', type=str, help= "Higher redshift limit")
args = parser.parse_args()

if args.lower and args.upper:
    lower = args.lower
    upper = args.upper
        
    print("Testing random sample from redshift bin")
    print("Lower:", args.lower)
    print("Upper:", args.upper)
    
else:
    raise ValueError("Your input is not valid. Please provide the lower and upper limits of the redshift bin you want to test.")
    

df=pd.read_csv(os.path.join("Random_Samples_Speagle/",'{}'.format(lower)+'sample.csv'))
snapshots=df['Snapshot']
ids=df['SubhaloID']

output = "Images/{}".format(lower)
pathlib.Path(output).mkdir(parents=True, exist_ok = True)

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"29fdb501f084ab3bde756d4827030bcb"}

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r




for i in tqdm(range(len(snapshots))):
    
    snap=snapshots[i]
    idd=ids[i]
    
    img=get('https://www.tng-project.org/api/TNG50-1/snapshots/'+str(snap)+'/subhalos/'+str(idd)+'/vis.png?partType=gas&size=50&sizeType=kpc&rasterPx=1100&rotation=face-on&axesUnits=kpc&plotStyle=open&labelZ=True&labelScale=True&labelSim=True')
    
    plt.axis('off')
    plt.imsave(os.path.join(output,'Snap{snap}-ID{idd}.png'.format(snap=snap, idd=idd)),mpimg.imread(img))


