import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import pathlib
import os
from astropy.constants import G
from astropy import units as u
from scipy.optimize import curve_fit
import emcee
import corner

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
output = "BTFR_Test/"
pathlib.Path(output).mkdir(parents=True, exist_ok = True)


def bin_data(x,y,nbins=6):
    
    x_bins=np.linspace(x.min(),x.max(),nbins,endpoint=True)
    bin_indices=np.digitize(x, x_bins, right=False)
    binned_data={bin_start:[] for bin_start in x_bins[0:-1]}

    for i, bin_start in enumerate(x_bins[0:]):
        mask=bin_indices==i+1
        binned_data[bin_start]=list(zip(x[mask],y[mask]))
   
    x1=np.array([(x_bins[i]+x_bins[i+1])/2 for i in range(len(x_bins)-1)])
    y1 = np.array([np.mean([i[1] for i in binned_data[bin_start]]) for bin_start in x_bins[0:-1]])
    dy = np.array([np.std([i[1] for i in binned_data[bin_start]])/len(binned_data[bin_start]) for bin_start
          in x_bins[0:-1]])
    return x1,y1,dy

m_baryonic=(df['$M_{stellar},M_{\odot}$']+df['$M_{gas},M_{\odot}$']).to_numpy()
G = G.to(u.kpc / u.Msun * (u.km/u.s)**2).value
v_circ=np.sqrt((G*df['M,M_{\odot}']/2/df['$R_{1/2M_*}(kpc)$']).to_numpy())

xdata=np.linspace(30,400,100)
ydata=47*np.power(xdata,4)


fig,ax=plt.subplots(figsize=(6,5))
plt.scatter(np.log10(v_circ),np.log10(m_baryonic),s=4,color='blue',label='${lower}<z<{upper}$ Disks from TNG50'.format(lower=lower, upper=upper),alpha=0.5)
plt.plot(np.log10(xdata),np.log10(ydata),label='McGaugh+12',linewidth=1,color='red')
plt.fill_between(np.log10(xdata),np.log10(53*np.power(xdata,4)),np.log10(41*np.power(xdata,4)),alpha=0.3,color='red')
plt.xlabel('$log V_{circ}[km/s]$')
plt.ylabel('$M_b,M_{\odot}$')
plt.legend(frameon=False)
plt.savefig(os.path.join(output,'{}_btfr_data.png'.format(lower)))

x,y,std=bin_data(np.log10(v_circ),np.log10(m_baryonic))

def line(theta,x1):
    m,c=theta
    return m*x1+c

def log_prior(theta): 
    
    return 0.0

def log_likelihood(theta, x1,y1,dy1):
    
    yM=line(theta,x1)    
    return -0.5 * (np.sum(2*np.pi*dy1**2)+np.sum(((y1 - yM)/dy1) ** 2))

def log_posterior(theta, x1,y1,dy1):
    
    lp=log_prior(theta)
    
    if not np.isfinite(lp):
        return -np.inf
    
    return lp + log_likelihood(theta,x1,y1,dy1)

nwalkers=50
nsteps=2000
p1=np.random.uniform(0,10,50)
p2=np.random.uniform(-10,10,50)
p0=np.vstack((p1,p2)).T
data=(x,y,std)

def main(p0,nwalkers,nsteps,ndim,log_posterior,data):
    
    sampler=emcee.EnsembleSampler(nwalkers,ndim,log_posterior,args=data)
    print('Running burn-in...')
    pos,_,_=sampler.run_mcmc(p0,1000,progress=True)
    sampler=sampler.reset()
   
    sampler=emcee.EnsembleSampler(nwalkers,ndim,log_posterior,args=data)
    print('Running production...')
    pos,prob,state=sampler.run_mcmc(pos,nsteps,progress=True)
    return sampler,pos,prob,state

sampler,pos,prob,state = main(p0,nwalkers,nsteps,2,log_posterior,data)
samples=sampler.flatchain
theta_max=samples[np.argmax(sampler.flatlnprobability)]
best_fit_model=line(theta_max,x)

labels=['b','logA']
fig = corner.corner(samples,labels=labels,plot_datapoints=True,
     show_titles=True,use_math_text=True,truths=theta_max,
    color='black',truth_color='red',quantiles=[0.68,0.95,0.99])
plt.savefig(os.path.join(output,'{}_btfr_corner.png'.format(lower)))


def ar(x):
    return 10**((np.log10(x)+1.3)/0.333)

def UPar(x):
    return 10**((np.log10(x)+1.3+0.063)/0.333)

def DOWNar(x):
    return 10**((np.log10(x)+1.3-0.063)/0.333)


ar_y=ar(xdata)

plt.figure()

plt.plot(xdata,ar_y,label='Avila-Reese+2008',linewidth=3)
plt.scatter(v_circ,m_baryonic,s=4,color='blue',label='${lower}<z<{upper}$ Disks from TNG50'.format(lower=lower, upper=upper),alpha=0.5)
plt.plot(xdata,ydata,label='McGaugh+12',linewidth=1,color='red')
plt.fill_between(xdata,53*np.power(xdata,4),41*np.power(xdata,4),alpha=0.3,color='red')
plt.fill_between(xdata,UPar(xdata),DOWNar(xdata),alpha=0.2,color='green')
plt.errorbar(10**x,10**y,10**std,label='Mean of Binned Data',fmt='o',elinewidth=3,color='green',capsize=4)
plt.plot(xdata,xdata**theta_max[0]*10**theta_max[1],label='Fit to Binned Data',linewidth=2,color='black')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$V_{circ}[km/s]$',fontsize=14,fontweight='bold')
plt.ylabel('$M_{baryonic},M_{\odot}$',fontsize=14,fontweight='bold')
ax.spines['top'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)
plt.xlim(30,400)
plt.legend(frameon=False)
plt.savefig(os.path.join(output,'{}_btfr_fit.png'.format(lower)))
