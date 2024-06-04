import illustris_python as il
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import requests
import argparse
import pathlib
from astropy.cosmology import Planck15 as cosmo
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline

parser = argparse.ArgumentParser()
parser.add_argument('--lower', type=str, help= "Lower redshift limit")
parser.add_argument('--upper', type=str, help= "Higher redshift limit")
args = parser.parse_args()

if args.lower and args.upper:
    lower = args.lower
    upper = args.upper
        
    print("Extracting disks for redshift bin")
    print("Lower:", args.lower)
    print("Upper:", args.upper)
    
else:
    raise ValueError("Your input is not valid. Please provide the lower and upper limits of the redshift bin.")

headers={"api-key":"29fdb501f084ab3bde756d4827030bcb"}
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


#Load snapshot numbers and corresponding redshifts
baseUrl = 'http://www.tng-project.org/api/TNG50-1/snapshots/'
snaps=get(baseUrl)
snaprs=np.array([(int(snaps[i]['number']),snaps[i]['redshift']) for i in range(len(snaps))])

#Extract snapshots within the required range
redshifts=np.array([i[1] for i in snaprs])
index=np.where((redshifts>float(lower))& (redshifts<float(upper)))
snaps=snaprs[index[0]]


print("The snapshots in the required redshift bin are:")

#Create Output Folder
output = lower+"_other/"
pathlib.Path(output).mkdir(parents=True, exist_ok = True)

h=0.6774

def assign_fields(data,redshift,snap_no):
    
    '''
    Segregate snapshot data into lists of relevant data
    '''
    print('Preprocessing data...')
    a=1/(1+redshift)
    
    #Stellar mass
    ms=data['SubhaloMassInRadType'][:,4]*1e10/h
    
    #Effective radius
    r=data['SubhaloHalfmassRadType'][:,4]*a/h
    
    #Star formation rate
    sfr=data['SubhaloSFRinRad']
    
    #Redshifts
    z=redshift*np.ones(len(ms))
    
    #Snapshot No
    s=snap_no*np.ones(len(ms))
    
    flags=data['SubhaloFlag']
    ids=np.array(list(range(len(ms))))
    
    #Retain only galaxies of cosmological origin
    good_ones=np.where((flags==True) &(ms > 0) & (sfr > 0))
    fields=[ms,r,sfr,z,flags,ids,s]
    newfields=[field[good_ones] for field in fields]
    ms,r,sfr,z,flags,ids,s=newfields
    
    log_ms=np.log10(ms)
    log_sfr=np.log10(sfr)
    
    #Retain galaxies within the require range
    req_range=np.where((log_ms>8)&(log_ms<11.5))

    fields=[ms,r,sfr,z,flags,ids,s]
    newfields=[field[req_range] for field in fields]
    ms,r,sfr,z,flags,ids,s=newfields
    
    return ms,r,sfr,z,flags,ids,s

def speagle(logmass,z):
    t=cosmo.age(z).value
    slope=0.84-0.026*t
    intercept=6.51-0.11*t
    return slope*(logmass)-intercept

def renzini(logmass,z):
    #z<0.085
    return 0.76*logmass-7.64 

def pearson(logmass,z):
    #0.2<z<6
    slope=1.1+0.5*np.log(0.03+z)
    intercept=0.38+0.12*z-10.5*slope
    return slope*logmass+intercept 

def kurczynski(logmass,z):
    if z<=1:
        return 0.919*logmass-8.394
    elif z<=1.5:
        return 0.825*logmass-7.474
    elif z<=2:
        return 0.867*logmass-7.484
    elif z<=2.5:
        return 0.849*logmass-7.513
    else:
        return 0.899*logmass-7.729


def bins(log_ms,log_sfr,nbins=13):
    
    '''
    13-Bin log_ms and assign corresponding datapoints, calculate <logSFR> in each bin
    '''
    mass_bins=np.linspace(log_ms.min(),log_ms.max(),nbins,endpoint=True)
    bin_indices=np.digitize(log_ms, mass_bins, right=False)
    binned_data={bin_start:[] for bin_start in mass_bins[0:-1]}

    for i, bin_start in enumerate(mass_bins[0:]):
        mask=bin_indices==i+1
        binned_data[bin_start]=list(zip(log_ms[mask],log_sfr[mask]))
   
    x=[(mass_bins[i]+mass_bins[i+1])/2 for i in range(len(mass_bins)-1)]
    y = [np.mean([i[1] for i in binned_data[bin_start]]) for bin_start in mass_bins[0:-1]]
    dy = [np.std([i[1] for i in binned_data[bin_start]])/len(binned_data[bin_start]) for bin_start
          in mass_bins[0:-1]]
    return x,y,dy

def line(x,m,c):
    return m*x+c

def fit(func,scat,label,z,log_ms_binned,log_sfr_binned, std_logsfr,log_ms,log_sfr, snap_no):
    
    print('Fitting..')
    
    '''
    Fit the binned(filtered) data to a straight line and compare 
    '''
    
    x = np.linspace(8, 11.5, len(log_ms))
    
    cs1=CubicSpline(x,func(x,z)+scat)
    cs2=CubicSpline(x,func(x,z)-scat)

    try:
        
        popt,pcov=curve_fit(line,log_ms_binned,log_sfr_binned,sigma=std_logsfr)
        p_err=np.sqrt(np.diag(pcov))
        print('Optimal m,c to binned data:',popt)
        print('Errors:',p_err)
        
    except:
        
        print('Unable to fit alongwith errors, fitting without errors...')
        popt,pcov=curve_fit(line,log_ms_binned,log_sfr_binned)
        p_err=np.sqrt(np.diag(pcov))
        print('Optimal m,c to binned data:',popt)
        print('Errors:',p_err)

    plt.figure(figsize=(5,4))
    plt.scatter(log_ms,log_sfr, s=3, label='Data')
    
    plt.plot(x,func(x,z), label=label, linewidth=2, color='black')
    plt.plot(x,line(x,popt[0],popt[1]), label='Best Fit', linewidth=2, color='fuchsia')
    plt.errorbar(log_ms_binned, log_sfr_binned, std_logsfr,fmt='o',capsize=5,color='orange',capthick=2,label='Binned Data')
    plt.fill_between(x,cs1(x),cs2(x),alpha=0.5,color='green')

    plt.xlabel('log Galaxy Stellar Mass [$M_\\odot$, $r < 2R_{eff}$]')
    plt.ylabel('log Star Formation Rate [yr$^{-1}$,$r < 2R_{eff}$]')
    
    plt.legend(frameon=False)
    plt.savefig(os.path.join(output,'{}'.format(snap_no)+'_MSfit.png'))
    print('{}'.format(snap_no)+'MSfit saved.')
    
def extract(func,scat,label,rs,log_ms,log_sfr,fields,snap_no):
    
    good_ones=np.where((abs(func(log_ms,rs)-log_sfr)<2*scat))
    newfields=[field[good_ones] for field in fields]
    ms2,r2,sfr2,z2,flags2,ids2,snaps12=newfields
    
    log_ms=np.log10(ms2)
    log_sfr=np.log10(sfr2)
    log_ms_binned,log_sfr_binned,std_logsfr=bins(log_ms,log_sfr)
    
    fit(func,scat,label,rs,log_ms_binned,log_sfr_binned,std_logsfr, log_ms,log_sfr, snap_no)
    
    return ms2,r2,sfr2,z2,flags2,ids2,snaps12

def get_func_deets(rs):
    function_map = {
      0.05: (renzini, 0.3, "Renzini&Peng15"),
      0.45: (speagle, 0.3, "Speagle+14"),
      1.05: (pearson, 0.3, "Pearson+18"),
      1.5: (kurczynski, 0.383, "Kucrzynski+16"),
      2: (kurczynski, 0.354, "Kucrzynski+16"),
      2.5: (kurczynski, 0.399, "Kucrzynski+16")
  }
    
    default_func = "kurczynski"
    default_scat = 0.369
    default_label = "Kucrzynski+16"

    for threshold, (func, scat, label) in function_map.items():
        if rs < threshold:
            return func, scat, label

    return default_func, default_scat, default_label

#Empty arrays to store required data
all_ms=np.array([])
all_r=np.array([])
all_sfr=np.array([])
all_z=np.array([])
all_id=np.array([])
all_mg=np.array([])
all_mdm=np.array([])
all_mt=np.array([])
all_v=np.array([])
all_s=np.array([])
all_rmax=np.array([])
all_snaps=np.array([])


basePath = 'sims.TNG/TNG50-1/output'
field=['SubhaloHalfmassRadType','SubhaloFlag','SubhaloSFRinRad','SubhaloMassInRadType']
   
for snap_no,rs in snaps[::-1]:
    
    print('Processing for Snapshot# {}...'.format(int(snap_no)))
    print('Loading TNG Data...')
    data=il.groupcat.loadSubhalos(basePath,snap_no,fields=field)
    
    ms,r,sfr,z,flags,ids,snaps1=assign_fields(data,rs,snap_no)
    log_ms=np.log10(ms)
    log_sfr=np.log10(sfr)
    
    rs=np.round(z[0],3)
    
    fields=[ms,r,sfr,z,flags,ids,snaps1]
    func,scat,label=get_func_deets(rs)        
    ms,r,sfr,z,flags,ids,snaps1=extract(func,scat,label,rs,log_ms,log_sfr,fields,snap_no)

    all_ms=np.concatenate((all_ms,ms))
    all_r=np.concatenate((all_r,r))
    all_sfr=np.concatenate((all_sfr,sfr))
    all_z=np.concatenate((all_z,z))
    all_id=np.concatenate((all_id,ids))
    all_snaps=np.concatenate((all_snaps,snaps1))
    
    
#Extraction based on mass size relation for locals

if float(lower)==0:
    
    print('Working on Mass Size Relation...')
    log_ms=np.log10(all_ms)
    log_r=np.log10(all_r)

    x=np.array([9.00507614213198, 9.302030456852792, 9.598984771573605, 
       9.898477157360405, 10.200507614213198, 10.50253807106599, 
       10.802030456852792, 11.101522842639593, 11.403553299492385])
    
    y=np.array([0.4129205921938089,0.4979811574697174,0.570121130551817, 
       0.619650067294751,0.682099596231494,0.7574697173620457,
       0.8242261103633917,0.9189771197846568,1.0201884253028264])

    z=np.poly1d(np.polyfit(x,y,4))

    good_ones=np.where(abs(z(log_ms)-log_r)<=0.1)
    fields=[all_ms,all_r,all_sfr,all_z,all_id,all_snaps]
    newfields=[field[good_ones] for field in fields]
    all_ms,all_r,all_sfr,all_z,all_id,all_snaps=newfields

    log_ms=np.log10(all_ms)
    log_r=np.log10(all_r)
    
    logms_to_plot,logr_to_plot,std=bins(log_ms,log_r)
    poly_bins=np.poly1d(np.polyfit(logms_to_plot,logr_to_plot,4))
    
    fig=plt.figure(figsize=(5,4))
    ax=fig.add_subplot(111)
    x_values=np.linspace(log_ms.min(),log_ms.max(),100)
    ax.plot(log_ms,log_r,'.',markersize=5)
    ax.plot(x_values,z(x_values),label='Lapi et al')
    ax.fill_between(x_values,z(x_values)+0.05,z(x_values)-0.05,color='pink',alpha=0.5,zorder=2)
    ax.plot(x_values,poly_bins(x_values),color='red',zorder=3,label='Polyfit to binned data')
    ax.errorbar(logms_to_plot,logr_to_plot,yerr=std,color='red',zorder=3,fmt='o',markersize=2,
            label='Mean of binned data')
    ax.set_xlabel('log [Galaxy Stellar Mass [$M_\odot$, $r < 2*R_{eff}]$]')
    ax.set_ylabel('log  $R_{eff}$')
    ax.legend()
    plt.savefig(os.path.join(output,'{}'.format(snap_no)+'_Rfit.png'))
    print('{}'.format(snap_no)+'Rfit saved.')
    
else:
        
    log_ms=np.log10(all_ms)
    log_r=np.log10(all_r) 
    log_ms_scaled=log_ms-np.log10(5e+10)
    x=np.linspace(min(log_ms_scaled),max(log_ms_scaled),100)

    def mass_size(lower):
        if lower==0.05:
            return [0.25,0.86,0.20]
        elif lower==0.45:
            return [0.22,0.78,0.21]
        elif lower==1.05:
            return [0.22,0.7,0.21]
        elif lower==1.45:
            return [0.23,0.65,0.24]
        elif lower==2.05:
            return [0.22,0.55,0.23]
        elif lower==2.45:
            return [0.18,0.51,0.24]

    m1,c1,scat=mass_size(float(lower))
    
    fig=plt.figure(figsize=(5,4))
    ax=fig.add_subplot(111)
    ax.plot(log_ms_scaled,log_r,'.',markersize=2)
    ax.set_xlabel('log [Galaxy Stellar Mass [$M_\odot$, $r < R_{eff}]/5*10^{10}M_\odot$]')
    ax.set_ylabel('log  $R_{eff}$')    
    ax.plot(x,m1*x+c1,label='late-type',linewidth=4)
    ax.plot(x,m1*x+c1-scat,'--',color='black',linewidth=2)
    ax.plot(x,m1*x+c1+scat,'--',color='black',linewidth=2)
    ax.set_xlim(-2.7,0.802)
    ax.legend(frameon=False)
    plt.savefig(os.path.join(output,'{lower}-{upper}MSR.png'.format(lower=lower,upper=upper)))
    
    if float(lower)<1.45:
        good_ones=np.where(abs(m1*(log_ms_scaled)+c1-log_r)<2*scat)
        fields=[all_ms,all_r,all_sfr,all_z,all_id,all_snaps]
        newfields=[field[good_ones] for field in fields]
        all_ms,all_r,all_sfr,all_z,all_id,all_snaps=newfields

        log_ms=np.log10(all_ms)
        log_r=np.log10(all_r)
        log_ms_scaled=log_ms-np.log10(5e+10)
        logms_scaled_binned,logreff_binned,std=bins(log_ms_scaled,log_r)

        try:
            m2,c2=curve_fit(line,logms_scaled_binned,logreff_binned,sigma=std)[0]
    
        except:
            logms_scaled_binned,logreff_binned,std=bins(log_ms_scaled,log_r,10)
            m2,c2=curve_fit(line,logms_scaled_binned,logreff_binned,sigma=std)[0]
    
        fig=plt.figure(figsize=(5,4))
        ax=fig.add_subplot(111)
        ax.plot(log_ms_scaled,log_r,'.',markersize=2)
        ax.plot(x,m1*x+c1,label='van der Wel et al')
        ax.plot(x,x*m2+c2,label='best fit',zorder=3,color='pink',linewidth=3)
        ax.errorbar(logms_scaled_binned,logreff_binned,yerr=std,color='red',zorder=3,label='Mean of binned data',fmt='o')
        ax.fill_between(x,m1*x+c1-scat,m1*x+c1+scat,color='yellow',alpha=0.5,zorder=2)
        ax.set_xlabel('log [Galaxy Stellar Mass [$M_\odot$, $r < 2*R_{eff}]/5*10^{10}M_\odot$]')
        ax.set_ylabel('log  $R_{eff}$')
        ax.set_ylim(-1,1.5)
        ax.legend(frameon=False)
        plt.savefig(os.path.join(output,'{lower}-{upper}MSRFit.png'.format(lower=lower,upper=upper)))
    
df=pd.DataFrame({'Snapshot':all_snaps.astype(int),'SubhaloID':all_id.astype(int),
                 '$z$':np.round(all_z,3)})
df.to_csv(os.path.join(output,'DISK_ID:{lower}-{upper}.csv'.format(lower=lower, upper=upper)))

def dump_to_csv(snap):
    
    rs = np.round(snap[1], 3)
    sub_ids = df[df['$z$'] == rs]['SubhaloID'].to_numpy().astype(int)
    snapnos = df[df['$z$'] == rs]['Snapshot'].to_numpy().astype(int)
    a = 1 / (1 + rs)

    full_data = il.groupcat.loadSubhalos(basePath, snap[0], 
                                         fields=['SubhaloMassInRad', 
                                                 'SubhaloMassInRadType',
                                                 'SubhaloHalfmassRadType', 
                                                 'SubhaloSFRinRad', 
                                                 'SubhaloVmax', 
                                                 'SubhaloVmaxRad',
                                                 'SubhaloVelDisp'])
    
    print('Data for snapshot#{} Loaded'.format(snap[0]))

    mg = full_data['SubhaloMassInRadType'][sub_ids][:, 0] * 1e10 / h
    ms = full_data['SubhaloMassInRadType'][sub_ids][:, 4] * 1e10 / h
    mdm = full_data['SubhaloMassInRadType'][sub_ids][:, 1] * 1e10 / h
    shmr = full_data['SubhaloHalfmassRadType'][sub_ids][:, 4] * a / h
    sfr = full_data['SubhaloSFRinRad'][sub_ids]
    v = full_data['SubhaloVmax'][sub_ids]
    s = full_data['SubhaloVelDisp'][sub_ids]
    tm = full_data['SubhaloMassInRad'][sub_ids] * 1e10 / h
    vmr= full_data['SubhaloVmaxRad'][sub_ids] * a/ h
    rs = np.full((len(mg)), rs)

    df2 = pd.DataFrame({'z': rs, 
                        'Snapshot': snapnos,
                        'SubhaloID': sub_ids, 
                        'M,M_{\odot}': tm, 
                        '$M_{gas},M_{\odot}$': mg, 
                        '$M_{stellar},M_{\odot}$': ms, 
                        '$M_{DM},M_{\odot}$': mdm, 
                        'SFR($M_{\odot}/yr$)': sfr, 
                        '$V_{max}$(km/s)': v, 
                        '$R_{max}$(kpc)': vmr,
                        '$\sigma_v$(km/s)': s,
                        '$R_{1/2M_*}(kpc)$': shmr})

    csv_file_path = '{lower}-{upper}.csv'.format(lower=lower, upper=upper)

    if not os.path.isfile(os.path.join(output,csv_file_path)):
        df2.to_csv(os.path.join(output,csv_file_path), mode='w', index=False, header=True)
    else:
        df2.to_csv(os.path.join(output,csv_file_path), mode='a', index=False, header=False)

    print('Catalog for snapshot#{} Saved'.format(snap[0]))
    
for snap in snaps[::-1]:
    dump_to_csv(snap)