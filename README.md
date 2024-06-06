# Disk Galaxies in IllustrisTNG

This project aims to extract central disk galaxies from the IllustrisTNG (TNG50-1) Simulation by applying sigma-clipping procedures with reference to the Star Forming Main Sequence and the Mass-Size Relation. The catalogs span across 7 bins in the redshift range $0 \lt z \lt 3$:

$0 \lt z \lt 0.05$<br>
$0.05 \lt z \lt 0.45$<br>
$0.45 \lt z \lt 1.05$<br>
$1.05 \lt z \lt 1.45$<br>
$1.45 \lt z \lt 2.05$<br>
$2.05 \lt z \lt 2.45$<br>
$2.45 \lt z \lt 3.05$<br>

Catalogs carry the following galaxy properties:<br>
<br>
<br>
z:Redshift<br>
Snapshot: Snapshot Number<br>
SubhaloID:<br> 
$M,M_{\odot}:$             Total Subhalo Mass in twice the effective radius(with a fiducial definition of stellar half mass radius),$2R_{eff}$<br> 
$M_{gas},M_{\odot}$:       Gas Mass in $2R_{eff}$<br>
$M_{stellar},M_{\odot}$:   Stellar Mass in $2R_{eff}$<br>
$M_{DM},M_{\odot}$:        Dark Matter Halo Mass in $2R_{eff}$<br>
SFR($M_{\odot}/yr$):       Star Formation Rate in $2R_{eff}$<br> 
$V_{max}$(km/s):           Asymptotic Rotation Curve Velocity<br>
$R_{max}$(kpc):            Distance from Galactic Centre upon achieving $V_{max}$<br>
$\sigma_v$(km/s):          Particle Velocity Dispersion<br>
$R_{1/2M_*}(kpc)$:         Stellar Half Mass Radius

Note that this set of programs requires the `illustris_python` module and simulation data, hence can be used in its present form only within the IllustrisTNG JupyterLab Environment. These files are too large to download thus these scripts most certainly can't be run locally.

### Requirements
`numpy`,`pandas`,`matplotlib`,`astropy`,`argparse`,`pathlib`,`tqdm`,`emcee`,`corner` <br>
Replace [ 29fdb501f084ab3bde756d4827030bcb ] with your own JupyterLab API Key in `get_other_disks.py`,`Generate_Images.py` and `get_disks.py`.


## 1. Extracting disks

To extract based on Speagle+14, run `get_disks.py` from the command line like so:

```
python3 get_disks.py --lower <lower z limit> --upper <upper z>
```
For example, `python3 get_disks.py --lower 0 --upper 0.05` will generate folders named by `lower`, each one containing the following files:


- `0-0.05.csv`: Catalog
- `DISK_ID 0-0.05.csv`: List of SubhaloIDs
- `<Snapshot Number>_compare`: A plot comparing the SFMS determined by different observational papers against TNG Data
- `<Snapshot Number>_fit`: A plot comparing the fit to data after sigma-clipping against Speagle+14 Relation
- `0-0.05MSR`: A plot comparing the data to mass size relation
- `0-0.05MSRFit`: Fit to the data after sigma-clipping w.r.t. Mass Size Relation(not available in high redshift bins because observations disagree with simulations.

Alternatively, to extract based on the best-fitting observational relation in each bin, run `get_other_disks.py`. Appropriate changes need to be made in line 27 of `RandomSample.py`, line 31 of `BTFR_Test.py` and line 31 of `Generate_Images.py` where the output folder names needs to be changed.

## Generate Random Samples

Generates ~100 random samples from the catalogs to run tests on.
```
python3 Random_Sample.py --lower <lower z limit> --upper <upper z>
```
gives the file `<lower>_sample` in the folder `Random Samples Speagle`. 

## Testing with Baryonic Tully Fisher Relation
Test the random samples with the BTFR.
```
python3 BTFR_Test.py --lower <lower z limit> --upper <upper z>
```
Run this to render files in the folder `BTFR_Test_Results/`:

- `<lower>_btfr_data` : A plot of datapoints against the BTFR in McGaugh12
- `<lower>_btfr_fit` : MCMC Fit to datapoints against BTFR in McGaugh12 and Avila-Reese+08
- `<lower>_btfr_corner` : Posterior Distributions of sampled distribution in the above fit.

## Generate Images

Generate images of the samples to check that they are indeed disks. 
```
python3 Generate_Images.py --lower <lower z limit> --upper <upper z>
```
Running this stores PNG Images named `Snap<Snapshot Number>-ID<SubhaloID>` in the folder `Images/<lower>`

The project is a result of a Bachelor's Thesis under Dr. Gauri Sharma and Prof. Shantanu Desai. In case of queries, contact: ep21btech11007@iith.ac.in
