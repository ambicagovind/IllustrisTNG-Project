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
$M,M_{\odot}:$ Total Subhalo Mass in twice the effective radius(with a fiducial definition of stellar half mass radius),$2R_{eff}$<br> 
$M_{gas},M_{\odot}$: Gas Mass in $2R_{eff}$<br>
$M_{stellar},M_{\odot}$: Stellar Mass in $2R_{eff}$<br>
$M_{DM},M_{\odot}$: Dark Matter Halo Mass in $2R_{eff}$<br>
SFR($M_{\odot}/yr$): Star Formation Rate in $2R_{eff}$<br> 
$V_{max}$(km/s): Asymptotic Rotation Curve Velocity<br>
$R_{max}$(kpc): Distance from Galactic Centre upon achieving $V_{max}$<br>
$\sigma_v$(km/s): Particle Velocity Dispersion<br>
$R_{1/2M_*}(kpc)$: Stellar Half Mass Radius

Note that this set of programs requires the `illustris_python` module and simulation data, hence can be used in its present form only within the IllustrisTNG JupyterLab Environment. 

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
- `<lower>-<upper>MSR`: A plot comparing the data to mass size relation
- `<lower>-<upper>MSRFit`: Fit to the data after sigma-clipping w.r.t. Mass Size Relation(not available in high redshift bins because observations disagree with simulations.

\end{itemize}

