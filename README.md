# MCMC-17O-plant-water

# Bayesian Triple Oxygen Isotope Parameter Estimation for plant waters

## Overview
This repository contains the MATLAB script **Bayesian_Equisetum.m**, which implements a Bayesian Markov Chain Monte Carlo (MCMC) approach for estimating parameters related to the triple oxygen isotope composition of plant water in _Equisetum_. The methodology is based on Bayesian inference to model the isotopic fractionation processes affecting water compositio.

The script is associated with the paper:

**Sharp, Z., et al. (2025). Extreme triple oxygen isotope fractionation in Equisetum. *PNAS, in review*.**

## Features
- Implements Bayesian parameter estimation using the MCMC toolbox for MATLAB.
- Incorporates external functions for plant water isotopes.
- Supports analysis of measured data.
- Allows estimation of key environmental parameters (T, RH, d18O of soil water), as a check on model applicability.

## Dependencies
This script requires the following MATLAB toolboxes and external functions:

### Required MATLAB Toolboxes:
- **MCMC toolbox** by Marko Laine: [https://mjlaine.github.io/mcmcstat/](https://mjlaine.github.io/mcmcstat/)

### Recommended MATLAB Toolboxes:
- **DERIVESTsuite** for numerical differentiation.

### External Functions (included in the repository, must be placed by the user in their working directory):
- `Bayesian_Equisetum.m`
- `ssfun.m`
- `model_equisetum3.m`
- `model_equisetum3_plot.m`
- `prediction_plot_from_n.m`

(Future versions of this will simplify this code. Stay tuned!)

## Input Data
The script requires a dataset in MATLAB `.mat` format with the following structure:
- `data.ydata`: Measured **δ¹⁸O, Δ'¹⁷O** values (permil/per meg)
- `data.xdata`: an index vector
- `data.yerr`: Uncertainty in **δ¹⁸O, Δ'¹⁷O** (1σ standard deviation)
- `data.fractional_positions`: fractional positions along the plant length

Example dataset: `equisetum.mat`

## Parameter Definitions
The model parameters estimated using MCMC are:
(VJH TO WRITE)

## Running the Script
1. Load MATLAB and navigate to the script directory.
2. Ensure that the required `.mat` dataset is available in the working directory.
3. Run the script by executing:
   ```matlab
   Bayesian_Equisetum
4. The script will perform MCMC sampling and generate posterior distributions of the model parameters.

## Output
The script produces:
- **MCMC chains**: Stored in `chain` variable.
- **Plots**:
  - Trace plots of MCMC chains.
  - Posterior density estimates.
  - Parameter correlation scatter plots.

 
## Citation
If you use this script in your research, please cite:

Sharp, Z., et al. (2025). Extreme triple oxygen isotope fractionation in Equisetum. *PNAS, in review*.

## License
This code is provided under the MIT License. See `LICENSE` for details.

---

For any questions, please contact: **vincent.john.hare@gmail.com**

