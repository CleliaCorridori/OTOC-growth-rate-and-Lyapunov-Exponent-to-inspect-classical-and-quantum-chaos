In this folder are contained the Jupyter Notebook and the Python scripts to compute the Out-of-Time-Order-Correlator (OTOC)
for the quantum Kicked Rotator and the Lyapunov exponent for the classical Kicked Rotator.

"OTOCfunc.py" 
It contains the functions used to compute the OTOC for the quantum KR for different values of the Kick Strength, K, 
and for different values of h_eff

"CGR_LYAPfunc.py"
It contains the functions to compute:
1. The growth rate of the OTOC for the quantum KR, called CGR;
2. The growth rate of the OTOC for the classical KR;
3. The Lyapunov exponent of the classical KR.

"OTOC.ipynb"
In this notebook is computed the OTOC for different values of the Kick Strength, K, and it is plotted for a selected number of kicks.
Here is plotted also ln(OTOC)/2t for different values of the Kick Strength, K.ù

"OTOCheff.ipynb"
In this notebook is computed the OTOC for different values of h_eff, and it is plotted for a selected number of kicks.
Here is plotted also ln(OTOC)/2t for different values of h_eff.

"CGRandLE.ipynb"
In this notebook, for different values of the Kick Strength, K, in the range [0.01,100] are computed
1. The growth rate of the OTOC for the quantum KR, called CGR;
2. The growth rate of the OTOC for the classical KR;
3. The Lyapunov exponent of the classical KR.
The results are also compared to the analytical behaviour of the Lyapunov exponent for the classical KR.
These quantities are plotted and analized in this notebook.
