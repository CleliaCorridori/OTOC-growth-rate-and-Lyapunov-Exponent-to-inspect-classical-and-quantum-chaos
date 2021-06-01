# QuantumInformation_project 2020/2021
## Out-of-Time-Order-Correlator growth rate and Lyapunov Exponent to inspect classical and quantum chaos
Cristina Cicali, Clelia Corridori, Anna Steffinlongo

In this folder are contained the Jupyter Notebook and the Python scripts to compute the Out-of-Time-Order-Correlator (OTOC)
for the quantum Kicked Rotator and the Lyapunov exponent for the classical Kicked Rotator.

The aim of this project is to study how chaos arises in quantum systems, how it is quantified and how it is related to the definition of classical chaos, inspecting all these questions focusing on the Kicked Rotor, KR. 
In particular, we developed a Python code to numerically simulate the temporal evolution of a gaussian wave-packet whose dynamics is described by the KR Hamiltonian. This lets us compute the OTOC and study its growth rate, CGR.
We also implemented a code to numerically reproduce the classical KR, by using the Chirikov standard map. From this, we computed the Lyapunov Exponent, LE, and the growth rate of the classical version of the OTOC. \\
By comparing these quantities, we will conclude that, since the CGR is more sensitive than the LE to the presence of chaotic islands in the phase space, it could be used to detect local chaos better than global chaos.

Used scripts:
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
Here is plotted also ln(OTOC)/2t for different values of the Kick Strength, K.

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
