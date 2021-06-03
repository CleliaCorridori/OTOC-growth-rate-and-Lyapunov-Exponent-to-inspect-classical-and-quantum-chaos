# QuantumInformation_project 2020/2021
## Out-of-Time-Ordered Correlator's growth rate and Lyapunov Exponent to inspect classical and quantum chaos
Cristina Cicali, Clelia Corridori, Anna Steffinlongo

In this folder are contained the Jupyter Notebooks and the Python scripts to compute the Out-of-Time-Ordered Correlator's (OTOC)
for the quantum Kicked Rotor and the Lyapunov exponent for the classical Kicked Rotor.

The aim of this project is to study how chaos arises in quantum systems, how it is quantified and how it is related to the definition of classical chaos, inspecting all these questions focusing on the Kicked Rotor, KR. 
In particular, we developed a Python code to numerically simulate the temporal evolution of a Gaussian wave-packet whose dynamics is described by the KR Hamiltonian. This lets us compute the OTOC and study its growth rate, the CGR.
We also implemented a code to numerically reproduce the classical KR, by using the Chirikov standard map. From this, we computed the Lyapunov Exponent, LE, and the growth rate of the classical version of the OTOC.
By comparing these quantities, we will conclude that, since the CGR is more sensitive than the LE to the presence of chaotic islands in the phase space, it could be used to detect local chaos better than global chaos, differently from LE.
