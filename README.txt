This is the readme for the models associated with the paper:

Rahul Gupta, Melissa Reneaux, Karmeshu (2016) Role of Macromolecular
Crowding and Geometrical Irregularity at Central Excitatory Synapses
in Shaping Synaptic Transmission. PLOS ONE

MATLAB Scripts for the Numerical Simulation of Glutamate Diffusion,
Stochastic AMPA Receptor Activation and EPSC generation

This code was contributed by Rahul Gupta.

There are three folders:

(1) Glutamate_diffusion.
    It contains scripts for two seperate situations of glutamate
    diffusion viz. (a) Standard diffusion with zero sigma_Dglut and
    (b) Diffusion with non-zero sigma_Dglut. Both these scripts
    perform simulation of the Brownian dynamics of glutamate diffusion
    within the cleft space under the prescribed absorbing and
    reflecting boundary conditions.

(2) Stochastic_AMPAR_activation.
    It contains scripts for the stochastic activation of AMPA
    receptors located in the PSD in response to the temporal profile
    of the local glutamate concentration obtained from the numerical
    simulation of glutamate diffusion. The script
    "AMPAR_activation_GillspiAlgorithm.m" perform the stochastic
    activation of AMPA receptors using modified Gillespie algorithm
    proposed by Skaugen and Walloe (1979). During simulation, it
    repeatedly calls another script "AMPAR_activation_MNscheme.m",
    which contains kinetic parameters of transition among the various
    states of AMPA receptor under the Milstein-Nicoll scheme.

(3) EPSC_generation.
    It contains a single script to compute the resultant profile of
    EPSC generated due to stochastic AMPA receptor activation.
 
