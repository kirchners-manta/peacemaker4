What is Peacemaker?
###########
Peacemaker uses the laws of statistical thermodynamics to calculate the thermodynamic properties 
of pure liquids and liquid mixtures.
It is based on the **Q**uantum **C**luster **E**quilibrium (QCE) theory, which is the idea, that 
the liquid bulk system can be described as a dense distribution of statistically reoccurring
molecular cluster motifs. 

What you need to provide in order to use Peacemaker
***********
- A set of clusters, which are representative for the system you want to investigate
- The vibrational frequencies of the clusters
- The volumes for the monomers
- The molar amounts of the components in your system
- The adiabatic interaction energy of the clusters 
.. math::

    \Delta_{bild}\epsilon(t_i w_j) = \epsilon(t_i w_j) - i\epsilon(t_1) - j\epsilon(w_1)

What you get from Peacemaker
***********
- The Gibbs free energy of the system
- The Helmholtz free energy of the system
- The internal energy of the system
- The enthalpy of the system
- The entropy of the system
- The heat capacity at constant volume 
- The heat capacity at constant pressure
- The volume of the system
- The population each cluster motif
- The concentration of each cluster motif
- A set of contributions for each degree of freedom to the thermodynamic quantities
- A set of contributions for each cluster to the thermodynamic quantities

.. figure:: figures/pm.png
    :width: 400
    :align: center

    Peacemaker