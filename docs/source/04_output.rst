Peacemaker writes results for the best :math:`a_{mf}`, :math:`b_{xv}` parameter pair if it could be 
determined or for the first :math:`a_{mf}`, :math:`b_{xv}` pair otherwise.
In the following, the output files are briefly described.
All output files contain the temperature in :math:`K` in column 1.

.. line-block::

    **volume.dat**
        Contains the volume and related quantities: volume :math:`V` in :math:`dm^3`, exclusion volume :math:`V_\mathrm{excl}`` in :math:`dm^3`, volumetric expansion coefficient :math:`\alpha` in :math:`K^{-1}`.

    **thermo0.dat**
        Contains thermodynamic quantities that do not depend on any derivative: Helmholtz free energy :math:`A` in :math:`kJ`, Gibbs free energy :math:`G` in :math:`kJ`.
    
    **thermo1.dat**
        Contains thermodynamic quantities that depend on first derivatives: internal energy :math:`U` in :math:`kJ`, enthalpy :math:`H` in :math:`kJ`, entropy :math:`S` in :math:`kJ`.
        
    **thermo2.dat**
        Contains thermodynamic quantities that depend on second derivatives: heat capacity at constant volume :math:`c_V` in :math:`kJ`, heat capacity at constant pressure :math:`c_P` in :math:`kJ`.

    **xxx_clusters.dat**
        Contains the contributions of each cluster to the quantity denoted by xxx divided by its absolute population (meaning these are cluster specific quantities). Possible quantities are: the partition function and its derivatives, the indistinguishability contribution, Helmholtz free energy, internal energy, Gibbs energy, enthalpy entropy, heat capacity at constant volume and pressure.
        
    **populations.dat**
        Contains populations of each cluster in the order they were specified in the clusterset. Populations are monomer normalized. For example, in a binary system:

.. math::

    N^\prime_\wp = \frac{\left(i_\wp+j_\wp\right)N_\wp}{N_\text{1,tot} + N_\text{2,tot}}.

.. line-block::
        Generally, in a multi-component system:

.. math::

    \qquad N^\prime_\wp = \sum_\mathrm{c} \frac{i_\mathrm{c} \cdot N_\wp}{N_\mathrm{c,tot}}.
       

.. line-block::
    **concentrations.dat**
        Contains concentrations in :math:`mol L^{-1}` of each cluster in the order they were specified in the clusterset. Concentrations are not monomer normalized.

        