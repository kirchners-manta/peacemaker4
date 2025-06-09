Peacemaker performs :math:`a_{mf}`, :math:`b_{xv}` parameter sampling on a grid, which can be specified in the ``[qce]`` section of the input file.  
For each pair, the quality of the resulting isobar is compared to certain experimental quantities.  
The following options are available for this purpose: single density, isobar, temperature of phase transition.  
The isobar quality is computed according to the following equation:

.. math::

   \mathrm{error} =
     w_\mathrm{density}
     \left(\frac{\rho - \rho^\mathrm{exp}}{\rho^\mathrm{exp}}\right)^2
     + w_\mathrm{isobar}
     \frac{1}{N} \sum_{i=1}^{N}
     \left(\frac{V_i - V_i^\mathrm{exp}}{V_i^\mathrm{exp}}\right)^2
     + w_\mathrm{phase\ transition}
     \left(\frac{T_\mathrm{pt} - T_\mathrm{pt}^\mathrm{exp}}{T_\mathrm{pt}^\mathrm{exp}}\right)^2

Any combination of the experimental data above can be chosen.  
The relative importance of each quantity can be specified by the weight :math:`w`.

Isobars are specified by an *isobar file*.  
This file shall contain two columns of numbers: temperatures in K in column one and volumes in :math:`\mathrm{dm^3}` in column two.  
All temperatures must be within the temperature range specified in the ``[qce]`` section.  
There are no requirements on the order of the temperatures.  
Temperatures may be included multiple times to put special weight on them.  
If a reference temperature is not equal to the temperature specified by the temperature range in the ``[qce]`` section,  
linear interpolation between the two closest temperatures is performed.
