Macro Atoms
-------------------------

.. todo:: This page shoud contain a description of how macro atoms work. The below is copied from JM's thesis.

The macro-atom scheme was created by Leon Lucy and is outlined in 
his 2002/03 papers. It was implemented in Python by Stuart Sim, initially
for the study of recombination lines in YSOs (Sim et al. 2005).

Lucy (2002,2004) hereafter L02, L03) has shown that it is possible to calculate the emissivity of a gas in statistical equilibrium without approximation for problems with large departures from LTE. His `macro-atom' scheme allows for all possible transition paths from a given level, dispensing with the two-level approximation, and provides a full non-LTE solution for the level populations based on Monte Carlo estimators. The macro-atom technique has already been used to model Wolf-Rayet star winds (Sim 2004), AGN disc winds (Sim et al. 2008), supernovae (Kromer and Sim 2009, Kerzendorf and Sim 2014) and YSOs (Sim et al. 2005). A full description of the approach can be found in L02 and L03. 

Following L02, let us consider an atomic species interacting with a radiation field.
If the quantity :math:`\epsilon_j` represents the ionization plus excitation energy of 
a level :math:`j` then the rates at which the level absorbs and emits radiant energy 
are given by

.. math::

    \dot{A}_{j}^{R} = R_{l j} \epsilon_{j l} 

and

.. math::

    \dot{E}_{j}^{R} = R_{j l} \epsilon_{j l} 


where :math:`\epsilon_{j l} = \epsilon_j - \epsilon_l`.
Here, we have adopted Lucy's convention, in which the subscript 
:math:`l` denotes a summation over all lower states (:math:`<j`), and
(:math:`u`) a summation over all upper states (:math:`>j`).
Similarly, the rates corresponding to _kinetic_ (collisional)
energy transport can then be written as

.. math::

 \dot{A}_{j}^{C} = C_{l j} \epsilon_{j l}

and


.. math::

  \dot{E}_{j}^{C} = C_{j l} \epsilon_{j},

Let us define :math:`{\cal R}` as a total rate, such that
:math:`{\cal R}_{l j}  = R_{l j} + C_{l j}`.
If we now impose statistical equilibrium

.. math::

    ({\cal R}_{l j}-{\cal R}_{j l})+({\cal R}_{u j}-{\cal R}_{ju})=0 \;\;\;,
\end{equation}
we obtain 

.. math::

 \dot{E}_{j}^{R}+\dot{E}_{j}^{C}+{\cal R}_{ju}\epsilon_{j}+
 {\cal R}_{j \ell}\epsilon_{l}  \nonumber \\  
 = \dot{A}_{j}^{R}+\dot{A}_{j}^{C}+{\cal R}_{u j} \epsilon_{j}
 +{\cal R}_{l j} \epsilon_{l}.  

This equation is the starting point for the macro-atom scheme. It shows that, when assuming radiative equilibrium, the energy flows through a system depend only on the transition probabilities and atomic physics associated with the levels the energy flow interacts with. By quantising this energy flow into radiant (:math:`r-`) and kinetic (:math:`k-`) packets, we can simulate the energy transport through a plasma discretised into volume elements (macro-atoms),whose associated transition probabilities govern the interaction of radiant and kinetic energy with the ionization and excitation energy associated with the ions of the plasma.

Although the equation above assumes strict radiative equilbrium,it is trivial to adjust it to include non-radiative source and sink terms. For example, in an expanding parcel of plasma, adiabatic cooling may be included with a simple modification to the RHS.

A Hybrid Scheme
=============================

.. todo:: Describe.
