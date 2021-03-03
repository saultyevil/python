Bound Free
##########


Source
======

Obtained from The `Opacity project <http://cdsweb.u-strasbg.fr/topbase/topbase.html>`_. See also `Cunto et at 1993, A&A, 275, L5 <http://articles.adsabs.harvard.edu/full/1993A%26A...275L...5C>`_


Translation to Python format
============================

ksl - It's not clear that we are now making use of the topbase data in this way but my original attempt to incorporate topbase photoinization data into Python is contained in the directory topbase. Processing of these files was done by py_top_phot. My feeling is that we can replace these remarks with those that are more up to date, once Nick and James discuss this section, and delete any mention of my original attempt on this in the data-gen archive.



Data format
===========
Our original photoionization cross sectiions came from a combination of \cite{verner96}, supplemented by a set of older values from \cite{verner95}

The TopBase photoionization files have the following format::


  PhotTopS  1  1 200    1    13.605698  50
  PhotTop    13.605698 6.304e-18
  PhotTop    16.627193 3.679e-18

whereas the Macro Atoms look like::

  PhotMacS       1       1       1       1       13.598430     100
  PhotMac       13.598430   6.3039999e-18
  PhotMac       13.942675   5.8969998e-18

The meaning of the columns is the same in both cases here. It may be simply a historical accident that we have both formats, or probably we were worried we would need to change. Topbase is generally the source for this information.

The initial line contains (a) the element, (b) the ion, (c) the level, (d) the level up, (e), the energy threshold in eV, and (f)  the number of x-sections to be read in.
The following lines gives the photon energy in eV, and the cross-section in cm$^{2}$.  To be a valid file the photon x-section  must not appear below the energy threshold

"Level up" corresponds to how many levels the electron is moving in the transition: this is simply 1
 
For the simple atom case, the header line can be parsed as follows
 
* z:  the atomic number
* ionization stage, in the astronomical  convention of the lower level ion, the one that gets ionized
* ISLP
* ll: the lower level index in the ion that gets ionized (with that ISLP)
* e: the energy threshold in eV (only photons above this energy ionize)
* npts: the number of points where the cross-section is measured
 
For the simple atom case the combination ISLP and level is unique
 
For the macro-atom case the entries in the header line are
 
* z:  the atomic number
* ul: the upper level index in the ion after ionization (usually 1)
* e: the energy threshold in eV (only photons above this energy ionize)
* npts: the number of points where the cross-ection is measured
 
For the macro atom case, the indices relate to the energy levels, that is these must be connected up properly
 
The TopBase energies are inaccurate and so generally adjustments are made using Chianti or some other source to fix up the energy levels.

Python structure
================
Where the data is stored internally in Python

Comments
========

**Extrapolation to higher energies**

Around about Python 76 we discovered that some topbase cross-sections have maximum energy cutoffs, for no good reason.
This causes unrealistic edges to appear in spectra. To counter this, JM wrote a 
script to extrapolate the cross-section to higher energies. This is done by 
calculating the gradient in log-space at the maximum energy and extrapolating
to 100 keV. A number of cross-sections had unrealistic gradients at the original 
maximum 
energy, and were identified by eye and then forced to have a :math:`\nu^{-3}` shape.
This is the shape of a hydrogenic cross-section and whilst it is not accurate 
for non-hydrogenic ions, it is more realistic (and conservative) than some of 
the unphysically shallow gradients that were being found.
This is also briefly described in section~3.7.2 of Matthews PhD thesis.
The python scripts can be found in progs/extrapolate\_xs/ 
with docstrings describing their use.

