
===============
Photon_sampling
===============

Photon_sampling.approach
========================
Choice of whether and how to use stratified sampling in creating photons during the
ionization stage of the calculation.  

**Type:** Enum (Int)

**Value:** 0=T,1=(f1,f2),2=cv,3=yso,4=user_defined,5=cloudy_test,6=wide,7=AGN,8=logarithmic

**Parent(s):**
  parameter_: list e.g. [1, 2, 5]


**File:** bands.c


Photon_sampling.nbands
======================
Python uses stratified samplign to generate photons during the ionization phase.  This
parameter allows the user to define the number of bands for stratified sampling, if s/he
wants to customize the bands used for the generation of photons

**Type:** Int

**Unit:** None

**Value:** greater than 0

**Parent(s):**
  parameter_: This parameter is required whenever the user wants to customize stratified sampling


**File:** bands.c


Photon_sampling.high_energy_limit
=================================
Stratified sampling is used during ionization cycles to generate photons.  This parmeter
specires the high energy limit for the frequenies of photons to be generated..

**Type:** Double

**Unit:** eV

**Value:** greater than 0

**Parent(s):**
  parameter_: This paremeter is requied whenever the user wants to customize energy bands for stratified sampling


**File:** bands.c


Photon_sampling.band_min_frac
=============================
When specifying manually the bands used for generating photons during the ionization pahse, this
parameter specifies the The minimum fraction of photons to be generated in this energy band. 
The number of times this parameter will be reqested depends upon the number of bands.  The summ
of the fractuibs need not sum to 1, in which case the remaining photons will be distributed according
to the luminosity in the energy bands

**Type:** Double

**Unit:** None

**Value:** Greater than 0 and should sum to less than 1 over all bands

**Parent(s):**
  parameter_: This parameter is requested whenever the user manually specifies the bands.


**File:** bands.c


Photon_sampling.band_boundary
=============================
When the user specifies what bands are used for stratfied sampling, this parameter specifies the boundaries
between energy bands in which a minimum fraction of photons will be generated.  The number of times this
parameter is request depends upon the number of energies bands being used.

**Type:** Double

**Unit:** eV

**Value:** Greater than 0

**Parent(s):**
  parameter_: Needed whenever the user chooses to customize the band boundaries


**File:** bands.c


Photon_sampling.low_energy_limit
================================
During the ionization phase, stratified sampling is used to provide good coverage of the full ionizing spectrum. This
parameter sets the lowest envergy (frequency) of for phtoons to be generated whne the user wants to customize the
bands.

**Type:** Double

**Unit:** eV

**Value:** greater than 0

**Parent(s):**
  parameter_: This parameter is required whenever the user wants to customize the bands for stratified smapling in the ionzation phase


**File:** bands.c


