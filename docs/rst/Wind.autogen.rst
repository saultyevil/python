
====
Wind
====

Wind.old_windfile
=================
The rootname of a previously saved model in a calculation one wishest to
continue (with the possiblity of making changes to some of the details of
the radiation sources, or to extract spectra from different inclinations)

**Type:** String

**Parent(s):**
  parameter_: Needed whenever on wish to begin a calculation from a previously saved wind file


**File:** python.c


Wind.mdot
=========
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** M☉/year

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** knigge.c


Wind.t.init
===========
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** None

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_domains.c


Wind.coord_system
=================
The coordinate system used for a describing a component of the wind.

**Type:** Enum (Int)

**Values:**

spherical. Spherical

cylindrical. Cylindrical

polar. Spherical polar

cyl_var. Cylindrical varying z


**Parent(s):**
  parameter_: Asked for each component (domain) of the wind


**File:** setup_domains.c


Wind.ionization
===============
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

LTE_te. Multi-line description, must keep indentation.

LTE_tr. Multi-line description, must keep indentation.

ML93. Multi-line description, must keep indentation.

fixed. Multi-line description, must keep indentation.

matrix_bb. Multi-line description, must keep indentation.

matrix_pow. Multi-line description, must keep indentation.

on.the.spot. Multi-line description, must keep indentation.


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup.c


Wind.radiation
==============
Multi-line description, must keep indentation.

**Type:** Boolean (yes/no)

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


Wind.radmax
===========
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_domains.c


Wind.filling_factor
===================
The volume filling factor of the outflow. The implementation
of clumping (microclumping) is described in 
Matthews et al. (2016), 2016MNRAS.458..293M. Asked once per domain. 

**Type:** Enum (Int)

**Values:**

1. smooth

<1. clumped


**Parent(s):**
  parameter_: Asked whenever there is a wind


**File:** setup_domains.c


Wind.dim.in.z_or_theta.direction
================================
Winds are calulated on spherical, cylindrical, or polar grids.
This input variable gives the size of the grid in the z or theta
direction.  Because some grid cells are used as a buffer, the
actual wind cells are contained in a slightly smaller grid than
the number given.

Note that in some situations there may be more than one wind
component, known technically as a domain.  In that case the user
will be queried for this value mulitple times, one for each domain

**Type:** Int

**Value:** Greater than 0

**Parent(s):**
  Number.of.wind.components_: Greater than 0


**File:** setup_domains.c


Wind.number_of_components
=========================
While most simple description of a wind consist of a single region of space, Python can calculate
radiative transfer through more complicated structres, where one region of space is described with one
prescription and another region of space with a second prescription. For example, one might want to place
a disk atomospher between the disk and a wind.  This parameter describes the number of components (aka domains)
of the wind.

**Type:** Int

**Unit:** None

**Value:** greater than 0

**Parent(s):**
  parameter_: This question is always asked.


**File:** python.c


Wind.type
=========
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

SV. Multi-line description, must keep indentation.

corona. Multi-line description, must keep indentation.

homologous. Multi-line description, must keep indentation.

hydro. Multi-line description, must keep indentation.

imported. Multi-line description, must keep indentation.

kwd. Multi-line description, must keep indentation.

shell. Multi-line description, must keep indentation.

star. Multi-line description, must keep indentation.

yso. Multi-line description, must keep indentation.


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_domains.c


Wind.model2import
=================
The name of a file to containing a generic model to read in to python from an ascii file.  (Note
that this is not the same as reading in a model generated by python, but is intended to allow
one to read in a generic model in a variety of formats with only a limited amount of information
required).

**Type:** String

**Parent(s):**
  parameter_: Reqired when a model is to be read in from an external souce


**File:** import.c


Wind.dim.in.x_or_r.direction
============================
Winds are calulated on spherical, cylindrical, or polar grids.
This input variable gives the size of the grid in the x or r
direction.  Because some grid cells are used as a buffer, the
actual wind cells are contained in a slightly smaller grid than
the number given.

Note that in some situations there may be more than one wind
component, known technically as a domain.  In that case the user
will be queried for this value mulitple times, one for each domain

**Type:** Int

**Value:** Greater than 0

**Parent(s):**
  Number.of.wind.components_: Greater than 0


**File:** setup_domains.c


