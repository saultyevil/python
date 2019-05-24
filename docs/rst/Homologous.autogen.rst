
==========
Homologous
==========

Homologous.radmin
=================
The starting point of for madel of a homoloous flow, a model in
which the velocity at any radius is proportional to the radius

**Type:** Double

**Unit:** cm

**Value:** greater than 0

**Parent(s):**
  parameter_: A required parameter for defining a homologous flow


**File:** homologous.c


Homologous.radmax
=================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** homologous.c


Homologous.density_exponent
===========================
The power law exponent which defines the decline in density of
a homologous flow as a function of radious.

**Type:** Double

**Value:** greater than 0 for a density that declines with radius

**Parent(s):**
  parameter_: A basic paameter needed to define a homlogous flow


**File:** homologous.c


Homologous.vbase
================
Multi-line description, must keep indentation.

**Type:** Double

**Unit:** cm

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** homologous.c


Homologous.boundary_mdot
========================
The mass loss rate at the base of the wind in a homlogous flow model, a flow
in which the velocity is proporional to the radius.  In general, mdot will
decline with radius, depending on the exponent of the power law that describes
the trend in density.

**Type:** Double

**Unit:** Msol/yr

**Value:** greater than 0

**Parent(s):**
  parameter_: One of the basic parameters for defining a hologous flow in Python


**File:** homologous.c


