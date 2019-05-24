
===
AGN
===

AGN.power_law_cutoff
====================
An advanced option - this is a low frequency cutoff for an 
AGN power law spectrum. It prevents the powerlaw being 
applied to low frequencies and giving an odd SED.

**Type:** Double

**Unit:** Hz

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.power_law_index
===================
The exponent alpha in a power las SED applied to an AGN
central source of the form L_nu=K nu**alpha

**Type:** Double

**Unit:** None

**Value:** Any - but sign is not assumed, so for negative index use a negative value

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.bremsstrahlung_alpha
========================
The frequency exponent alpha in bremstrahlung SED of the form
L_nu=nu**alpha exp(-hnu/kT)

**Type:** Double

**Unit:** None

**Value:** Any - sign is not asssumed so use negative if you want negative

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


AGN.bremsstrahlung_temp
=======================
The temperature T in bremstrahlung SED of the form
L_nu=nu**alpha exp(-hnu/kT)

**Type:** Double

**Unit:** K

**Value:** Greater than 0

**Parent(s):**
  parameter_: None


**File:** setup_star_bh.c


AGN.lamp_post_height
====================
Multi-line description, must keep indentation.

**Type:** rddoub

**Unit:** co.gravitational_radius

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_star_bh.c


