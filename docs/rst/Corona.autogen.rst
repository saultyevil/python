
======
Corona
======

Corona.radmax
=============
The corona is a box-shaped region which sits immediately
above the disk.  radmax defines the outer edge of the box.

**Type:** Double

**Unit:** cm

**Value:** Greater than the radius of the central object

**Parent(s):**
  Wind.type_: Selected as one of the choces for this variable


**File:** corona.c


Corona.base_den
===============
The coronal model is defined in terms of a base density
and a scale height

**Type:** Double

**Unit:** number/cm**3

**Value:** Greater than 0

**Parent(s):**
  Wind.type_: Selected as one of the choces for this variable


**File:** corona.c


Corona.radmin
=============
The corona is a box-shaped region which sits immediately
above the disk.  radmin defines the inner edge of the box.

**Type:** Double

**Unit:** cm

**Value:** Greater than that the radius of the central object

**Parent(s):**
  Wind.type_: Selected as one of the choces for this variable


**File:** corona.c


Corona.zmax
===========
The corona is a box-shaped region which sits immediately
above the disk.  zmax defines the height of the box.

**Type:** Double

**Unit:** cm

**Value:** Greater than that the radius of the central object

**Parent(s):**
  Wind.type_: Selected as one of the choces for this variable


**File:** corona.c


Corona.scale_height
===================
The coronal model is defined in terms of a base density
and a scale height

**Type:** Double

**Unit:** cm

**Value:** Greater than 0

**Parent(s):**
  Wind.type_: Selected as one of the choces for this variable


**File:** corona.c


Corona.vel_frac
===============
For the coronal model, the azimuthal velocity is
given by the velocity of the underlying disk.  One
can also give the corona a radial velocity, which is
a fraction of the disk velocity.  (As coded, if this
number is positive, the velicty is the r direction is
toward the central object).

**Type:** Double

**Unit:** None

**Value:** Any, 0 implies no radial velocity.

**Parent(s):**
  Wind.type_: Selected as one of the choces for this variable


**File:** corona.c


