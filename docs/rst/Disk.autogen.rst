
====
Disk
====

Disk.rad_type_to_make_wind
==========================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb. Multi-line description, must keep indentation.

models. Multi-line description, must keep indentation.


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_disk.c


Disk.radmax
===========
The outer edge of the disk.  Photons inside this radius are
absorbed or re-radiated.  Photons which are outside this radius
pass through the disk plane.

**Type:** rddoub

**Unit:** cm

**Value:** Greater than 0

**Parent(s):**
  disk.type_: disktype must be 1 or 2, standard or vertically extended disk


**File:** setup_disk.c


Disk.radiation
==============
Multi-line description, must keep indentation.

**Type:** Boolean(yes/no)

**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** setup_disk.c


----------------------------------------

Disk.temperature.profile
------------------------
The choice of disk temperature profile

**Type:** Enum (Int)

**Values:**

standard. A Shakura - Sunyaev  disk, with a hard inner boundar

readin. Read the profile in from a file; the user will be queried for the name of the file

yso. YSO???

analytic. DEPRECATED??? A profile designed for the situation where the disk is being illuminated by star


**Parent(s):**
  Disk.radiation_: This input is requested for all disks that radiate


**File:** setup_disk.c


----------------------------------------

Disk.T_profile_file
^^^^^^^^^^^^^^^^^^^
When the user chooses to read in the temperature profile as a
function of radius, the user is asked the name of the file that
contains the desired profile.

**Type:** String

**Parent(s):**
  Disk.temperature.profile_: Unspecified


**File:** setup_disk.c


Disk.type
=========
Parameter defining whether there is a disk in the system

**Type:** Enum (Int)

**Values:**

none. no.disk

flat. standard.flat.disk

vertically.extended. vertically.extended.disk


**Parent(s):**
  parameter_: This question is always asked


**File:** setup_disk.c


----------------------------------------

Disk.z1
-------
For a vertically extended the disk, the height of the disk is
set to be Disk.rad_mask*(r/Disk.rad_max)**Disk.z1 where Disk.z1
is the power law index

**Type:** rddoub

**Unit:** None

**Value:** Condition e.g. greater than 0 or list e.g. [1, 2, 5]

**Parent(s):**
  Disk.type_: This question is ascked whenever the Disk.type is vertically extended


**File:** setup_disk.c


Disk.z0
=======
fractional.height.at.diskrad.  The physical height at the
outer disk will be this * disk.radmax

**Type:** rddoub

**Unit:** None

**Value:** Greater than 0

**Parent(s):**
  disk_type_: disk_type=vertically extended


**File:** setup_disk.c


Disk.mdot
=========
The mass transfer rate in the disk when considering a standard Shakura-disk.

**Type:** rddoub

**Unit:** msol/yr

**File:** setup_disk.c


Disk.rad_type_in_final_spectrum
===============================
Multi-line description, must keep indentation.

**Type:** Enumerator

**Values:**

bb. Multi-line description, must keep indentation.

models. Multi-line description, must keep indentation.

uniform. Multi-line description, must keep indentation.


**Parent(s):**
  parameter_: Condition e.g. greater than 0 or list e.g. [1, 2, 5]


**File:** python.c


