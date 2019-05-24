
======
Reverb
======

Reverb.matom_lines
==================
Number of macro-atom lines to track paths for individually. This many
reverb.matom_line entries are required, and the line associated with each has
the path of photons deexciting into it recorded in its own array. Note: This
doesn't give rise to any noticable differences to the pure wind mode in most
simulations.

**Type:** Int

**Value:** 0 or N

**Parent(s):**
  reverb.type_: 3


**File:** setup_reverb.c


Reverb.filter_lines
===================
Whether or not to filter any lines out of the output file. This is used to keep output
file sizes down, and avoid them overwhelming the user.

**Type:** Int

**Values:**

0. **No filtering**
   
   Include *all* photons that contribute to the spectra in the output
   file. Not recommended as it leads to gargantuan file sizes.

-1. **Filter continuum**
   
   Include all photons whose last interaction was scatter
   or emission in a line. Recommended setting for exploratory runs where you'd
   like to identify which lines are the easiest to process.

N. **Filter lines**
   
   Include N reverb.filter_line entries, each specifying one
   line to keep in the output file. If reverb.matom_lines is >0, all macro-atom
   lines of interest are automatically included in the filter list.


**Parent(s):**
  reverb.type_: Greater than 0


**File:** setup_reverb.c


Reverb.path_bins
================
Number of bins for photon paths. Reverb modes that record the distribution of
path lengths in every wind cell bin them in this number of bins. Bins are
logarithmically spaced between the minimum scale in the system (the smallest
'minimum radius' in any domain) and the 10 * the maximum scale in the system
(10 * the 'maximum radius' in any domain). Default value is 1000, going much
higher does not lead to qualitative differences in TF, going lower makes the
bin boundaries show up in the TF.

**Type:** Int

**Value:** Greater than 0

**Parent(s):**
  reverb.type_: 2, 3


**File:** setup_reverb.c


Reverb.dump_cell
================
Position for a cell, listed as a pair of R:Z coordinates. Will accept any
position that falls within a grid, will error out on ones that don't. This can
be slightly awkward and you may want to run a quick test then use py_wind to
idenfity where wind locations are.

**Type:** Float:Float


**Unit:** cm:cm


**Value:** >0:>0


**Parent(s):**
  reverb.dump_cells_: Greater than 0


**File:** setup_reverb.c


Reverb.type
===========
Whether to perform reverberation mapping. Reverberation mapping tracks the
path of photons emitted in the simulation as they travel through the geometry,
assuming that any delays from recombination etc. are negligible and all delays
are due to light travel time. For each final spectrum, all contributing
photons are output to a '.delay_dump' file that can then be processed using
our 'tfpy' Python (no relation) library.

**Type:** Enum (Int)

**Values:**

none. **Off**

photon. **Simple 'photon' mode**
   
   Each photon is assigned an initial path based on its distance from the
   central source (assuming emission in the disk and wind is correlated with
   emission from the CO).

wind. **Wind mode**
   
   CO photons are assigned paths as in Photon mode, disk photons are assigned
   paths as set by the reverb.disk_type parameter. Photons generated in the
   wind are assigned a path based on the *distribution* of paths of photons
   that have contributed to continuum absorption in that cell.

matom. **Macro-atom mode**
   
   This works as wind mode, but for a number of specified macro-atom lines
   paths are tracked for those photons who cause a deexcitation into a given
   line. When a photon is emitted in one of those lines, the path is drawn from
   that specific distribution. This distribution is build up not just from the
   last cycle of the simulation, but from all cycles after the wind achieves
   >90% convergence. This is necessary as some lines are poorly-sampled.
   
   This mode gives pretty much identical results to wind mode, but at least we
   made it to check rather than just assuming it would be fine.


**File:** setup_reverb.c


Reverb.matom_line
=================
Specifies a line associated with a given macro-atom transition. The species
and transition involved are specified. The internal line associated with this
transition will be printed to standard-out for use when processing outputs. A
line is specified as Element:Ion:Upper level:Lower level.

**Type:** Int:Int:Int:Int


**Value:** >0:>0:>1:>0


**Parent(s):**
  reverb.matom_lines_: Greater than 0


**File:** setup_reverb.c


Reverb.dump_cells
=================
Number of cells to dump. When dumping the path distribution info for a range
of cells, this specifies the number of lines of reverb.dump_cell that will be
provided.

**Type:** Int

**Value:** 0 or N

**Parent(s):**
  reverb.visualisation_: 2, 3


**File:** setup_reverb.c


Reverb.filter_line
==================
Line number of one line to include in the output .delay_dump file. This is
the python internal line number. It can be found using either the macro-atom
mode (which prints out the line number once it's found one) or by doing an
exploratory run with reverb.filter_lines = -1, then looking through the delay
dump file for photons of the right wavelength to see what their line is. This
should almost certainly be changed to be specified using a species and
wavelength!

**Type:** Int

**Value:** Any valid line index

**Parent(s):**
  reverb.filter_lines_: Greater than 0


**File:** setup_reverb.c


Reverb.visualisation
====================
Which type of visualisation to output, if any. Reverb modes that keep arrays
of photon paths per cell can output them either as averages in a 3d model, or
as a selection of flat text files with full bin-by-bin breakdowns. Useful for
diagnostics.

**Type:** Enum (Int)

**Values:**

none. None

vtk. **Mesh visualisation**
   
   Outputs mean incident path per cell, photon count per cell, and mean
   observed delay to '.vtk' format, readable using a range of programs including
   (my preferred option) VisIt, available at https://visit.llnl.gov/.

dump. **Dump cells**
   
   Outputs distributions of paths for continuum heating and each line to a range of 'dump cells'
   specified by X & Z position using the reverb.dump_cells/reverb.dump_cell options.

both. **Both**


**Parent(s):**
  reverb.type_: 2, 3


**File:** setup_reverb.c


Reverb.disk_type
================
Setting for how photons generated in the disk are treated when generating path
distributions for wind cells.

**Type:** Enum (Int)

**Values:**

correlated. **Correlated**
   
   This mode assumes that disk emission is correlated with the
   central source. Photons generated in the disk start with a delay equal to
   the direct distance to the central source. This assumes that the ionisation
   state and luminosity of the disk surface layer is mostly determined by
   unscattered photons from the central source.

uncorrelated. **Uncorrelated**
   
   This mode generates photons with a delay of 0 wherever in the
   disk they come from. This mode is of slightly questionable use and should be
   ignored in preference to 0 or 2. It will, in practise, generally work out
   similar to type 0 as most of the UV photons are generated close-in to the CO.

ignore. **Ignored**
   
   This mode assumes that disk photons do *not* correlate
   with the central source (i.e. disk surface  ionisation state and emissivity is
   driven not by irradiation from the CO but by the mass inflow). This means that
   whilst they contribute to heating the wind, they do not strongly contribute to
   the lags for a given line. Photons generated by the disk do not contribute to
   the path distributions in the wind in this mode.
   
   By removing the (generally) short-delay disk photons from the wind path
   distributions, this will slightly bias them towards the longer delays
   associated with wind self-heating/excitation.


**Parent(s):**
  reverb.type_: 2, 3


**File:** setup_reverb.c


Reverb.angle_bins
=================
Used when generating 3d .vtk output files for visualisation. Sets the number
of angle bins used in the output. Aesthetic only; bigger makes prettier meshes
with larger filesizes.

**Type:** Int

**Value:** Greater than 0

**Parent(s):**
  reverb.visualisation_: 1, 3


**File:** setup_reverb.c


