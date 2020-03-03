Importing Wind Models
#####################

In Python it is possible to import a wind model from an arbitrary grid. Currently
supported by Python are 1D spherical grids, or 2D cylindrical or polar grids.

Grid Conventions
================

In order for Python to be able to interpret and read in the provided grid, Python
expects the input files to be formatted in a certain standard. The input format
for the different coordinate systems is mostly standardised, with only small
differences relating to the coordinates for the different systems.

.. admonition :: Cell Quantity Locations

    The grid coordinates and velocities of a cell are defined at the inner vertex
    of cells. But, the electron and radiation temperature, as well as the mass
    density are defined at the center of cells.

Parameters
----------

To use an imported model, when choosing a wind type, you must select `imported`,
i.e.,

`Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported) imported`

The coordinate system must then be specified using the parameter,

`Wind.coord_system(spherical,cylindrical,polar,cyl_var) spherical`

Finally, the path to the model grid to read in must be specified using,

`Wind.model2import 1d_model.import.txt`

Coordinate Systems
------------------

Spherical
^^^^^^^^^

Using a spherical coordinate system, a 1D spherically symmetric model can be
read into Python.

To read in a grid of this type, the following columns are required for each cell:

* i                        :  the element number for each cell
* :math:`r`                :  the radial coordinate in CGS
* :math:`v_{r}`            :  the radial velocity in CGS
* :math:`\rho`             :  the mass density in CGS
* :math:`T_{e}` (optional) :  the electron temperature in Kelvin
* :math:`T_{r}` (optional) :  the radiation temperature in Kelvin

Cylindrical
^^^^^^^^^^^

Using cylindrical coordinates, a 2D model can be read into Python.

.. admonition :: Grid Coordinates

    Note that the grid coordinates and the velocity is specified in Cartesian
    coordinates.

To read in a grid of this type, the following columns are required for each cell:

* i                        :  the i element number (row)
* j                        :  the j element number (column)
* inwind                   :  a flag indicating whether the cell is in the wind or not
* :math:`x`                :  the x coordinate in CGS
* :math:`z`                :  the z coordinate in CGS
* :math:`v_x`              :  the velocity in the x direction in CGS
* :math:`v_y`              :  the velocity in the y direction in CGS
* :math:`v_z`              :  the velocity in the z direction in CGS
* :math:`\rho`             :  the mass density in CGS
* :math:`T_{e}` (optional) :  the electron temperature in Kelvin
* :math:`T_{r}` (optional) :  the radiation temperature in Kelvin


Polar
^^^^^

.. admonition :: Cartesian Velocity

    The velocity in for the polar grid is required to be in Cartesian
    coordinates due to conventions within the Python programming style. As such,
    any polar velocity components must first be projected into their Cartesian
    equivalent.


* i                        :  the i element number (row)
* j                        :  the j element number (column)
* inwind                   :  a flag indicating whether the cell is in the wind or not
* :math:`r`                :  the radial coordinate in CGS
* :math:`\theta`           :  the :math:`\theta` coordinate in CGS
* :math:`v_x`              :  the velocity in the x direction in CGS
* :math:`v_y`              :  the velocity in the y direction in CGS
* :math:`v_z`              :  the velocity in the z direction in CGS
* :math:`\rho`             :  the mass density in CGS
* :math:`T_{e}` (optional) :  the electron temperature in Kelvin
* :math:`T_{r}` (optional) :  the radiation temperature in Kelvin


Setting Values for `inwind`
---------------------------

The `inwind` flag is used to mark if a grid cell is either in the wind or not
in the wind. The following enumerator flags are used,

.. code :: c

    W_IGNORE      = -2   // ignore this grid cell
    W_NOT_INWIND  = -1   // this cell is not in the wind
    W_ALL_INWIND  = 0    // this cell is in the wind

Whilst it is possible to set in `inwind = 1` for a grid cell, that is that the
cell is partially in the wind, Python will instead set these cells with
`inwind = -2` and ignore these grid cells.

Setting Wind Temperatures
-------------------------

Reading in a temperature is optional when importing a model. However, if one
temperature value for a cell is provided, then Python assumes that this is
the electron temperature and the radiation temperature will be initialised as,

.. math ::
    T_{r} = 1.1 T_{e}.

However, if two temperature values are provided for the cells, then the first
temperature will be assumed as being the electron temperature and the second
will be the radiation temperature.

Guard Cells
-----------

In order to calculate velocity gradients at the outer edges of the wind, Python
requires the outer cell boundaries to have "guard cells". This translates into
the last row/column of cells being set as not being in the wind.

For example, if you are importing a spherical wind model with 100 grid cells,
then the final grid cells should not be in the wind, i.e. the radial coordinate
of the first guard cell should be less than the maximum wind radius.

Maximum Wind Radius
-------------------

.. todo :: I'm unclear if this is desired behaviour at the moment

In order for Python to be able to calculate when a Photon has escaped the wind,
a maximum wind radius must be specified as with the default models. This value
should be set whilst taking into account the coordinates of the guard cells.

Tools
=====

Some tools to convert Python `root.wind_save` files into models which can be
imported exist in `$PYTHON/py_progs` are are named,

* import_spherical.py
* import_cyl.py
* import_rtheta.py
