Getting Started
###############

What machines will python run on? We have run python various versions of linux and on Mac.
It is compiled using mpicc, with an option to compile with gcc. 

It uses the Gnu Scientific Libraries (gsl)

(Developers should also have cproto in their path in order to create new prototypes, and access to indent to insure that routines are formatted in a standard fashion. They will also want to make sure the py_progs routines are properly installed, as indicated below).

Installation
============

Python and the various routines associated are set up in a self-contained directory structure.
The basic directory structure and the data files that one needs to run Python need to be retrieved and compiled.

If you want to obtain a stable (!) release, go to the `Releases <https://github.com/agnwinds/python/releases/>`_ page.

If you want to download the latest dev version, you can zip up the git repository by clicking on the zip icon to the right of the GitHub page.
Alternatively, clone it directly as

.. code:: bash

    $ git clone https://github.com/agnwinds/python.git

If you anticipate contributing to development we suggest Forking the repository and submitting pull requests with any proposed changes.

Once you have the files, you need to cd to the new directory and set your environment variables

.. code:: bash

    $ export PYTHON = /path/to/python/
    $ cd $PYTHON
    $ ./configure
    $ make install
    $ make clean

One can run a more rigorous clean of GSL with :code:`make distclean`, or remove the compiled GSL libraries altogether with :code:`make rm_lib`.

note that export syntax is for bash- for csh use

.. code:: console

    $ setenv PYTHON /path/to/python/


The atomic data needed to run Python is included in the distribution.  

Running python
==============

To run python you need to add the following to your $PATH variable:

.. code:: bash

    $PYTHON/bin

You can then setup your symbolic links by running

.. code:: bash

    $ Setup_Py_Dir

and run the code by typing, e.g.

.. code:: bash

    $ py root.pf


Auxiliary programs
------------------

There are two programs that are useful for extracting information about models

* windsave2table generates a series of astropy tables that can be used to inspect elements of the various models, including densities of specific ions
* py_wind is a mainly interactive routine that prints similar infomation to the screen.

The two files are run as follows

.. code:: bash

    $ windsave2table root
    $ py_wind root

Brief descriptions of command line options for running these routines can obtained using a -h switch


Directory structure
-------------------

The python directory structure is fairly simple:

source
  Location of source code

bin
  Location of executables

docs 
  Location of documentation, including sphinx docs, doxygen, parameters and documentation for the python programs in py_progs.

data
  Location for all datafiles. Files that are mainly for reference should be gzipped to save space. Such files are not recreated in

bin
  The location of the executables. (It is a good idea to put this directory in your path)

software
  This directory contains libraries which are used in in python that must be recompiled when creating an installation on a new machine, primarily Bill Pence's cfitsio package and the GNU scientific library gsl

py_progs
  python programs for helping analyse the code. We recommend adding this directory to your PATH and PYTHON_PATH environment variables.

examples
  A directory with a few examples of python runs. (Note that the input files will have changed and so one may not be able to run these examples without some changes in the input files.)

Please help by reporting bugs in installation
---------------------------------------------

This can be done by submitting a bug under the `Issues <https://github.com/agnwinds/python/issues/>`_ page
