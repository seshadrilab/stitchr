Getting started
===============

Installation
------------

``stitchr`` runs on Python 3.9+ and installation requires pip. Activate a virtual environment with these available if you have one, then clone the repository into wherever you keep your GitHub repos and install:

.. code:: bash

   $ git clone https://github.com/seshadrilab/stitchr.git
   $ cd stitchr
   $ pip install .

This will install all required dependencies, including ``IMGTgeneDL``, ``biopython``, and ``PySide6`` (for the GUI).

``IMGTgeneDL`` can be used via the ``stitchrdl`` command to download suitably formatted data sets to the required directory like so:

``stitchrdl -s human``

See the :ref:`species-covered-label` section for details on the species for which data can be downloaded in this manner.

Running GUI-stitchr
^^^^^^^^^^^^^^^^^^^

After installation, launch the graphical interface with:

.. code:: bash

   $ gui_stitchr

Quick start example
-------------------

The only required fields are the minimal components describing a single rearranged TCR chain: V gene name, J gene name, and CDR3 sequence (either DNA or amino acids). Constant regions must also be specified for all non-human/non-mouse species.

.. code:: bash

   stitchr -v [IMGT V gene] -j [IMGT J gene] -cdr3 [CDR3aa]

   stitchr -v TRBV7-3*01 -j TRBJ1-1*01 -cdr3 CASSYLQAQYTEAFF

   stitchr -v TRAV1-2 -j TRAJ33 -cdr3 TGTGCTGTGCTGGATAGCAACTATCAGTTAATCTGG

See the :ref:`usage-label` section for more detailed usage instructions. ``stitchr`` can also be run in a high-throughput manner (see :ref:`thimble-label`), or via a simple graphical user interface (see :ref:`gui-label`).