APOETC
======
This package is used to calculate an exposure time based on a desired signal to noise ratio. See documentation `here`_.

.. _here: https://canasmh.github.io/APOETC


Limitations
------------
This package was intented only for the Astrophysical Research Consortium (ARC) 3.5m telescope at the Apache Point Observatory. As of right now, only the arctic instrument is available when using the arc.Instrument class.

Installation
------------
Use the package manager `pip`_ to install APOETC. First, ``cd`` into the main directory (containing the setup.py file), and then type in the command line

.. _pip: https://pip.pypa.io/en/stable/

``pip install ./``


Importing APOETC
----------------

One can easily import the package (once downloaded) by typing into python

.. code:: python

   import APOETC

While this should work, I recommend importing the modules from APOETC instead importing the actual package. i.e.,

.. code:: python

   from APOETC import *



This way, whenever you call a class or a function from the modules, you only need to write ```module.class``` instead of writing ```APOETC.module.class```.

Using APOETC
------------
As mentioned above, the only instrument that functions in this package is the Arctic imaging instrument. For starters, let's take a look at the arc module.

.. code:: python

   from APOETC import arc
   instrument = arc.Instrument('Arctic')


