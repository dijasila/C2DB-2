C2DB version 2
=============

This repository contains shared code to create figures and other
resources that might be useful to share


A C2DB database can be found in
-------------------------------
`c2db-nodata.db`. This database only incluedes key-value-pairs.


Procedure to include figure script
----------------------------------

- Put your code in separate folders with descriptive names.
- Preferably: Make a script ``plot.py`` with a function called
  ``plot`` that can be called to plot your figures.


Various utility functions can be found in
-----------------------------------------

`utils.py`. In order to use these functions you need to use the magic
code listed in the top of that module.
