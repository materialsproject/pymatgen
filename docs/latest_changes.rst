Change log
==========

v3.2.9
------
* Major PD stability improvements, especially for very high dim hulls with lots
  of entries.
* Improvements to Ewald summation to be close to GULP implementation.
* Deprecate physical constants module in favor of scipy's version.
* Remove many pyhull references to use scipy's ConvexHull implementation.
* Bug fix for sulfide correction.
