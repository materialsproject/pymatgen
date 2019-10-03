Change log
==========

v2019.10.3
----------
* Faster get_all_neighbors based on @chc273's improvements. get_all_neighbors
  now returns a Site-like object with nn_distance, image and index attrbutes.
  Much easier to use.
* Bug fix for XCrySDen parser (@stevetorr)
* Added optional mid_struct to direct interpolation (@jmmshn)
