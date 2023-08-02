---
layout: default
title: pymatgen.util.num.md
nav_exclude: true
---

# pymatgen.util.num module

This module provides utilities for basic math operations.


### pymatgen.util.num.abs_cap(val, max_abs_val=1)
Returns the value with its absolute value capped at max_abs_val.
Particularly useful in passing values to trigonometric functions where
numerical errors may result in an argument > 1 being passed in.


* **Parameters**


    * **val** (*float*) – Input value.


    * **max_abs_val** (*float*) – The maximum absolute value for val. Defaults to 1.



* **Returns**

    val if abs(val) < 1 else sign of val \* max_abs_val.



### pymatgen.util.num.make_symmetric_matrix_from_upper_tri(val)
Given a symmetric matrix in upper triangular matrix form as flat array indexes as:
[A_xx,A_yy,A_zz,A_xy,A_xz,A_yz]
This will generate the full matrix:
[[A_xx,A_xy,A_xz],[A_xy,A_yy,A_yz],[A_xz,A_yz,A_zz].


### pymatgen.util.num.maxloc(seq)
Return the index of the (first) maximum in seq.

```python
>>> assert maxloc([1,3,2,3]) == 1
```


### pymatgen.util.num.min_max_indexes(seq)
Uses enumerate, max, and min to return the indices of the values
in a list with the maximum and minimum value.


* **Parameters**

    **seq** – A sequence of numbers.



### pymatgen.util.num.minloc(seq)
Return the index of the (first) minimum in seq.

```python
>>> assert minloc(range(3)) == 0
```


### pymatgen.util.num.monotonic(values, mode='<', atol=1e-08)
True if values are monotonically (decreasing|increasing).
mode is “<” for a decreasing sequence, “>” for an increasing sequence.
Two numbers are considered equal if they differ less than atol.

<!-- warning:
Not very efficient for large data sets. -->
```python
>>> values = [1.2, 1.3, 1.4]
>>> monotonic(values, mode="<")
False
>>> monotonic(values, mode=">")
True
```


### pymatgen.util.num.non_decreasing(values)
True if values are not decreasing.


### pymatgen.util.num.non_increasing(values)
True if values are not increasing.


### pymatgen.util.num.round_to_sigfigs(num, sig_figs)
Rounds a number rounded to a specific number of significant
figures instead of to a specific precision.


### pymatgen.util.num.strictly_decreasing(values)
True if values are strictly decreasing.


### pymatgen.util.num.strictly_increasing(values)
True if values are strictly increasing.