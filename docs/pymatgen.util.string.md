---
layout: default
title: pymatgen.util.string.md
nav_exclude: true
---

# pymatgen.util.string module

This module provides utility classes for string operations.


### _class_ pymatgen.util.string.Stringify()
Bases: `object`

Mix-in class for string formatting, e.g. superscripting numbers and symbols or superscripting.


#### STRING_MODE(_ = 'SUBSCRIPT_ )

#### to_html_string()
Generates a HTML formatted string. This uses the output from to_latex_string to generate a HTML output.
:return: HTML formatted string.


#### to_latex_string()
Generates a LaTeX formatted string. The mode is set by the class variable STRING_MODE, which defaults to
“SUBSCRIPT”. E.g., Fe2O3 is transformed to Fe$_{2}$O$_{3}$. Setting STRING_MODE to “SUPERSCRIPT” creates
superscript, e.g., Fe2+ becomes Fe^{2+}. The initial string is obtained from the class’s __str__ method.


* **Returns**

    String for display as in LaTeX with proper superscripts and subscripts.



#### to_pretty_string()

* **Returns**

    A pretty string representation. By default, the __str__ output is used, but this method can be
    overridden if a different representation from default is desired.



#### to_unicode_string()

* **Returns**

    Unicode string with proper sub and superscripts. Note that this works only with systems where the sub
    and superscripts are pure integers.



### pymatgen.util.string.charge_string(charge, brackets=True, explicit_one=True)
Returns a string representing the charge of an Ion. By default, the
charge is placed in brackets with the sign preceding the magnitude, e.g.,
‘[+2]’. For uncharged species, the string returned is ‘(aq)’.


* **Parameters**


    * **charge** – the charge of the Ion


    * **brackets** – whether to enclose the charge in brackets, e.g. [+2]. Default: True


    * **explicit_one** – whether to include the number one for monovalent ions, e.g.
    +1 rather than +. Default: True



### pymatgen.util.string.disordered_formula(disordered_struct, symbols=('x', 'y', 'z'), fmt='plain')
Returns a formula of a form like AxB1-x (x=0.5)
for disordered structures. Will only return a
formula for disordered structures with one
kind of disordered site at present.


* **Parameters**


    * **disordered_struct** – a disordered structure


    * **symbols** – a tuple of characters to use for


    * **subscripts** (*'x'**, **'y'**, **'z'*) –


    * **is** (*by default this*) –


    * **disordered** (*but if you have more than three*) –


    * **added** (*species more symbols will need to be*) –


    * **fmt** (*str*) – ‘plain’, ‘HTML’ or ‘LaTeX’


Returns (str): a disordered formula string


### pymatgen.util.string.formula_double_format(afloat, ignore_ones=True, tol: float = 1e-08)
This function is used to make pretty formulas by formatting the amounts.
Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.


* **Parameters**


    * **afloat** (*float*) – a float


    * **ignore_ones** (*bool*) – if true, floats of 1 are ignored.


    * **tol** (*float*) – Tolerance to round to nearest int. i.e. 2.0000000001 -> 2



* **Returns**

    A string representation of the float for formulas.



### pymatgen.util.string.htmlify(formula)
Generates a HTML formatted formula, e.g. Fe2O3 is transformed to
Fe<sub>2</sub>O</sub>3</sub>.

Note that Composition now has a to_html_string() method that may
be used instead.


* **Parameters**

    **formula** –



* **Returns**




### pymatgen.util.string.latexify(formula)
Generates a LaTeX formatted formula. E.g., Fe2O3 is transformed to
Fe$_{2}$O$_{3}$.

Note that Composition now has a to_latex_string() method that may
be used instead.


* **Parameters**

    **formula** (*str*) – Input formula.



* **Returns**

    Formula suitable for display as in LaTeX with proper subscripts.



### pymatgen.util.string.latexify_spacegroup(spacegroup_symbol)
Generates a latex formatted spacegroup. E.g., P2_1/c is converted to
P2$_{1}$/c and P-1 is converted to P$\\overline{1}$.

Note that SymmetryGroup now has a to_latex_string() method that may
be called instead.


* **Parameters**

    **spacegroup_symbol** (*str*) – A spacegroup symbol



* **Returns**

    A latex formatted spacegroup with proper subscripts and overlines.



### pymatgen.util.string.str_delimited(results, header=None, delimiter='\\t')
Given a tuple of tuples, generate a delimited string form.
>>> results = [[“a”,”b”,”c”],[“d”,”e”,”f”],[1,2,3]]
>>> print(str_delimited(results,delimiter=”,”))
a,b,c
d,e,f
1,2,3.


* **Parameters**


    * **results** – 2d sequence of arbitrary types.


    * **header** – optional header


    * **delimiter** – Defaults to “t” for tab-delimited output.



* **Returns**

    Aligned string output in a table-like format.



### pymatgen.util.string.stream_has_colours(stream)
True if stream supports colours. Python cookbook, #475186.


### pymatgen.util.string.transformation_to_string(matrix, translation_vec=(0, 0, 0), components=('x', 'y', 'z'), c='', delim=',')
Convenience method. Given matrix returns string, e.g. x+2y+1/4
:param matrix
:param translation_vec
:param components: either (‘x’, ‘y’, ‘z’) or (‘a’, ‘b’, ‘c’)
:param c: optional additional character to print (used for magmoms)
:param delim: delimiter
:return: xyz string.


### pymatgen.util.string.unicodeify(formula)
Generates a formula with unicode subscripts, e.g. Fe2O3 is transformed
to Fe₂O₃. Does not support formulae with decimal points.

Note that Composition now has a to_unicode_string() method that may
be used instead.


* **Parameters**

    **formula** –



* **Returns**




### pymatgen.util.string.unicodeify_spacegroup(spacegroup_symbol)
Generates a unicode formatted spacegroup. E.g., P2$_{1}$/c is converted to
P2₁/c and P$\\overline{1}$ is converted to P̅1.

Note that SymmetryGroup now has a to_unicode_string() method that
may be called instead.


* **Parameters**

    **spacegroup_symbol** (*str*) – A spacegroup symbol as LaTeX



* **Returns**

    A unicode spacegroup with proper subscripts and overlines.



### pymatgen.util.string.unicodeify_species(specie_string)
Generates a unicode formatted species string, with appropriate
superscripts for oxidation states.

Note that Species now has a to_unicode_string() method that
may be used instead.


* **Parameters**

    **specie_string** (*str*) – Species string, e.g. O2-



* **Returns**

    Species string, e.g. O²⁻