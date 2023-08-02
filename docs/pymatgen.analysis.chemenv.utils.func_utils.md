---
layout: default
title: pymatgen.analysis.chemenv.utils.func_utils.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.func_utils module

This module contains some utility functions and classes that are used in the chemenv package.


### _class_ pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction(function, options_dict=None)
Bases: `object`

Abstract class for all ratio functions.

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {_ )

#### evaluate(value)
Evaluate the ratio function for the given value.


* **Parameters**

    **value** – Value for which ratio function has to be evaluated.



* **Returns**

    Ratio function corresponding to the value.



#### _classmethod_ from_dict(dd)
Construct ratio function from dict.


* **Parameters**

    **dd** – Dict representation of the ratio function



* **Returns**

    Ratio function object.



#### setup_parameters(options_dict)
Set up the parameters for this ratio function.


* **Parameters**

    **options_dict** – Dictionary containing the parameters for the ratio function.



* **Returns**

    None.



### _class_ pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions applied to the continuous symmetry measure (CSM).

Uses “finite” ratio functions.

See the following reference for details:
ChemEnv: a fast and robust coordination environment identification tool,
D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'power2_decreasing_exp': ['max_csm', 'alpha'], 'smootherstep': ['lower_csm', 'upper_csm'], 'smoothstep': ['lower_csm', 'upper_csm']_ )

#### fractions(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



#### mean_estimator(data)
Get the weighted CSM using this CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate the weighted CSM.



* **Returns**

    Weighted CSM from this ratio function.



#### power2_decreasing_exp(vals)
Get the evaluation of the ratio function f(x)=exp(-a\*x)\*(x-1)^2.

The CSM values (i.e. “x”), are scaled to the “max_csm” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### ratios(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



#### smootherstep(vals)
Get the evaluation of the smootherstep ratio function: f(x)=6\*x^5-15\*x^4+10\*x^3.

The CSM values (i.e. “x”), are scaled between the “lower_csm” and “upper_csm” parameters.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### smoothstep(vals)
Get the evaluation of the smoothstep ratio function: f(x)=3\*x^2-2\*x^3.

The CSM values (i.e. “x”), are scaled between the “lower_csm” and “upper_csm” parameters.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



### _class_ pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions applied to the continuous symmetry measure (CSM).

Uses “infinite” ratio functions.

See the following reference for details:
ChemEnv: a fast and robust coordination environment identification tool,
D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'power2_inverse_decreasing': ['max_csm'], 'power2_inverse_power2_decreasing': ['max_csm']_ )

#### fractions(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



#### mean_estimator(data)
Get the weighted CSM using this CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate the weighted CSM.



* **Returns**

    Weighted CSM from this ratio function.



#### power2_inverse_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x.

The CSM values (i.e. “x”), are scaled to the “max_csm” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### power2_inverse_power2_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x^2.

The CSM values (i.e. “x”), are scaled to the “max_csm” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – CSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the CSM values.



#### ratios(data)
Get the fractions from the CSM ratio function applied to the data.


* **Parameters**

    **data** – List of CSM values to estimate fractions.



* **Returns**

    Corresponding fractions for each CSM.



### _class_ pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions applied to differences of
continuous symmetry measures (DeltaCSM).

Uses “finite” ratio functions.

See the following reference for details:
ChemEnv: a fast and robust coordination environment identification tool,
D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'smootherstep': ['delta_csm_min', 'delta_csm_max']_ )

#### smootherstep(vals)
Get the evaluation of the smootherstep ratio function: f(x)=6\*x^5-15\*x^4+10\*x^3.

The DeltaCSM values (i.e. “x”), are scaled between the “delta_csm_min” and “delta_csm_max” parameters.


* **Parameters**

    **vals** – DeltaCSM values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the DeltaCSM values.



### _class_ pymatgen.analysis.chemenv.utils.func_utils.RatioFunction(function, options_dict=None)
Bases: `AbstractRatioFunction`

Concrete implementation of a series of ratio functions.

Constructor for AbstractRatioFunction.


* **Parameters**


    * **function** – Ration function name.


    * **options_dict** – Dictionary containing the parameters for the ratio function.



#### ALLOWED_FUNCTIONS(_: ClassVar[dict[str, list]_ _ = {'inverse_smootherstep': ['lower', 'upper'], 'inverse_smoothstep': ['lower', 'upper'], 'power2_decreasing_exp': ['max', 'alpha'], 'power2_inverse_decreasing': ['max'], 'power2_inverse_power2_decreasing': ['max'], 'smootherstep': ['lower', 'upper'], 'smoothstep': ['lower', 'upper']_ )

#### inverse_smootherstep(vals)
Get the evaluation of the “inverse” smootherstep ratio function: f(x)=1-(6\*x^5-15\*x^4+10\*x^3).

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### inverse_smoothstep(vals)
Get the evaluation of the “inverse” smoothstep ratio function: f(x)=1-(3\*x^2-2\*x^3).

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### power2_decreasing_exp(vals)
Get the evaluation of the ratio function f(x)=exp(-a\*x)\*(x-1)^2.

The values (i.e. “x”), are scaled to the “max” parameter. The “a” constant
correspond to the “alpha” parameter.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### power2_inverse_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x.

The values (i.e. “x”), are scaled to the “max” parameter.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### power2_inverse_power2_decreasing(vals)
Get the evaluation of the ratio function f(x)=(x-1)^2 / x^2.

The values (i.e. “x”), are scaled to the “max” parameter.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### smootherstep(vals)
Get the evaluation of the smootherstep ratio function: f(x)=6\*x^5-15\*x^4+10\*x^3.

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.



#### smoothstep(vals)
Get the evaluation of the smoothstep ratio function: f(x)=3\*x^2-2\*x^3.

The values (i.e. “x”), are scaled between the “lower” and “upper” parameters.


* **Parameters**

    **vals** – Values for which the ratio function has to be evaluated.



* **Returns**

    Result of the ratio function applied to the values.