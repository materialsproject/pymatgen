# Introduction

This docmentation provides a guide for pymatgen administrators. The following assumes you are using miniconda or Anaconda.

# Releases

The general procedure to releasing pymatgen comprises three steps.

1. Wait for all unittests to pass on CircleCI.
2. Release PyPI versions + doc.
3. Release conda versions.
4. Release Dash documentation.

## Initial setup

Pymatgen uses [invoke](http://www.pyinvoke.org/) to automate releases. You will also need sphinx and dash2doc. Install these using:

```
pip install --upgrade invoke sphinx doc2dash
```

For 2018, we will release both py27 and py37 versions of pymatgen. Create environments for py27 and py37 using conda.

```
conda create --yes -n py37 python=3.7
conda create --yes -n py27 python=2.7
```

For each env, install some packages using conda followed by dev install for pymatgen.

```
source activate py37
conda install --yes numpy scipy sympy matplotlib
pip install invoke sphinx doc2dash
python setup.py develop
source activate py27
conda install --yes numpy scipy sympy matplotlib
pip install invoke sphinx doc2dash
python setup.py develop
```

## Doing the release

```
source activate py37
invoke release --notest --nodoc
source activate py27
python setup.py bdist_wheel
twine upload dist/*p27*
invoke update-doc
source deactivate
python setup.py develop
```