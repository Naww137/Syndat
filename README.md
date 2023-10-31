## [Full Syndat Documentation](https://naww137.github.io/Syndat/build/html/index.html)


If you use this software or the methodology behind it, please cite the following journal article:

Walton, N., Brown, J., Fritsch, W., Brown, D., Nobre, G., Sobes, V., "Methodology for physics-informed generation of synthetic neutron time-of-flight measurement data," Computer Physics Communications., vol. 294, pg. 108927, ISSN 0010-4655, https://doi.org/10.1016/j.cpc.2023.108927.

A pre-print manuscript is also available at:
https://doi.org/10.48550/arXiv.2303.09698.


## Quick and easy install instructions

You can simply install Syndat using pip:
```
pip install git+https://github.com/Naww137/nuc_syndat
```
If a new version of Syndat has been released, you can update it by running:
```
pip install --upgrade git+https://github.com/Naww137/nuc_syndat
```

## If you prefer to install from a cloned directory

Clone this repository in the location of your choice, this will be refered to as the toplevel directory. From the toplevel directory run the following command:

`pip install .`
or
`python -m pip install .`
if your python environment is not active.

This command will install the default branch of nuc_syndat to whatever environment you are working in. You will be required to rebuild/reinstall in order to reflect updates. 

To install using [develoment mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) use the `-e` option with pip install. This will install an editable version of the packages with links to the source code s.t. you do not have to re-build to reflect updates made to the source.
```
pip install -e .
```
or 
`python -m pip install -e .`

## A user example can be found in the examples folder

test
