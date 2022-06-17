# nuc_syndat

Insert desciption of this project.

## Installation & Setup

### Create a virtual environment (optional) 
 - If you would like to separate this package from your root python installation, you can create a virtual environment as a workspace. If you are unfamiliar with this process or why to do it, see [this article](https://towardsdatascience.com/virtual-environments-for-absolute-beginners-what-is-it-and-how-to-create-one-examples-a48da8982d4b). 
   - To create an alias to easily activate this venv, see this [link](https://wpbeaches.com/make-an-alias-in-bash-or-zsh-shell-in-macos-with-terminal/).

### Clone and install this repository 

#### Option 1 - 

This option uses pip to read and install the nuc-syndat package and other dependent packages at once. This is a bare-bones set of packages and does not include a text editor or ide such as spyder. This option is recomended for the follow use-cases:
                  - Users only as the install of nuc_syndat is not a development installation
                  - Using this package in it's own venv as it will install the listed version of each dependent package

The necessary packages and dependencies to use nuc_syndat are listed below:
** need to add pandas package to this

```
cycler==0.11.0
fonttools==4.33.3
kiwisolver==1.4.3
matplotlib==3.5.2
nuc-syndat @ git+https://github.com/Naww137/nuc_syndat@188ee1aa959c84746833c0022c68b4f92546b13a
numpy==1.22.4
packaging==21.3
Pillow==9.1.1
pyparsing==3.0.9
python-dateutil==2.8.2
six==1.16.0
```

Copy this text to a file named `requirements.txt`. Ensure that this file is in the root directory of the environment that you want to work in. the run `pip install -r requirements.txt` from that directory. 

#### Option 2 -

The following command will clone and install the nuc_syndat package to whatever environment you are working in. You will be required to rebuild/reinstall in order to reflect updates.
```
pip install git+https://github.com/Naww137/nuc_syndat
```

To install using [develoment mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) use the `-e` option with pip install. This will install an editable version of the packages with links to the source code s.t. you do not have to re-build to reflect updates made to the source.
```
pip install -e git+https://github.com/Naww137/nuc_syndat#egg=nuc_syndat
```

### Documentation with sphinx
This project uses sphinx with autodoc to read documentation embeded with pandas [docstrings](https://pandas.pydata.org/docs/development/contributing_docstring.html#plots-in-examples).

Sphinx can be installed with `pip install -U sphinx`

The next step is to build with `sphinx-quickstart` or `sphinx-build -b html sourcedir builddir`
HTML files, with the style and other settings in conf.py, can be generated with `make html`.


 
## Other Resources
The basis for the python package structure used in this project can be found [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

The basis for the sphinx documentation used in this project can be found [here](https://betterprogramming.pub/auto-documenting-a-python-project-using-sphinx-8878f9ddc6e9)
  

