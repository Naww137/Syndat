# nuc_syndat

Insert desciption of this project.


After following the instructions to clone and install this project, some of the functionality of nuc-syndat is highlighted in the examples.py script. This script shows the intended workflow and how to use the modules in syndat.


## Setup your environment

### Create a virtual environment (optional) 
If you would like to separate this package from your root python installation, you can create a virtual environment as a workspace. If you are unfamiliar with this process or why to do it, see [this article](https://towardsdatascience.com/virtual-environments-for-absolute-beginners-what-is-it-and-how-to-create-one-examples-a48da8982d4b). 
To create an alias to easily activate this venv, see these links for [bash](https://wpbeaches.com/make-an-alias-in-bash-or-zsh-shell-in-macos-with-terminal/) or [zsh](https://linuxhint.com/configure-use-aliases-zsh/) shells.

### Install dependencies

If you have created a clean, virtual environment, you may use requirements to install all dependencies necessary to run nuc_syndat. If you have not created a virtual environment, this is not recommended as it may change the version of some python modules in your root environment while nuc-syndat will likely run on the most updated versions of these dependent packages. If you want to run nuc_syndat in your root environment, skip this step and go to step 2 but be aware that you may need to install some dependent packages (i.e. scipy, pandas).

This option uses pip to read and install the nuc-syndat package and other dependent packages at once. This is a bare-bones set of packages and does not include a text editor or ide such as spyder. This option is recomended for the follow use-cases:
                  - Users only as the install of nuc_syndat is not a development installation
                  - Using this package in it's own venv as it will install the listed version of each dependent package

The necessary dependencies to use nuc_syndat are listed below:

```
cycler==0.11.0
fonttools==4.33.3
kiwisolver==1.4.3
matplotlib==3.5.2
numpy==1.22.4
packaging==21.3
pandas==1.4.2
Pillow==9.1.1
pyparsing==3.0.9
python-dateutil==2.8.2
pytz==2022.1
scipy==1.8.1
six==1.16.0
```

Copy this text to a file named `requirements.txt`. Ensure that this file is in the root directory of the environment that you want to work in. the run `pip install -r requirements.txt` from that directory. 

## Installing the nuc-syndat package

#### Option - 1
Clone this repository in the location of your choice, this will be refered to as the toplevel directory. From the toplevel directory run the following command:

`pip install .`
or
`python -m pip install .`
if your python environment is not active.

This command will install the default branch of nuc_syndat to whatever environment you are working in. You will be required to rebuild/reinstall in order to reflect updates.

#### Option - 2
To install using [develoment mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) use the `-e` option with pip install. This will install an editable version of the packages with links to the source code s.t. you do not have to re-build to reflect updates made to the source.
```
pip install -e .
```
or 
`python -m pip install -e .`


## Documentation with sphinx
This project uses sphinx with autodoc to read documentation embeded with pandas [docstrings](https://pandas.pydata.org/docs/development/contributing_docstring.html#plots-in-examples).

Sphinx can be installed with `pip install -U sphinx`

The next step is to build with `sphinx-quickstart` or `sphinx-build -b html sourcedir builddir`
HTML files, with the style and other settings in conf.py, can be generated with `make html`.


 
## Other Resources
The basis for the python package structure used in this project can be found [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

The basis for the sphinx documentation used in this project can be found [here](https://betterprogramming.pub/auto-documenting-a-python-project-using-sphinx-8878f9ddc6e9)


  
## Examples

For examples on how to use this package, see the examples folder from the main directory. Within that folder, the "exec" ipython notebook steps through the generation and visualization of a sythetic experimental dataset. 

