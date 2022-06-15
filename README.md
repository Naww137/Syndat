# nuc_syndat

Insert desciption of this project.

## Installation & Setup

#### Setup & Use nuc_syndat Package

1. Create a virtual environment
 - It is optional, yet recomended to create a virtual environment as a workspace for this package. If you are unfamiliar with this process or why to do it, see [this article](https://towardsdatascience.com/virtual-environments-for-absolute-beginners-what-is-it-and-how-to-create-one-examples-a48da8982d4b). 
   - To create an alias to easily activate this venv, see this [link](https://wpbeaches.com/make-an-alias-in-bash-or-zsh-shell-in-macos-with-terminal/).
 - The necessary packages and dependencies to use nuc_syndat are listed in the requirements.txt file. Ensure that this file is in the root directory of your virtual environment the run `pip install -r requirements.txt` from that directory.
 - 
 

2. 



to install editable version of the package use `-e` s.t. you do not have to re-build to reflect updates.

#### Documentation with sphinx
This project uses sphinx with autodoc to read documentation embeded with pandas [docstrings](https://pandas.pydata.org/docs/development/contributing_docstring.html#plots-in-examples).

Sphinx can be installed with `pip install -U sphinx`

The next step is to build with `sphinx-quickstart` or `sphinx-build -b html sourcedir builddir`



## Other Resources
The basis for the python package structure used in this project can be found [here](https://packaging.python.org/en/latest/tutorials/packaging-projects/)

The basis for the sphinx documentation used in this project can be found [here](https://betterprogramming.pub/auto-documenting-a-python-project-using-sphinx-8878f9ddc6e9)
  

