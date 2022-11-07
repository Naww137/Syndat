Installation
============



Setup Your Environment
----------------------

If you would like to separate this package from your root python installation, you can create a virtual environment as a workspace. 
If you are unfamiliar with this process or why to do it, see [this article](https://towardsdatascience.com/virtual-environments-for-absolute-beginners-what-is-it-and-how-to-create-one-examples-a48da8982d4b). 
To create an alias to easily activate this venv, see these links for [bash](https://wpbeaches.com/make-an-alias-in-bash-or-zsh-shell-in-macos-with-terminal/) or [zsh](https://linuxhint.com/configure-use-aliases-zsh/) shells.



Quick & Easy Install
--------------------

To use **Syndat**, you can simply install it using python's command line utility pip:

:code:`pip install git+https://github.com/Naww137/nuc_syndat`

This will install the **Syndat** package in the site-packages directory for the currently active python environment you are in.

If a new version of Syndat has been released, you can update Syndat by running:

:code:`pip install --upgrade git+https://github.com/Naww137/nuc_syndat`

If you would like to run in `development mode <https://setuptools.pypa.io/en/latest/userguide/development_mode.html>`_, install using:

:code:`pip install -e`


Install From a Cloned directory
-------------------------------

Clone this repository in the location of your choice, this will be refered to as the toplevel directory. From the toplevel directory run the following command:

:code:`pip install .`

or

:code:`python -m pip install .`

if your python environment is not active.

This command will install the default branch of nuc_syndat to whatever environment you are working in. You will be required to rebuild/reinstall in order to reflect updates.

To install using [develoment mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) use the `-e` option with pip install. This will install an editable version of the packages with links to the source code s.t. you do not have to re-build to reflect updates made to the source.

:code:```pip install -e .```

or 

:code:`python -m pip install -e .`




Documentation
-------------

This project comes with significant documentation. In order to devlop and compile this documentation with `sphinx <https://www.sphinx-doc.org/en/master/>`_
follow the next few steps. If you are only a user, all of the Syndat documentation is hosted on this site.

First make sure you have sphinx installed as well as the cloud theme

:code:`pip install sphinx`

:code:`pip install cloud-sptheme`

From the directory where you cloned ADAM run:

:code:`sphinx-build -b html docs/source/ docs/build/html`

To update documetation from the restructured text files without having to rebuild, navigate to the docs/ directory and run:

:code:`make html`

Any changes to the conf.py file will not take effect unless you rebuild.
The HTML files will be created in docs/build/html. The main page can be reached by openning the index.html file. 
This will launch an HTML page and any of the subsequent pages can be reached from here.

