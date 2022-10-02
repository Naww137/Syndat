How to Use Syndat
=================


Installation
------------
To use **Syndat**, you can install install it using pip:

:code:`pip install git+https://github.com/Naww137/nuc_syndat`

If a new version of Syndat has been released, you can update Syndat by running:

:code:`pip install --upgrade git+https://github.com/Naww137/nuc_syndat`

If you would like to run in `development mode <https://setuptools.pypa.io/en/latest/userguide/development_mode.html>`_, install using:

:code:`pip install -e`

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


Methodology
-----------
Write up some of the methodology here.



.. _usersguide:

Users Guide
-----------

