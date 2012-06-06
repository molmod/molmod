Installation instructions
=========================


Preparing your mind
-------------------

MolMod is developed and tested in modern Linux environments. The
installation and usage will therefore be relatively easy on Linux. If you want
to use MolMod on other operating systems such as Windows or OSX, you should
have a minimal computer geek status to get it working. We are always interested
in hearing from your installation adventures.

Preparing your Linux system
---------------------------

Some software packages should be installed before MolMod can be installed or
used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed for MolMod:

* Python 2.5, 2.6 or 2.7 (including the header files): http://www.python.org/doc/
* Numpy 1.0 or later: http://numpy.scipy.org/
* A Fortran and a C compiler supported by the F2PY module in Numpy, e.g.
  gfortran and gcc: http://gcc.gnu.org/
* Git: http://git-scm.com/

Most Linux distributions can install this software with just a single command
on the command line by the administrator. They are listed below for several
popular Linux distributions:

* Ubuntu 10.4::

    sudo apt-get install python python-dev python-numpy gfortran gcc git-core

* Debian 5. You first have to become root because the sudo program is not
  configured by default. ::

    su -
    apt-get install python python-dev python-numpy gfortran gcc git-core
    exit

* Fedora 17::

    sudo yum install python-devel numpy numpy-f2py gcc-gfortran gcc git
    exit

* Suse 11.2::

    sudo zypper install python-devel python-numpy gcc gcc-fortran git

  There seems to be something odd going on with the default Python configuration
  on Suse installations. You have to edit the file
  ``/usr/lib64/python2.4/distutils/distutils.cfg`` or
  ``/usr/lib32/python2.4/distutils/distutils.cfg``, depending on the CPU
  architecture, to comment out the line ``prefix=/usr/local`` with a ``#``
  symbol. Otherwise it is impossible to install Python packages in the home
  directory, as we will do below.

In order to enable the installation and usage of Python packages in the home
directory, as we will do in the next section, one must configure a few
environment variables:

* Bash users: add the following two lines to your ``~/.bashrc`` file::

    export PYTHONPATH=$HOME/lib/python:$HOME/lib64/python:$PYTHONPATH
    export PATH=$HOME/bin:$PATH

* TC Shell users: add the lines to your ``~/.tcshrc`` file::

    setenv PYTHONPATH $HOME/lib/python:$HOME/lib64/python:$PYTHONPATH
    setenv PATH $HOME/bin:$PATH

If you don't know which shell you are using, you are probably using Bash. Note
that some of these lines may already be present. **These settings are only
loaded in new terminal sessions, so close your terminal and open a new one
before proceeding.**


Installing the bleeding edge version of MolMod
----------------------------------------------

The following series of commands will download the latest versions of the
MolMod package, and will then install it into your home directory. Make sure you
execute these commands in some sort of temporary directory. ::

    git clone git://molmod.ugent.be/git/molmod.git
    (cd molmod; ./setup.py install --home=~)

You are now ready to start using MolMod!

A few quick checks
------------------

It may be interesting to double check your installation before proceeding,
unless you `feel lucky`. The MolMod files are installed in the following
directories:

* Scripts: ``~/bin``
* Modules: ``~/lib/python`` or ``~/lib64/python``
* Data: ``~/share``

There should be at least some files present in these directories.

The MolMod modules should be accessible from any Python session. This can be
checked by starting Python interactively and loading the modules manually. There
should be no errors when importing the modules::

    toon@poony ~> python
    Python 2.6.5 (r265:79063, Apr 16 2010, 13:57:41)
    [GCC 4.4.3] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import molmod
    >>> quit()
    toon@poony ~>


Upgrading to the bleeding edge version of MolMod
------------------------------------------------

In case you want to upgrade MolMod to the latests development version after a
previous install, then execute the following commands (in the same directory)::

    (cd molmod; git pull; rm -r ~/lib*/python/molmod*; ./setup.py install --home=~)


Testing your installation
-------------------------

For the development and testing one needs to install one additional package:

* Nosetests >= 0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/

Most Linux distributions can install this software with just a single command on
the command line by the administrator.

* Ubuntu 10.4::

    sudo apt-get install python-nose

* Debian 5 does not have Python 2.6. ::

    su -
    apt-get install python-nose
    exit

* Fedora 17::

    sudo yum install python-nose

* Suse 11.2. ::

    sudo zypper install python-nose

Once these dependecies are installed, go to the temporary directory where the
MolMod source code was downloaded and execute the following commands::

    cd molmod
    ./cleanfiles.sh
    ./setup.py build_ext -i
    nosetests -v

This will run a series of tests to check the validity of the outcomes generated
by MolMod. If some tests fail, post the output of the tests on the `MolMod
mailing list <https://groups.google.com/forum/#!forum/molmodlib>`_. (There are
currently three test that are skipped.)
