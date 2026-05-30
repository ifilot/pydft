.. index:: installation

Installation
============

:program:`PyDFT` is available via PyPI. We recommend installing it inside a
virtual environment, which keeps the package and its dependencies separate from
the rest of your system. This prevents accidental upgrades or conflicts that can
break existing setups, and makes it easy to remove the installation later
without touching your system-wide Python.

To create an environment, install :program:`PyDFT`, and activate it, run::

    python3 -m venv venv
    source venv/bin/activate
    pip install --upgrade pip
    pip install pydft

If you also need **isosurface functionality**, install the optional
:program:`pytessel` package::

    pip install pytessel

When you are done working, you can leave the environment with::

    deactivate

If you prefer not to use a virtual environment, you can install for the current
user only with ``pip install --user pydft``, which avoids any system-wide
changes.

.. warning::
    We **do not recommend** using ``sudo pip install ...``. Installing packages
    with administrative privileges can overwrite or conflict with system Python
    components and may break tools your operating system relies on.
