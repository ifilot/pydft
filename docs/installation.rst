.. _installation:
.. index:: Installation

Installation
============

:program:`PyDFT` is available via both Anaconda as well as PyPi. If you are using
Windows and have little experience with Python, we recommend to use Anaconda.
Anaconda is a complete Python package containing many modules and useful programs.
Anaconda can be obtained `via this link <https://www.anaconda.com/download>`_.

.. note::

	If you run into trouble installing :program:`PyDFT` in Anaconda, an easy
	solution can be to create a separate environment for Anaconda. Please consult
	the "Troubleshooting" section as seen below.

Anaconda
--------

Open a Anaconda command prompt and run the following command:

.. code:: bash

	conda install -c ifilot pydft pyqint pylebedev pytessel
	conda install -c conda-forge mendeleev

PyPi
----

Open a terminal and run the following command:

.. code:: bash

	pip install pydft pyqint pylebedev pytessel mendeleev

Testing
-------

To test that your installation is working, you can run the following snippet
of code

.. code:: python

	import pydft
	print(pydft.__version__)

The version number should be returned.

Simple calculation
------------------

To perform a more involved calculation, one can run the following small
example code:

.. code:: python

	from pydft import MoleculeBuilder,DFT

	co = MoleculeBuilder().get_molecule("CO")
	dft = DFT(co, basis='sto3g')
	en = dft.scf(1e-4)
	print("Total electronic energy: %f Ht" % en)

which should return the following result:

.. code::

	Total electronic energy: -111.147096 Ht

Troubleshooting
---------------

The Anaconda packaging system can sometimes be quite finicky and 
packages conflict with each other. A common scenario wherein this happens is
when you have installed a combination of Conda and PyPi packages.
A way to work around this issue is to create a separate environment and only 
use that environment for the electronic resources associated with this project.

.. note::
	
	The result of each command is explicitly shown. Please do not let this
	overwhelm you. You do not have to explicitly check every line. These
	results are merely shown so that you know what to expect.

To create the new environment (called :code:`pydft-env`), run::

    conda create -n pydft-env python=3.11

which will give the following result::

	## Package Plan ##

	  environment location: C:\Users\iawfi\anaconda3\envs\pydft-env

	  added / updated specs:
	    - python=3.11


	The following packages will be downloaded:

	    package                    |            build
	    ---------------------------|-----------------
	    bzip2-1.0.8                |       hcfcfb64_5         122 KB  conda-forge
	    ca-certificates-2023.11.17 |       h56e8100_0         151 KB  conda-forge
	    libexpat-2.5.0             |       h63175ca_1         135 KB  conda-forge
	    libffi-3.4.2               |       h8ffe710_5          41 KB  conda-forge
	    libsqlite-3.44.2           |       hcfcfb64_0         833 KB  conda-forge
	    libzlib-1.2.13             |       hcfcfb64_5          54 KB  conda-forge
	    openssl-3.2.0              |       hcfcfb64_0         7.8 MB  conda-forge
	    pip-23.3.1                 |     pyhd8ed1ab_0         1.3 MB  conda-forge
	    python-3.11.6              |h2628c8c_0_cpython        17.3 MB  conda-forge
	    setuptools-68.2.2          |     pyhd8ed1ab_0         454 KB  conda-forge
	    tk-8.6.13                  |       h5226925_1         3.3 MB  conda-forge
	    tzdata-2023c               |       h71feb2d_0         115 KB  conda-forge
	    ucrt-10.0.22621.0          |       h57928b3_0         1.2 MB  conda-forge
	    vc-14.3                    |      h64f974e_17          17 KB  conda-forge
	    vc14_runtime-14.36.32532   |      hdcecf7f_17         722 KB  conda-forge
	    vs2015_runtime-14.36.32532 |      h05e6639_17          17 KB  conda-forge
	    wheel-0.41.3               |     pyhd8ed1ab_0          57 KB  conda-forge
	    xz-5.2.6                   |       h8d14728_0         213 KB  conda-forge
	    ------------------------------------------------------------
	                                           Total:        33.8 MB

	The following NEW packages will be INSTALLED:

	  bzip2              conda-forge/win-64::bzip2-1.0.8-hcfcfb64_5
	  ca-certificates    conda-forge/win-64::ca-certificates-2023.11.17-h56e8100_0
	  libexpat           conda-forge/win-64::libexpat-2.5.0-h63175ca_1
	  libffi             conda-forge/win-64::libffi-3.4.2-h8ffe710_5
	  libsqlite          conda-forge/win-64::libsqlite-3.44.2-hcfcfb64_0
	  libzlib            conda-forge/win-64::libzlib-1.2.13-hcfcfb64_5
	  openssl            conda-forge/win-64::openssl-3.2.0-hcfcfb64_0
	  pip                conda-forge/noarch::pip-23.3.1-pyhd8ed1ab_0
	  python             conda-forge/win-64::python-3.11.6-h2628c8c_0_cpython
	  setuptools         conda-forge/noarch::setuptools-68.2.2-pyhd8ed1ab_0
	  tk                 conda-forge/win-64::tk-8.6.13-h5226925_1
	  tzdata             conda-forge/noarch::tzdata-2023c-h71feb2d_0
	  ucrt               conda-forge/win-64::ucrt-10.0.22621.0-h57928b3_0
	  vc                 conda-forge/win-64::vc-14.3-h64f974e_17
	  vc14_runtime       conda-forge/win-64::vc14_runtime-14.36.32532-hdcecf7f_17
	  vs2015_runtime     conda-forge/win-64::vs2015_runtime-14.36.32532-h05e6639_17
	  wheel              conda-forge/noarch::wheel-0.41.3-pyhd8ed1ab_0
	  xz                 conda-forge/win-64::xz-5.2.6-h8d14728_0


	Proceed ([y]/n)? y


	Downloading and Extracting Packages

	Preparing transaction: done
	Verifying transaction: done
	Executing transaction: done
	#
	# To activate this environment, use
	#
	#     $ conda activate pydft-env
	#
	# To deactivate an active environment, use
	#
	#     $ conda deactivate

Next, we will active the environment::

    conda activate pydft-env

You will see that your command line now starts with :code:`(pydft-env)` instead
of :code:`(base)`.

We can now install the required packages into the environment::

    conda install -c ifilot pyqint pylebedev pydft pytessel

You will see a response similar to the one as seen below::

	## Package Plan ##

	  environment location: C:\Users\iawfi\anaconda3\envs\pydft-env

	  added / updated specs:
	    - pydft
	    - pylebedev
	    - pyqint
	    - pytessel


	The following packages will be downloaded:

	    package                    |            build
	    ---------------------------|-----------------
	    colorama-0.4.6             |     pyhd8ed1ab_0          25 KB  conda-forge
	    intel-openmp-2023.2.0      |   h57928b3_50497         2.4 MB  conda-forge
	    libblas-3.9.0              |     20_win64_mkl         4.8 MB  conda-forge
	    libcblas-3.9.0             |     20_win64_mkl         4.8 MB  conda-forge
	    libhwloc-2.9.3             |default_haede6df_1009         2.5 MB  conda-forge
	    libiconv-1.17              |       h8ffe710_0         698 KB  conda-forge
	    liblapack-3.9.0            |     20_win64_mkl         4.8 MB  conda-forge
	    libxml2-2.11.6             |       hc3477c8_0         1.6 MB  conda-forge
	    mkl-2023.2.0               |   h6a75c08_50497       138.0 MB  conda-forge
	    numpy-1.26.2               |  py311h0b4df5a_0         6.8 MB  conda-forge
	    pthreads-win32-2.9.1       |       hfa6e2cd_3         141 KB  conda-forge
	    pydft-0.2.4.1              |     pyh4f56d60_0          52 KB  ifilot
	    pyqint-0.14.0.1            |  py311hcfd9ee6_0         286 KB  ifilot
	    pytessel-1.1.0             |  py311h1d48e73_0          55 KB  ifilot
	    python_abi-3.11            |          4_cp311           7 KB  conda-forge
	    scipy-1.11.4               |  py311h0b4df5a_0        14.2 MB  conda-forge
	    tbb-2021.10.0              |       h91493d7_2         153 KB  conda-forge
	    tqdm-4.66.1                |     pyhd8ed1ab_0          87 KB  conda-forge
	    ------------------------------------------------------------
	                                           Total:       181.1 MB

	The following NEW packages will be INSTALLED:

	  colorama           conda-forge/noarch::colorama-0.4.6-pyhd8ed1ab_0
	  intel-openmp       conda-forge/win-64::intel-openmp-2023.2.0-h57928b3_50497
	  libblas            conda-forge/win-64::libblas-3.9.0-20_win64_mkl
	  libcblas           conda-forge/win-64::libcblas-3.9.0-20_win64_mkl
	  libhwloc           conda-forge/win-64::libhwloc-2.9.3-default_haede6df_1009
	  libiconv           conda-forge/win-64::libiconv-1.17-h8ffe710_0
	  liblapack          conda-forge/win-64::liblapack-3.9.0-20_win64_mkl
	  libxml2            conda-forge/win-64::libxml2-2.11.6-hc3477c8_0
	  mkl                conda-forge/win-64::mkl-2023.2.0-h6a75c08_50497
	  numpy              conda-forge/win-64::numpy-1.26.2-py311h0b4df5a_0
	  pthreads-win32     conda-forge/win-64::pthreads-win32-2.9.1-hfa6e2cd_3
	  pydft              ifilot/noarch::pydft-0.2.4.1-pyh4f56d60_0
	  pylebedev          ifilot/noarch::pylebedev-1.0.0-pyh1d129d4_0
	  pyqint             ifilot/win-64::pyqint-0.14.0.1-py311hcfd9ee6_0
	  pytessel           ifilot/win-64::pytessel-1.1.0-py311h1d48e73_0
	  python_abi         conda-forge/win-64::python_abi-3.11-4_cp311
	  scipy              conda-forge/win-64::scipy-1.11.4-py311h0b4df5a_0
	  tbb                conda-forge/win-64::tbb-2021.10.0-h91493d7_2
	  tqdm               conda-forge/noarch::tqdm-4.66.1-pyhd8ed1ab_0


	Proceed ([y]/n)? y


	Downloading and Extracting Packages

	Preparing transaction: done
	Verifying transaction: done
	Executing transaction: done

Finally, you can install the IDE Spyder using::

    conda install spyder matplotlib scipy pandas openpyxl mendeleev

This might take a while (the environment needs to resolve all dependencies), 
but you should see a response similar to the one below::

	## Package Plan ##

	  environment location: C:\Users\iawfi\anaconda3\envs\pydft-env

	  added / updated specs:
	    - matplotlib
	    - openpyxl
	    - pandas
	    - scipy
	    - spyder


	The following packages will be downloaded:

	    package                    |            build
	    ---------------------------|-----------------
	    alabaster-0.7.13           |     pyhd8ed1ab_0          18 KB  conda-forge
	    arrow-1.3.0                |     pyhd8ed1ab_0          98 KB  conda-forge
	    astroid-3.0.1              |  py311h1ea47a8_0         498 KB  conda-forge
	    asttokens-2.4.1            |     pyhd8ed1ab_0          28 KB  conda-forge
	    atomicwrites-1.4.1         |     pyhd8ed1ab_0          12 KB  conda-forge
	    attrs-23.1.0               |     pyh71513ae_1          54 KB  conda-forge
	    autopep8-2.0.4             |     pyhd8ed1ab_0          45 KB  conda-forge
	    babel-2.13.1               |     pyhd8ed1ab_0         6.6 MB  conda-forge
	    bcrypt-4.0.1               |  py311hc37eb10_1         144 KB  conda-forge
	    beautifulsoup4-4.12.2      |     pyha770c72_0         112 KB  conda-forge
	    binaryornot-0.4.4          |             py_1         370 KB  conda-forge
	    black-23.10.1              |  py311h1ea47a8_0         369 KB  conda-forge
	    bleach-6.1.0               |     pyhd8ed1ab_0         128 KB  conda-forge
	    brotli-1.1.0               |       hcfcfb64_1          19 KB  conda-forge
	    brotli-bin-1.1.0           |       hcfcfb64_1          20 KB  conda-forge
	    brotli-python-1.1.0        |  py311h12c1d0e_1         315 KB  conda-forge
	    certifi-2023.11.17         |     pyhd8ed1ab_0         155 KB  conda-forge
	    cffi-1.16.0                |  py311ha68e1ae_0         290 KB  conda-forge
	    chardet-5.2.0              |  py311h1ea47a8_1         278 KB  conda-forge
	    charset-normalizer-3.3.2   |     pyhd8ed1ab_0          46 KB  conda-forge
	    click-8.1.7                | win_pyh7428d3b_0          83 KB  conda-forge
	    cloudpickle-3.0.0          |     pyhd8ed1ab_0          24 KB  conda-forge
	    comm-0.1.4                 |     pyhd8ed1ab_0          11 KB  conda-forge
	    contourpy-1.2.0            |  py311h005e61a_0         201 KB  conda-forge
	    cookiecutter-2.5.0         |     pyhca7485f_0          97 KB  conda-forge
	    cryptography-41.0.5        |  py311h28e9c30_0         1.1 MB  conda-forge
	    cycler-0.12.1              |     pyhd8ed1ab_0          13 KB  conda-forge
	    debugpy-1.8.0              |  py311h12c1d0e_1         3.7 MB  conda-forge
	    decorator-5.1.1            |     pyhd8ed1ab_0          12 KB  conda-forge
	    defusedxml-0.7.1           |     pyhd8ed1ab_0          23 KB  conda-forge
	    diff-match-patch-20230430  |     pyhd8ed1ab_0          40 KB  conda-forge
	    dill-0.3.7                 |     pyhd8ed1ab_0          86 KB  conda-forge
	    docstring-to-markdown-0.13 |     pyhd8ed1ab_0          31 KB  conda-forge
	    docutils-0.20.1            |  py311h1ea47a8_2         950 KB  conda-forge
	    entrypoints-0.4            |     pyhd8ed1ab_0           9 KB  conda-forge
	    et_xmlfile-1.1.0           |     pyhd8ed1ab_0          10 KB  conda-forge
	    exceptiongroup-1.2.0       |     pyhd8ed1ab_0          20 KB  conda-forge
	    executing-2.0.1            |     pyhd8ed1ab_0          27 KB  conda-forge
	    flake8-6.1.0               |     pyhd8ed1ab_0         109 KB  conda-forge
	    fonttools-4.45.1           |  py311ha68e1ae_0         2.3 MB  conda-forge
	    freetype-2.12.1            |       hdaf720e_2         498 KB  conda-forge
	    gettext-0.21.1             |       h5728263_0         5.3 MB  conda-forge
	    glib-2.78.1                |       h12be248_1         495 KB  conda-forge
	    glib-tools-2.78.1          |       h12be248_1         141 KB  conda-forge
	    gst-plugins-base-1.22.7    |       h001b923_0         1.9 MB  conda-forge
	    gstreamer-1.22.7           |       hb4038d2_0         1.8 MB  conda-forge
	    icu-72.1                   |       h63175ca_0        12.6 MB  conda-forge
	    idna-3.5                   |     pyhd8ed1ab_0          48 KB  conda-forge
	    imagesize-1.4.1            |     pyhd8ed1ab_0          10 KB  conda-forge
	    importlib-metadata-6.8.0   |     pyha770c72_0          25 KB  conda-forge
	    importlib_metadata-6.8.0   |       hd8ed1ab_0           9 KB  conda-forge
	    importlib_resources-6.1.1  |     pyhd8ed1ab_0          29 KB  conda-forge
	    inflection-0.5.1           |     pyh9f0ad1d_0           9 KB  conda-forge
	    intervaltree-3.1.0         |     pyhd8ed1ab_1          27 KB  conda-forge
	    ipykernel-6.26.0           |     pyha63f2e9_0         114 KB  conda-forge
	    ipython-8.18.0             |     pyh5737063_0         577 KB  conda-forge
	    isort-5.12.0               |     pyhd8ed1ab_1          72 KB  conda-forge
	    jaraco.classes-3.3.0       |     pyhd8ed1ab_0          11 KB  conda-forge
	    jedi-0.19.1                |     pyhd8ed1ab_0         822 KB  conda-forge
	    jellyfish-1.0.3            |  py311hc37eb10_0         191 KB  conda-forge
	    jinja2-3.1.2               |     pyhd8ed1ab_1          99 KB  conda-forge
	    jsonschema-4.20.0          |     pyhd8ed1ab_0          70 KB  conda-forge
	    jsonschema-specifications-2023.11.1|     pyhd8ed1ab_0          15 KB  conda-forge
	    jupyter_client-8.6.0       |     pyhd8ed1ab_0         103 KB  conda-forge
	    jupyter_core-5.5.0         |  py311h1ea47a8_0         108 KB  conda-forge
	    jupyterlab_pygments-0.3.0  |     pyhd8ed1ab_0          18 KB  conda-forge
	    keyring-24.3.0             |  py311h1ea47a8_0          92 KB  conda-forge
	    kiwisolver-1.4.5           |  py311h005e61a_1          55 KB  conda-forge
	    krb5-1.20.1                |       heb0366b_0         701 KB  conda-forge
	    lcms2-2.15                 |       he9d350c_2         486 KB  conda-forge
	    lerc-4.0.0                 |       h63175ca_0         190 KB  conda-forge
	    libbrotlicommon-1.1.0      |       hcfcfb64_1          69 KB  conda-forge
	    libbrotlidec-1.1.0         |       hcfcfb64_1          32 KB  conda-forge
	    libbrotlienc-1.1.0         |       hcfcfb64_1         241 KB  conda-forge
	    libclang-16.0.6            |default_heb8d277_2          35 KB  conda-forge
	    libclang13-16.0.6          |default_hc80b9e7_2        22.1 MB  conda-forge
	    libdeflate-1.19            |       hcfcfb64_0         150 KB  conda-forge
	    libglib-2.78.1             |       h16e383f_1         2.5 MB  conda-forge
	    libjpeg-turbo-2.1.5.1      |       hcfcfb64_1         672 KB  conda-forge
	    libogg-1.3.4               |       h8ffe710_1          34 KB  conda-forge
	    libpng-1.6.39              |       h19919ed_0         336 KB  conda-forge
	    libsodium-1.0.18           |       h8d14728_1         697 KB  conda-forge
	    libspatialindex-1.9.3      |       h39d44d4_4         437 KB  conda-forge
	    libtiff-4.6.0              |       h4554b19_1         766 KB  conda-forge
	    libvorbis-1.3.7            |       h0e60522_0         267 KB  conda-forge
	    libwebp-1.3.2              |       hcfcfb64_1          69 KB  conda-forge
	    libwebp-base-1.3.2         |       hcfcfb64_0         263 KB  conda-forge
	    libxcb-1.15                |       hcd874cb_0         947 KB  conda-forge
	    m2w64-gcc-libgfortran-5.3.0|                6         342 KB  conda-forge
	    m2w64-gcc-libs-5.3.0       |                7         520 KB  conda-forge
	    m2w64-gcc-libs-core-5.3.0  |                7         214 KB  conda-forge
	    m2w64-gmp-6.1.0            |                2         726 KB  conda-forge
	    m2w64-libwinpthread-git-5.0.0.4634.697f757|                2          31 KB  conda-forge
	    markdown-it-py-3.0.0       |     pyhd8ed1ab_0          63 KB  conda-forge
	    markupsafe-2.1.3           |  py311ha68e1ae_1          29 KB  conda-forge
	    matplotlib-3.8.2           |  py311h1ea47a8_0           9 KB  conda-forge
	    matplotlib-base-3.8.2      |  py311h6e989c2_0         7.3 MB  conda-forge
	    matplotlib-inline-0.1.6    |     pyhd8ed1ab_0          12 KB  conda-forge
	    mccabe-0.7.0               |     pyhd8ed1ab_0          11 KB  conda-forge
	    mdurl-0.1.0                |     pyhd8ed1ab_0          13 KB  conda-forge
	    mistune-3.0.2              |     pyhd8ed1ab_0          64 KB  conda-forge
	    more-itertools-10.1.0      |     pyhd8ed1ab_0          52 KB  conda-forge
	    msys2-conda-epoch-20160418 |                1           3 KB  conda-forge
	    munkres-1.1.4              |     pyh9f0ad1d_0          12 KB  conda-forge
	    mypy_extensions-1.0.0      |     pyha770c72_0          10 KB  conda-forge
	    nbclient-0.8.0             |     pyhd8ed1ab_0          63 KB  conda-forge
	    nbconvert-7.11.0           |     pyhd8ed1ab_0           8 KB  conda-forge
	    nbconvert-core-7.11.0      |     pyhd8ed1ab_0         183 KB  conda-forge
	    nbconvert-pandoc-7.11.0    |     pyhd8ed1ab_0           7 KB  conda-forge
	    nbformat-5.9.2             |     pyhd8ed1ab_0          98 KB  conda-forge
	    nest-asyncio-1.5.8         |     pyhd8ed1ab_0          11 KB  conda-forge
	    numpydoc-1.5.0             |     pyhd8ed1ab_0          46 KB  conda-forge
	    openjpeg-2.5.0             |       h3d672ee_3         231 KB  conda-forge
	    openpyxl-3.1.2             |  py311ha68e1ae_1         635 KB  conda-forge
	    packaging-23.2             |     pyhd8ed1ab_0          48 KB  conda-forge
	    pandas-2.1.3               |  py311hf63dbb6_0        13.2 MB  conda-forge
	    pandoc-3.1.3               |       h57928b3_0        17.8 MB  conda-forge
	    pandocfilters-1.5.0        |     pyhd8ed1ab_0          11 KB  conda-forge
	    paramiko-3.3.1             |     pyhd8ed1ab_0         155 KB  conda-forge
	    parso-0.8.3                |     pyhd8ed1ab_0          69 KB  conda-forge
	    pathspec-0.11.2            |     pyhd8ed1ab_0          38 KB  conda-forge
	    pcre2-10.42                |       h17e33f8_0         860 KB  conda-forge
	    pexpect-4.8.0              |     pyh1a96a4e_2          48 KB  conda-forge
	    pickleshare-0.7.5          |          py_1003           9 KB  conda-forge
	    pillow-10.0.1              |  py311hd926f49_1        44.7 MB  conda-forge
	    pkgutil-resolve-name-1.3.10|     pyhd8ed1ab_1          11 KB  conda-forge
	    platformdirs-4.0.0         |     pyhd8ed1ab_0          19 KB  conda-forge
	    pluggy-1.3.0               |     pyhd8ed1ab_0          22 KB  conda-forge
	    ply-3.11                   |             py_1          44 KB  conda-forge
	    prompt-toolkit-3.0.41      |     pyha770c72_0         264 KB  conda-forge
	    prompt_toolkit-3.0.41      |       hd8ed1ab_0           7 KB  conda-forge
	    psutil-5.9.5               |  py311ha68e1ae_1         504 KB  conda-forge
	    pthread-stubs-0.4          |    hcd874cb_1001           6 KB  conda-forge
	    ptyprocess-0.7.0           |     pyhd3deb0d_0          16 KB  conda-forge
	    pure_eval-0.2.2            |     pyhd8ed1ab_0          14 KB  conda-forge
	    pycodestyle-2.11.1         |     pyhd8ed1ab_0          34 KB  conda-forge
	    pycparser-2.21             |     pyhd8ed1ab_0         100 KB  conda-forge
	    pydocstyle-6.3.0           |     pyhd8ed1ab_0          39 KB  conda-forge
	    pyflakes-3.1.0             |     pyhd8ed1ab_0          57 KB  conda-forge
	    pygments-2.17.2            |     pyhd8ed1ab_0         840 KB  conda-forge
	    pylint-3.0.2               |     pyhd8ed1ab_0         338 KB  conda-forge
	    pylint-venv-3.0.3          |     pyhd8ed1ab_0          11 KB  conda-forge
	    pyls-spyder-0.4.0          |     pyhd8ed1ab_0          10 KB  conda-forge
	    pynacl-1.5.0               |  py311hd53affc_3         1.2 MB  conda-forge
	    pyparsing-3.1.1            |     pyhd8ed1ab_0          87 KB  conda-forge
	    pyqt-5.15.9                |  py311h125bc19_5         3.7 MB  conda-forge
	    pyqt5-sip-12.12.2          |  py311h12c1d0e_5          78 KB  conda-forge
	    pyqtwebengine-5.15.9       |  py311h5a77453_5         123 KB  conda-forge
	    pysocks-1.7.1              |     pyh0701188_6          19 KB  conda-forge
	    python-dateutil-2.8.2      |     pyhd8ed1ab_0         240 KB  conda-forge
	    python-fastjsonschema-2.19.0|     pyhd8ed1ab_0         221 KB  conda-forge
	    python-lsp-black-1.3.0     |     pyhd8ed1ab_0          12 KB  conda-forge
	    python-lsp-jsonrpc-1.1.2   |     pyhd8ed1ab_0          14 KB  conda-forge
	    python-lsp-server-1.9.0    |     pyhd8ed1ab_0           7 KB  conda-forge
	    python-lsp-server-base-1.9.0|     pyhd8ed1ab_0          60 KB  conda-forge
	    python-slugify-8.0.1       |     pyhd8ed1ab_2          15 KB  conda-forge
	    python-tzdata-2023.3       |     pyhd8ed1ab_0         140 KB  conda-forge
	    pytoolconfig-1.2.5         |     pyhd8ed1ab_0          21 KB  conda-forge
	    pytz-2023.3.post1          |     pyhd8ed1ab_0         183 KB  conda-forge
	    pywin32-306                |  py311h12c1d0e_2         5.8 MB  conda-forge
	    pywin32-ctypes-0.2.2       |  py311h1ea47a8_1          56 KB  conda-forge
	    pyyaml-6.0.1               |  py311ha68e1ae_1         171 KB  conda-forge
	    pyzmq-25.1.1               |  py311h9250fbb_2         480 KB  conda-forge
	    qdarkstyle-3.2             |     pyhd8ed1ab_0         612 KB  conda-forge
	    qstylizer-0.2.2            |     pyhd8ed1ab_0          17 KB  conda-forge
	    qt-main-5.15.8             |      h2c8576c_12        56.9 MB  conda-forge
	    qt-webengine-5.15.8        |       h5b1ea0b_0        62.1 MB  conda-forge
	    qtawesome-1.2.3            |     pyhd8ed1ab_0         1.4 MB  conda-forge
	    qtconsole-5.5.1            |     pyhd8ed1ab_0           7 KB  conda-forge
	    qtconsole-base-5.5.1       |     pyha770c72_0          98 KB  conda-forge
	    qtpy-2.4.1                 |     pyhd8ed1ab_0          60 KB  conda-forge
	    referencing-0.31.0         |     pyhd8ed1ab_0          37 KB  conda-forge
	    requests-2.31.0            |     pyhd8ed1ab_0          55 KB  conda-forge
	    rich-13.7.0                |     pyhd8ed1ab_0         180 KB  conda-forge
	    rope-1.11.0                |     pyhd8ed1ab_1         145 KB  conda-forge
	    rpds-py-0.13.1             |  py311hc37eb10_0         178 KB  conda-forge
	    rtree-1.1.0                |  py311hcacb13a_0          62 KB  conda-forge
	    sip-6.7.12                 |  py311h12c1d0e_0         581 KB  conda-forge
	    six-1.16.0                 |     pyh6c4a22f_0          14 KB  conda-forge
	    snowballstemmer-2.2.0      |     pyhd8ed1ab_0          57 KB  conda-forge
	    sortedcontainers-2.4.0     |     pyhd8ed1ab_0          26 KB  conda-forge
	    soupsieve-2.5              |     pyhd8ed1ab_1          36 KB  conda-forge
	    sphinx-7.2.6               |     pyhd8ed1ab_0         1.2 MB  conda-forge
	    sphinxcontrib-applehelp-1.0.7|     pyhd8ed1ab_0          29 KB  conda-forge
	    sphinxcontrib-devhelp-1.0.5|     pyhd8ed1ab_0          24 KB  conda-forge
	    sphinxcontrib-htmlhelp-2.0.4|     pyhd8ed1ab_0          32 KB  conda-forge
	    sphinxcontrib-jsmath-1.0.1 |     pyhd8ed1ab_0          10 KB  conda-forge
	    sphinxcontrib-qthelp-1.0.6 |     pyhd8ed1ab_0          26 KB  conda-forge
	    sphinxcontrib-serializinghtml-1.1.9|     pyhd8ed1ab_0          28 KB  conda-forge
	    spyder-5.5.0               |  py311h1ea47a8_3        11.2 MB  conda-forge
	    spyder-kernels-2.5.0       | win_pyh7428d3b_0          80 KB  conda-forge
	    stack_data-0.6.2           |     pyhd8ed1ab_0          26 KB  conda-forge
	    text-unidecode-1.3         |     pyhd8ed1ab_1          64 KB  conda-forge
	    textdistance-4.5.0         |     pyhd8ed1ab_0          28 KB  conda-forge
	    three-merge-0.1.1          |     pyh9f0ad1d_0           8 KB  conda-forge
	    tinycss2-1.2.1             |     pyhd8ed1ab_0          23 KB  conda-forge
	    toml-0.10.2                |     pyhd8ed1ab_0          18 KB  conda-forge
	    tomli-2.0.1                |     pyhd8ed1ab_0          16 KB  conda-forge
	    tomlkit-0.12.3             |     pyha770c72_0          36 KB  conda-forge
	    tornado-6.3.3              |  py311ha68e1ae_1         826 KB  conda-forge
	    traitlets-5.13.0           |     pyhd8ed1ab_0         107 KB  conda-forge
	    types-python-dateutil-2.8.19.14|     pyhd8ed1ab_0          21 KB  conda-forge
	    typing-extensions-4.8.0    |       hd8ed1ab_0          10 KB  conda-forge
	    typing_extensions-4.8.0    |     pyha770c72_0          34 KB  conda-forge
	    ujson-5.8.0                |  py311h12c1d0e_0          47 KB  conda-forge
	    urllib3-2.1.0              |     pyhd8ed1ab_0          83 KB  conda-forge
	    watchdog-3.0.0             |  py311h1ea47a8_1         151 KB  conda-forge
	    wcwidth-0.2.12             |     pyhd8ed1ab_0          32 KB  conda-forge
	    webencodings-0.5.1         |     pyhd8ed1ab_2          15 KB  conda-forge
	    whatthepatch-1.0.5         |     pyhd8ed1ab_0          17 KB  conda-forge
	    win_inet_pton-1.1.0        |     pyhd8ed1ab_6           8 KB  conda-forge
	    xorg-libxau-1.0.11         |       hcd874cb_0          50 KB  conda-forge
	    xorg-libxdmcp-1.1.3        |       hcd874cb_0          66 KB  conda-forge
	    yaml-0.2.5                 |       h8ffe710_2          62 KB  conda-forge
	    yapf-0.40.1                |     pyhd8ed1ab_0         172 KB  conda-forge
	    zeromq-4.3.5               |       h63175ca_0         4.0 MB  conda-forge
	    zipp-3.17.0                |     pyhd8ed1ab_0          19 KB  conda-forge
	    zstd-1.5.5                 |       h12be248_0         335 KB  conda-forge
	    ------------------------------------------------------------
	                                           Total:       318.3 MB

	The following NEW packages will be INSTALLED:

	  alabaster          conda-forge/noarch::alabaster-0.7.13-pyhd8ed1ab_0
	  arrow              conda-forge/noarch::arrow-1.3.0-pyhd8ed1ab_0
	  astroid            conda-forge/win-64::astroid-3.0.1-py311h1ea47a8_0
	  asttokens          conda-forge/noarch::asttokens-2.4.1-pyhd8ed1ab_0
	  atomicwrites       conda-forge/noarch::atomicwrites-1.4.1-pyhd8ed1ab_0
	  attrs              conda-forge/noarch::attrs-23.1.0-pyh71513ae_1
	  autopep8           conda-forge/noarch::autopep8-2.0.4-pyhd8ed1ab_0
	  babel              conda-forge/noarch::babel-2.13.1-pyhd8ed1ab_0
	  bcrypt             conda-forge/win-64::bcrypt-4.0.1-py311hc37eb10_1
	  beautifulsoup4     conda-forge/noarch::beautifulsoup4-4.12.2-pyha770c72_0
	  binaryornot        conda-forge/noarch::binaryornot-0.4.4-py_1
	  black              conda-forge/win-64::black-23.10.1-py311h1ea47a8_0
	  bleach             conda-forge/noarch::bleach-6.1.0-pyhd8ed1ab_0
	  brotli             conda-forge/win-64::brotli-1.1.0-hcfcfb64_1
	  brotli-bin         conda-forge/win-64::brotli-bin-1.1.0-hcfcfb64_1
	  brotli-python      conda-forge/win-64::brotli-python-1.1.0-py311h12c1d0e_1
	  certifi            conda-forge/noarch::certifi-2023.11.17-pyhd8ed1ab_0
	  cffi               conda-forge/win-64::cffi-1.16.0-py311ha68e1ae_0
	  chardet            conda-forge/win-64::chardet-5.2.0-py311h1ea47a8_1
	  charset-normalizer conda-forge/noarch::charset-normalizer-3.3.2-pyhd8ed1ab_0
	  click              conda-forge/noarch::click-8.1.7-win_pyh7428d3b_0
	  cloudpickle        conda-forge/noarch::cloudpickle-3.0.0-pyhd8ed1ab_0
	  comm               conda-forge/noarch::comm-0.1.4-pyhd8ed1ab_0
	  contourpy          conda-forge/win-64::contourpy-1.2.0-py311h005e61a_0
	  cookiecutter       conda-forge/noarch::cookiecutter-2.5.0-pyhca7485f_0
	  cryptography       conda-forge/win-64::cryptography-41.0.5-py311h28e9c30_0
	  cycler             conda-forge/noarch::cycler-0.12.1-pyhd8ed1ab_0
	  debugpy            conda-forge/win-64::debugpy-1.8.0-py311h12c1d0e_1
	  decorator          conda-forge/noarch::decorator-5.1.1-pyhd8ed1ab_0
	  defusedxml         conda-forge/noarch::defusedxml-0.7.1-pyhd8ed1ab_0
	  diff-match-patch   conda-forge/noarch::diff-match-patch-20230430-pyhd8ed1ab_0
	  dill               conda-forge/noarch::dill-0.3.7-pyhd8ed1ab_0
	  docstring-to-mark~ conda-forge/noarch::docstring-to-markdown-0.13-pyhd8ed1ab_0
	  docutils           conda-forge/win-64::docutils-0.20.1-py311h1ea47a8_2
	  entrypoints        conda-forge/noarch::entrypoints-0.4-pyhd8ed1ab_0
	  et_xmlfile         conda-forge/noarch::et_xmlfile-1.1.0-pyhd8ed1ab_0
	  exceptiongroup     conda-forge/noarch::exceptiongroup-1.2.0-pyhd8ed1ab_0
	  executing          conda-forge/noarch::executing-2.0.1-pyhd8ed1ab_0
	  flake8             conda-forge/noarch::flake8-6.1.0-pyhd8ed1ab_0
	  fonttools          conda-forge/win-64::fonttools-4.45.1-py311ha68e1ae_0
	  freetype           conda-forge/win-64::freetype-2.12.1-hdaf720e_2
	  gettext            conda-forge/win-64::gettext-0.21.1-h5728263_0
	  glib               conda-forge/win-64::glib-2.78.1-h12be248_1
	  glib-tools         conda-forge/win-64::glib-tools-2.78.1-h12be248_1
	  gst-plugins-base   conda-forge/win-64::gst-plugins-base-1.22.7-h001b923_0
	  gstreamer          conda-forge/win-64::gstreamer-1.22.7-hb4038d2_0
	  icu                conda-forge/win-64::icu-72.1-h63175ca_0
	  idna               conda-forge/noarch::idna-3.5-pyhd8ed1ab_0
	  imagesize          conda-forge/noarch::imagesize-1.4.1-pyhd8ed1ab_0
	  importlib-metadata conda-forge/noarch::importlib-metadata-6.8.0-pyha770c72_0
	  importlib_metadata conda-forge/noarch::importlib_metadata-6.8.0-hd8ed1ab_0
	  importlib_resourc~ conda-forge/noarch::importlib_resources-6.1.1-pyhd8ed1ab_0
	  inflection         conda-forge/noarch::inflection-0.5.1-pyh9f0ad1d_0
	  intervaltree       conda-forge/noarch::intervaltree-3.1.0-pyhd8ed1ab_1
	  ipykernel          conda-forge/noarch::ipykernel-6.26.0-pyha63f2e9_0
	  ipython            conda-forge/noarch::ipython-8.18.0-pyh5737063_0
	  isort              conda-forge/noarch::isort-5.12.0-pyhd8ed1ab_1
	  jaraco.classes     conda-forge/noarch::jaraco.classes-3.3.0-pyhd8ed1ab_0
	  jedi               conda-forge/noarch::jedi-0.19.1-pyhd8ed1ab_0
	  jellyfish          conda-forge/win-64::jellyfish-1.0.3-py311hc37eb10_0
	  jinja2             conda-forge/noarch::jinja2-3.1.2-pyhd8ed1ab_1
	  jsonschema         conda-forge/noarch::jsonschema-4.20.0-pyhd8ed1ab_0
	  jsonschema-specif~ conda-forge/noarch::jsonschema-specifications-2023.11.1-pyhd8ed1ab_0
	  jupyter_client     conda-forge/noarch::jupyter_client-8.6.0-pyhd8ed1ab_0
	  jupyter_core       conda-forge/win-64::jupyter_core-5.5.0-py311h1ea47a8_0
	  jupyterlab_pygmen~ conda-forge/noarch::jupyterlab_pygments-0.3.0-pyhd8ed1ab_0
	  keyring            conda-forge/win-64::keyring-24.3.0-py311h1ea47a8_0
	  kiwisolver         conda-forge/win-64::kiwisolver-1.4.5-py311h005e61a_1
	  krb5               conda-forge/win-64::krb5-1.20.1-heb0366b_0
	  lcms2              conda-forge/win-64::lcms2-2.15-he9d350c_2
	  lerc               conda-forge/win-64::lerc-4.0.0-h63175ca_0
	  libbrotlicommon    conda-forge/win-64::libbrotlicommon-1.1.0-hcfcfb64_1
	  libbrotlidec       conda-forge/win-64::libbrotlidec-1.1.0-hcfcfb64_1
	  libbrotlienc       conda-forge/win-64::libbrotlienc-1.1.0-hcfcfb64_1
	  libclang           conda-forge/win-64::libclang-16.0.6-default_heb8d277_2
	  libclang13         conda-forge/win-64::libclang13-16.0.6-default_hc80b9e7_2
	  libdeflate         conda-forge/win-64::libdeflate-1.19-hcfcfb64_0
	  libglib            conda-forge/win-64::libglib-2.78.1-h16e383f_1
	  libjpeg-turbo      conda-forge/win-64::libjpeg-turbo-2.1.5.1-hcfcfb64_1
	  libogg             conda-forge/win-64::libogg-1.3.4-h8ffe710_1
	  libpng             conda-forge/win-64::libpng-1.6.39-h19919ed_0
	  libsodium          conda-forge/win-64::libsodium-1.0.18-h8d14728_1
	  libspatialindex    conda-forge/win-64::libspatialindex-1.9.3-h39d44d4_4
	  libtiff            conda-forge/win-64::libtiff-4.6.0-h4554b19_1
	  libvorbis          conda-forge/win-64::libvorbis-1.3.7-h0e60522_0
	  libwebp            conda-forge/win-64::libwebp-1.3.2-hcfcfb64_1
	  libwebp-base       conda-forge/win-64::libwebp-base-1.3.2-hcfcfb64_0
	  libxcb             conda-forge/win-64::libxcb-1.15-hcd874cb_0
	  m2w64-gcc-libgfor~ conda-forge/win-64::m2w64-gcc-libgfortran-5.3.0-6
	  m2w64-gcc-libs     conda-forge/win-64::m2w64-gcc-libs-5.3.0-7
	  m2w64-gcc-libs-co~ conda-forge/win-64::m2w64-gcc-libs-core-5.3.0-7
	  m2w64-gmp          conda-forge/win-64::m2w64-gmp-6.1.0-2
	  m2w64-libwinpthre~ conda-forge/win-64::m2w64-libwinpthread-git-5.0.0.4634.697f757-2
	  markdown-it-py     conda-forge/noarch::markdown-it-py-3.0.0-pyhd8ed1ab_0
	  markupsafe         conda-forge/win-64::markupsafe-2.1.3-py311ha68e1ae_1
	  matplotlib         conda-forge/win-64::matplotlib-3.8.2-py311h1ea47a8_0
	  matplotlib-base    conda-forge/win-64::matplotlib-base-3.8.2-py311h6e989c2_0
	  matplotlib-inline  conda-forge/noarch::matplotlib-inline-0.1.6-pyhd8ed1ab_0
	  mccabe             conda-forge/noarch::mccabe-0.7.0-pyhd8ed1ab_0
	  mdurl              conda-forge/noarch::mdurl-0.1.0-pyhd8ed1ab_0
	  mistune            conda-forge/noarch::mistune-3.0.2-pyhd8ed1ab_0
	  more-itertools     conda-forge/noarch::more-itertools-10.1.0-pyhd8ed1ab_0
	  msys2-conda-epoch  conda-forge/win-64::msys2-conda-epoch-20160418-1
	  munkres            conda-forge/noarch::munkres-1.1.4-pyh9f0ad1d_0
	  mypy_extensions    conda-forge/noarch::mypy_extensions-1.0.0-pyha770c72_0
	  nbclient           conda-forge/noarch::nbclient-0.8.0-pyhd8ed1ab_0
	  nbconvert          conda-forge/noarch::nbconvert-7.11.0-pyhd8ed1ab_0
	  nbconvert-core     conda-forge/noarch::nbconvert-core-7.11.0-pyhd8ed1ab_0
	  nbconvert-pandoc   conda-forge/noarch::nbconvert-pandoc-7.11.0-pyhd8ed1ab_0
	  nbformat           conda-forge/noarch::nbformat-5.9.2-pyhd8ed1ab_0
	  nest-asyncio       conda-forge/noarch::nest-asyncio-1.5.8-pyhd8ed1ab_0
	  numpydoc           conda-forge/noarch::numpydoc-1.5.0-pyhd8ed1ab_0
	  openjpeg           conda-forge/win-64::openjpeg-2.5.0-h3d672ee_3
	  openpyxl           conda-forge/win-64::openpyxl-3.1.2-py311ha68e1ae_1
	  packaging          conda-forge/noarch::packaging-23.2-pyhd8ed1ab_0
	  pandas             conda-forge/win-64::pandas-2.1.3-py311hf63dbb6_0
	  pandoc             conda-forge/win-64::pandoc-3.1.3-h57928b3_0
	  pandocfilters      conda-forge/noarch::pandocfilters-1.5.0-pyhd8ed1ab_0
	  paramiko           conda-forge/noarch::paramiko-3.3.1-pyhd8ed1ab_0
	  parso              conda-forge/noarch::parso-0.8.3-pyhd8ed1ab_0
	  pathspec           conda-forge/noarch::pathspec-0.11.2-pyhd8ed1ab_0
	  pcre2              conda-forge/win-64::pcre2-10.42-h17e33f8_0
	  pexpect            conda-forge/noarch::pexpect-4.8.0-pyh1a96a4e_2
	  pickleshare        conda-forge/noarch::pickleshare-0.7.5-py_1003
	  pillow             conda-forge/win-64::pillow-10.0.1-py311hd926f49_1
	  pkgutil-resolve-n~ conda-forge/noarch::pkgutil-resolve-name-1.3.10-pyhd8ed1ab_1
	  platformdirs       conda-forge/noarch::platformdirs-4.0.0-pyhd8ed1ab_0
	  pluggy             conda-forge/noarch::pluggy-1.3.0-pyhd8ed1ab_0
	  ply                conda-forge/noarch::ply-3.11-py_1
	  prompt-toolkit     conda-forge/noarch::prompt-toolkit-3.0.41-pyha770c72_0
	  prompt_toolkit     conda-forge/noarch::prompt_toolkit-3.0.41-hd8ed1ab_0
	  psutil             conda-forge/win-64::psutil-5.9.5-py311ha68e1ae_1
	  pthread-stubs      conda-forge/win-64::pthread-stubs-0.4-hcd874cb_1001
	  ptyprocess         conda-forge/noarch::ptyprocess-0.7.0-pyhd3deb0d_0
	  pure_eval          conda-forge/noarch::pure_eval-0.2.2-pyhd8ed1ab_0
	  pycodestyle        conda-forge/noarch::pycodestyle-2.11.1-pyhd8ed1ab_0
	  pycparser          conda-forge/noarch::pycparser-2.21-pyhd8ed1ab_0
	  pydocstyle         conda-forge/noarch::pydocstyle-6.3.0-pyhd8ed1ab_0
	  pyflakes           conda-forge/noarch::pyflakes-3.1.0-pyhd8ed1ab_0
	  pygments           conda-forge/noarch::pygments-2.17.2-pyhd8ed1ab_0
	  pylint             conda-forge/noarch::pylint-3.0.2-pyhd8ed1ab_0
	  pylint-venv        conda-forge/noarch::pylint-venv-3.0.3-pyhd8ed1ab_0
	  pyls-spyder        conda-forge/noarch::pyls-spyder-0.4.0-pyhd8ed1ab_0
	  pynacl             conda-forge/win-64::pynacl-1.5.0-py311hd53affc_3
	  pyparsing          conda-forge/noarch::pyparsing-3.1.1-pyhd8ed1ab_0
	  pyqt               conda-forge/win-64::pyqt-5.15.9-py311h125bc19_5
	  pyqt5-sip          conda-forge/win-64::pyqt5-sip-12.12.2-py311h12c1d0e_5
	  pyqtwebengine      conda-forge/win-64::pyqtwebengine-5.15.9-py311h5a77453_5
	  pysocks            conda-forge/noarch::pysocks-1.7.1-pyh0701188_6
	  python-dateutil    conda-forge/noarch::python-dateutil-2.8.2-pyhd8ed1ab_0
	  python-fastjsonsc~ conda-forge/noarch::python-fastjsonschema-2.19.0-pyhd8ed1ab_0
	  python-lsp-black   conda-forge/noarch::python-lsp-black-1.3.0-pyhd8ed1ab_0
	  python-lsp-jsonrpc conda-forge/noarch::python-lsp-jsonrpc-1.1.2-pyhd8ed1ab_0
	  python-lsp-server  conda-forge/noarch::python-lsp-server-1.9.0-pyhd8ed1ab_0
	  python-lsp-server~ conda-forge/noarch::python-lsp-server-base-1.9.0-pyhd8ed1ab_0
	  python-slugify     conda-forge/noarch::python-slugify-8.0.1-pyhd8ed1ab_2
	  python-tzdata      conda-forge/noarch::python-tzdata-2023.3-pyhd8ed1ab_0
	  pytoolconfig       conda-forge/noarch::pytoolconfig-1.2.5-pyhd8ed1ab_0
	  pytz               conda-forge/noarch::pytz-2023.3.post1-pyhd8ed1ab_0
	  pywin32            conda-forge/win-64::pywin32-306-py311h12c1d0e_2
	  pywin32-ctypes     conda-forge/win-64::pywin32-ctypes-0.2.2-py311h1ea47a8_1
	  pyyaml             conda-forge/win-64::pyyaml-6.0.1-py311ha68e1ae_1
	  pyzmq              conda-forge/win-64::pyzmq-25.1.1-py311h9250fbb_2
	  qdarkstyle         conda-forge/noarch::qdarkstyle-3.2-pyhd8ed1ab_0
	  qstylizer          conda-forge/noarch::qstylizer-0.2.2-pyhd8ed1ab_0
	  qt-main            conda-forge/win-64::qt-main-5.15.8-h2c8576c_12
	  qt-webengine       conda-forge/win-64::qt-webengine-5.15.8-h5b1ea0b_0
	  qtawesome          conda-forge/noarch::qtawesome-1.2.3-pyhd8ed1ab_0
	  qtconsole          conda-forge/noarch::qtconsole-5.5.1-pyhd8ed1ab_0
	  qtconsole-base     conda-forge/noarch::qtconsole-base-5.5.1-pyha770c72_0
	  qtpy               conda-forge/noarch::qtpy-2.4.1-pyhd8ed1ab_0
	  referencing        conda-forge/noarch::referencing-0.31.0-pyhd8ed1ab_0
	  requests           conda-forge/noarch::requests-2.31.0-pyhd8ed1ab_0
	  rich               conda-forge/noarch::rich-13.7.0-pyhd8ed1ab_0
	  rope               conda-forge/noarch::rope-1.11.0-pyhd8ed1ab_1
	  rpds-py            conda-forge/win-64::rpds-py-0.13.1-py311hc37eb10_0
	  rtree              conda-forge/win-64::rtree-1.1.0-py311hcacb13a_0
	  sip                conda-forge/win-64::sip-6.7.12-py311h12c1d0e_0
	  six                conda-forge/noarch::six-1.16.0-pyh6c4a22f_0
	  snowballstemmer    conda-forge/noarch::snowballstemmer-2.2.0-pyhd8ed1ab_0
	  sortedcontainers   conda-forge/noarch::sortedcontainers-2.4.0-pyhd8ed1ab_0
	  soupsieve          conda-forge/noarch::soupsieve-2.5-pyhd8ed1ab_1
	  sphinx             conda-forge/noarch::sphinx-7.2.6-pyhd8ed1ab_0
	  sphinxcontrib-app~ conda-forge/noarch::sphinxcontrib-applehelp-1.0.7-pyhd8ed1ab_0
	  sphinxcontrib-dev~ conda-forge/noarch::sphinxcontrib-devhelp-1.0.5-pyhd8ed1ab_0
	  sphinxcontrib-htm~ conda-forge/noarch::sphinxcontrib-htmlhelp-2.0.4-pyhd8ed1ab_0
	  sphinxcontrib-jsm~ conda-forge/noarch::sphinxcontrib-jsmath-1.0.1-pyhd8ed1ab_0
	  sphinxcontrib-qth~ conda-forge/noarch::sphinxcontrib-qthelp-1.0.6-pyhd8ed1ab_0
	  sphinxcontrib-ser~ conda-forge/noarch::sphinxcontrib-serializinghtml-1.1.9-pyhd8ed1ab_0
	  spyder             conda-forge/win-64::spyder-5.5.0-py311h1ea47a8_3
	  spyder-kernels     conda-forge/noarch::spyder-kernels-2.5.0-win_pyh7428d3b_0
	  stack_data         conda-forge/noarch::stack_data-0.6.2-pyhd8ed1ab_0
	  text-unidecode     conda-forge/noarch::text-unidecode-1.3-pyhd8ed1ab_1
	  textdistance       conda-forge/noarch::textdistance-4.5.0-pyhd8ed1ab_0
	  three-merge        conda-forge/noarch::three-merge-0.1.1-pyh9f0ad1d_0
	  tinycss2           conda-forge/noarch::tinycss2-1.2.1-pyhd8ed1ab_0
	  toml               conda-forge/noarch::toml-0.10.2-pyhd8ed1ab_0
	  tomli              conda-forge/noarch::tomli-2.0.1-pyhd8ed1ab_0
	  tomlkit            conda-forge/noarch::tomlkit-0.12.3-pyha770c72_0
	  tornado            conda-forge/win-64::tornado-6.3.3-py311ha68e1ae_1
	  traitlets          conda-forge/noarch::traitlets-5.13.0-pyhd8ed1ab_0
	  types-python-date~ conda-forge/noarch::types-python-dateutil-2.8.19.14-pyhd8ed1ab_0
	  typing-extensions  conda-forge/noarch::typing-extensions-4.8.0-hd8ed1ab_0
	  typing_extensions  conda-forge/noarch::typing_extensions-4.8.0-pyha770c72_0
	  ujson              conda-forge/win-64::ujson-5.8.0-py311h12c1d0e_0
	  urllib3            conda-forge/noarch::urllib3-2.1.0-pyhd8ed1ab_0
	  watchdog           conda-forge/win-64::watchdog-3.0.0-py311h1ea47a8_1
	  wcwidth            conda-forge/noarch::wcwidth-0.2.12-pyhd8ed1ab_0
	  webencodings       conda-forge/noarch::webencodings-0.5.1-pyhd8ed1ab_2
	  whatthepatch       conda-forge/noarch::whatthepatch-1.0.5-pyhd8ed1ab_0
	  win_inet_pton      conda-forge/noarch::win_inet_pton-1.1.0-pyhd8ed1ab_6
	  xorg-libxau        conda-forge/win-64::xorg-libxau-1.0.11-hcd874cb_0
	  xorg-libxdmcp      conda-forge/win-64::xorg-libxdmcp-1.1.3-hcd874cb_0
	  yaml               conda-forge/win-64::yaml-0.2.5-h8ffe710_2
	  yapf               conda-forge/noarch::yapf-0.40.1-pyhd8ed1ab_0
	  zeromq             conda-forge/win-64::zeromq-4.3.5-h63175ca_0
	  zipp               conda-forge/noarch::zipp-3.17.0-pyhd8ed1ab_0
	  zstd               conda-forge/win-64::zstd-1.5.5-h12be248_0


	Proceed ([y]/n)? y


	Downloading and Extracting Packages

	Preparing transaction: done
	Verifying transaction: done
	Executing transaction: done


.. note::

	* When you open the Spyder IDE, make sure you select the right one. If you
	  have multiple installations of Spyder, the specific environment of Spyder
	  will be mentioned between parentheses.
	* You can only open one Spyder window at the same time. Please make sure you
	  have closed all other Spyder windows (specifically those corresponding
	  to a different environment) before opening :code:`Spyder (pydft-env)`.
	* After installing :program:`PyDFT` into its own environment, please revisit
	  the "Simple calculation" section on this page and check that your
	  installation is working.