# create_goVirt

`create_goVirt` powers the deployment of the new virtual-site implementation of GōMartini which can be combined with the latest iteration of the [Martini](http://cgmartini.nl/) model, [Martini 3](https://doi.org/10.1038/s41592-021-01098-3). 

## Installation

`create_goVirt` requires python 3.6 or greater. It is distributed via PyPi, and can be installed using the pip command:
``````
pip install XXXXXXX (STILL NEED TO PACKAGE IT.)
``````

This installs the last released version. You can update an existing installation by running `pip install -U XXXXXXX`. In some cases you may want to experiment with running the latest development version. You can install this version with the following command:
``````
pip install git+https://github.com/marrink-lab/vermouth-martinize.git#vermouth
``````

Note that development versions, may contain bugs that cause it to produce incorrect topologies. Check the produced output carefully!

The behavior of the pip command can vary depending of the specificity of your python installation. See the documentation on installing a python package to learn more.

## Basic usage