# ccCluster v. 2.0
copyright 2015-2030

Welcome to ccCluster
Developed at the ESRF by Gianluca Santoni.
Currently developed and maintained by Tiankun Zhou
tiankun.zhou@esrf.fr

This program is used to run hierarchycal cluster analysis on protein diffraction data.
When using, please cite:

Santoni, G., Zander, U., Mueller-Dieckmann, C., Leonard, G. & Popov, A. (2017). J. Appl. Cryst. 50,
https://doi.org/10.1107/S1600576717015229.

## Installation with conda for your own version if you are not at ESRF
You need to have XDS and ccp4 installed first
Also, you need to have conda installed for creating new conda environments

To install ccCluster, you would need to create a virtual environment with libraries installed:
```
conda create --prefix /Path/To/Your/environment/HCA_env -c conda-forge python numpy scipy matplotlib pyqt pyside6 pathos pip cctbx-base
```

After, activate the environment and build the program:
```
conda activate /Path/To/Your/environment/HCA_env

cd ccCluster

python setup.py build 

pip install .
```
Now ccCluster commands will be available every time you activate the virtual environment.
You should have all python libraries that needed for ccCluster.
If not, you can use:

```
conda install libraries_name -c conda forge
```
or
```
pip install libraries_name
```

## Basic Usage
At first, you must run ccCalc to generate the distance files.
ccCalc must receive, the first time you run it for a project, a list of HKL files.
To do this, you can simply call
```
ccCluster-gui
```
Then you can perform the ccCalc, ccCluster and check results easily

Or if you would like to do it in the terminal:
```
ccCalc -f <FILE1>.HKL ... <FILEn>.HKL
```

if no file is specified, it will walk all the subdirectories of the current folder and look for HKL files.
This will produce a file named ccClusterLog.txt
One the run is done, you can open ccCluster.
The most basic run can be launched with
```
ccCluster -p -i ccClusterlog.txt -o output_dir -t threshold -ref reference_HKL -clu cluster1 cluster2 cluster3 ...
```
To run in a fully user-independent mode, you can cal the HCAPipeline command

```
HCAPipeline -f file1 file2 file3 ... -o output_dir -t threshold -ano -ref reference_HKL -clu cluster1 cluster2 cluster3 ...
```

