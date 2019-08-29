# GALAPAGO - Framework

This framework is devoted to make the plots of the IFCA-LongLived Analysis. It works with plain trees, previously created by getting the relevant collections from MINIAOD.

This framework is writen in Python.

## Instructions to install

Initialize a CMSSW working directory where the framework is going to be located (the recommendation is to use the CMSSW_9_4_4 release) and clone the code in some dir created inside the src dir:

```
cmsrel CMSSW_9_4_4

cd CMSSW_9_4_4/src

scram b

cmsenv

mkdir MyAnalysis

cd MyAnalysis

git clone https://github.com/longlivedpeople/Galapago-Framework.git

```

## How this framework works

The running files are located inside this directory, all of them are .py files. All the modules that are required to make the plots are located inside the include/ directory.

The plots combine different datasets by grouping them in blocks and properly scaling them according to the luminosity, pile up and other relevant weights.

The datasets that are used to make the plots (with all the information that is required) are declared in the ```samples.dat``` file. In this file it is specified the block of each dataset, its color, name, label (of the legend), lovation of the NTuples, cross section and type (DATA or MC). Each user has the option to select which datasets of the ```samples.dat``` file wants to include in the plots in the running files.

<p>The plots creation process is divided in two steps: Loop and plotting:</p>
<ul> 
  <li> In the loop part the framework run over the events of the NTuples, applying the cuts and filling a set of predefined histograms. This step works with instances of the Sample.py, CutManager.py and processHandler.py classes, located in the include/ directory. The filled histograms are stored in a .root file (that the user can specify if desired)</li>
  <li> The plotting step works with instances of both the Sample.py and Canvas.py classes. It accesses the histograms stored in the previously created .root file and make the plots. For this step to succeed the .root file must exist, to the previous step needs to be run at least once time before. The user has the option to skip the loop step if the file already exists.</li>
</ul>

## Instructions to run 

To create the distributions the user has to run the plotDistributions.py file by typing the command:

```python plotDistributions.py --samples samples.dat --output outputHistos.root (--avoidLoop)```

where --samples indicates the samples.dat file that is used, --output indicates the .root file where the histograms are stores and --avoidLoop is an option that allows the user to skip the loop part and plot the distributions directly.
