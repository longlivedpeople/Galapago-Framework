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

The framework combine different datasets by grouping them in blocks and properly scaling them according to the luminosity, pile up and other relevant weights.

The datasets that are used to make the plots (with all the information that is required) are declared in the ```samples.dat``` file. In this file it is specified the block of each dataset, its color, name, label (of the legend), location of the NTuples, cross section and type (DATA or MC). Each user has the option to select which datasets of the ```samples.dat``` file wants to include in the plots in the running files.

<p>The plots creation process is divided in two steps: Loop and plotting:</p>
<ul> 
  <li> In the loop part the framework run over the events of the NTuples, applying the cuts and filling a set of predefined histograms. This step works with instances of the Sample.py, CutManager.py and processHandler.py classes, located in the include/ directory. The filled histograms are stored in a .root file (that the user can specify if desired)</li>
  <li> The plotting step works with instances of both the Sample.py and Canvas.py classes. It accesses the histograms stored in the previously created .root file and make the plots. For this step to succeed the .root file must exist, to the previous step needs to be run at least once time before. The user has the option to skip the loop step if the file already exists.</li>
</ul>

### The samples.dat file

<p> In this <code>.dat</code> file all the details of the datasets reserved for creating the NTuples are specified. The information needed (in the correct order) to be filled is </p>
<ul>
  <li> <strong>Block</strong>: It specifies a <em>block</em> of similar datasets that are drawn jointly e.g. WW, WZ and WZ are drawn together under the block name *Diboson*.</li>
  <li> <strong>Color</strong>: The color of each dataset. </li>
  <li> <strong>Name</strong>: This is the name used to identify the sample. The output histograms used to make the plots will have this identified. It also serves to select in the running file if one dataset is included or not in the plotting.</li>
  <li> <strong>Label</strong>: The text that is put in the legend for each sample (or block). </li>
  <li> <strong>Location of friends</strong>: Absolute path of the .root files of the NTuples. </li>
  <li> <strong>Xsec</strong>: Cross-section of the dataset, used to compute the luminosity weights. Should be filled in fb^{-1}.</li>
  <li> <strong>IsData</strong>: 1 if the sample is a data sample or 0 if it is Monte Carlo.</li>
</ul>

### How to read all (or some) datasets of the samples.dat file

To read the ```samples.dat``` file the user has to define an instance of the class ```Tree``` (defined in ```include/Samples.py```). The command to do so is:

```tree = Sample.Tree('samples.dat', 'Name of tree instance', isData, loopFile)```

where ```'Name of the tree instance'``` is an identifier of the instance created with this dataset, ```isData``` indicates if the dataset is data (1) or Monte Carlo (0) (necessary to apply weights, for example) and ```loopFile``` is the name of the root file where the histograms of this dataset are stored.

<strong>All the datasets described in the ```samples.dat``` file are read.</strong> If the user wants to select some, but not all, he/she can use the function ```selectSamples()``` defined in the ```ìnclude/helper.py``` module. This function takes as an input the ```samples.dat``` file and a list of the datasets identified by the name. It creates a temporary ```.tmp_sampleFileMC.txt``` file with just the selected datasets that is given to the Tree instance instead. The sequence would be: 

```
listOfDatasets = []
listOfDatasets.append('Dataset1')
listOfDatasets.append('Dataset2')
listOfDatasets.append('Dataset2')

tree = Sample.Tree(helper.selectSamples('samples.dat', listOfDatasets, sType), 'Name of tree instance', isData, loopFile)
```


### How to read the trees of the datasets and create the histograms 

Once the TTree's of the datasets are read the histograms can be created by means of the function ```Loop()``` of the class Tree declared in the ```include/Sample.py``` module. All the events of the trees defined in the Tree instance are read and a set of predefined histograms are filled and stored in the ```loopFile``` given as an input when the Tree instance was created first. The command is

```
tree.Loop(luminosity, filename = False, maxNumber = False)
```

where ```luminosity``` stands for the luminosity value to which the histograms are normalized, and ```filename``` and ```maxNumber``` are optional parameters, usually set to ```False```.

The set of predefined histograms is declared in the ```include/processHandler.py``` module. If the user wants to add a histogram to the list, it must be defined here. One histogram is created per dataset (specified in the ```samples.dat``` file) separately.


### Plotting

Once the individual histograms of the datasets are filled and stored in the ```.root``` file, they can be combined and drawn in the final plots. This is done by means of functions defined in the running files. In ```plotDistributions.py``` there are available ```makePlot()```, to combine and compare Monte Carlo backgrounds and signal distributions, and ```makeSensitivity``` to compute the sensitivity S/sqrt(S+B) for a given set of signal samples.


## Instructions to run 

To create the distributions the user has to run the plotDistributions.py file by typing the command:

```python plotDistributions.py --samples samples.dat --output outputHistos.root (--avoidLoop)```

where --samples indicates the samples.dat file that is used, --output indicates the .root file where the histograms are stores and --avoidLoop is an option that allows the user to skip the loop part and plot the distributions directly.
