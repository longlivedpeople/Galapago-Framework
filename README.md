# GALAPAGO - Framework

This framework is devoted to make the plots of the IFCA-LongLived Analysis. It works with plain trees, previously created by getting the relevant collections from MINIAOD.

This framework is writen in Python.

## Instructions to install

Initialize a CMSSW working directory where the framework is going to be located (the recommendation is to use the CMSSW_10_6_20 release) and clone the code in some dir created inside the src dir:

```
cmsrel CMSSW_10_6_20

cd CMSSW_10_6_20/src

scram b

cmsenv

mkdir Analysis

cd Analysis

git clone https://github.com/longlivedpeople/Galapago-Framework.git

```


## How to run

Galapago is designed to work with dedicated Ntuples for displaced vertex analyses. These Ntuples are produced by the analyzer located in:

<p> Galapago works with Ntuples created with the IFCALongLived analyzer, that can be accesed here: https://github.com/longlivedpeople/IFCALongLivedAnalysis.</p>

The starting point are the NTuples. The NTuples that are used to make the plots and the information required is given to the framework in the form of a .dat file.

<p>The plots creation process is divided in two steps: Filling and harvesting:</p>
<ul> 
  <li> In the <strong>filling step</strong> the framework runs over the events of the NTuples, applying the cuts and filling a set of predefined histograms. This step works with instances of the Sample.py, CutManager.py and processHandler.py classes, located in the include/ directory. The filled histograms are stored in several .root files that are given as an input to the harvesting step.</li>
  <li> The <strong>harvesting step</strong> works with instances of both the Sample.py and Canvas.py classes. It accesses the histograms stored in the previously created .root files and make the plots. For this step to succeed the .root files must exist, so the previous step needs to be run at least once time before.</li>
</ul>


### Filling step

To run the <strong>filling step</strong> with the default configuration just run the following command:
```
python fillPlots.py -o [output_name] -d [datacard_name] -m [mode] -c [config_name] (-q)
```
This command will fill the set of predefined histograms and store them in several ```.root``` files (one output ```.root``` file per ```.root``` sample file given as an input).
 
<p> The available options to this command line are </p>
<ul>
  <li> ```-o``` indicates the name of the directory where the ```.root``` files are saved.</li>
  <li> ```-d``` indicates the name of the datacard used to read the ntuples (in dat directory).</li>
  <li> ```-c``` indicates the name of the config file used to apply the selections and regions (explained below).</li>
  <li> ```-m``` indicates the mode in which the histograms are filled among the ones available (explained below).</li>
  <li> ```-q``` is an optional parameter to choose whether the filling step is launhed to CONDOR (preferred) or not.</li>
</ul> 

Before executing the harvesting step, we can <strong>check that all files run correctly in the filling step</strong> just by using the ```checkCompleted.py``` script:
```
python checkCompleted.py -d dir_name
```

#### Dat file

<p> In this framework the samples are managed through the dat files, where all the details needed to create the plots are specified. Usually, each sample will be composed by several .root files (NTuples) that are located in the same directory. The information in the .dat file include: </p>
<ul>
  <li> <strong>Block</strong>: It specifies a <em>block</em> of similar datasets that are drawn jointly e.g. WW, WZ and WZ are drawn together under the block name *Diboson*.</li>
  <li> <strong>Color</strong>: The color of each dataset. </li>
  <li> <strong>Name</strong>: This is the name used to identify each sample. The output histograms used to make the plots will have this identifier. It also serves to select in the running file if one sample is included or not in the plotting.</li>
  <li> <strong>Label</strong>: The text that is put in the legend for each sample (or block). </li>
  <li> <strong>Location of friends</strong>: Absolute path of the directory containing the NTuples (several .root files). </li>
  <li> <strong>Xsec</strong>: Cross-section of the dataset, used to compute the luminosity weights. Should be filled in fb^{-1}.</li>
  <li> <strong>IsData</strong>: 1 if the sample is a data sample or 0 if it is Monte Carlo.</li>
</ul>

#### Config file

(to be completed)

#### Running modes

(to be completed)

### Harvesting step

Once finished, the harvesting step can be run by running whatever ```harvesting_*.py``` file available, giving the previous output directory as the input directory e.g.:
```
python harvesting_BackgroundValidation.py -i dir_name 
```
that will create several png/pdf files with the harvested plots.

### Examples

(to be completed)

## Galapago notebooks (deprecated)

To get a deeper undestanding of Galapago and learn how to tune both steps, a set of Python Notebooks are provided in the ```notebooks/``` folder:
<ul>
  <li> Nb0_Quick-start.ipynb: It explains how to fill and harvest one simple histogram from the begining to the end.</li>
  <li> Nb1_Filling-with-processes.ipynb: It explains how to tune the ```include/processHandler.py``` file to modify the filling step.</li>
  <li> Nb2_Cuts-with-CutManager.ipynb: It explains how to apply cuts by using the ```include/CutManager.py``` Galapago module.</li>
</ul>

To run the notebooks, just init session in https://swan002.cern.ch/ (no additional configuration is required) and open a terminal. To load the framework just clone the repository from the last branch updated:
```
git clone https://github.com/longlivedpeople/Galapago-Framework.git -b winterImplementations
```
and every notebook should run correctly.

## Macros

Several macros that make use of Galapago processing sequences are stored in folder macros/. They are scripts that perform several measurement e.g. specific efficiency plots.

### Trigger-studies

Scripts developed to measure the analysis trigger efficiencies in both data and simulation and extract scale factors. The datasets are included in MET dataset and their events were collected with MET triggers. The simulation is the same dileptonic ttbar simulation used in the analysis.

How to run:
```
python MuonTrigger-efficiency.py -e [2016, 2016APV, 2018]
python PhotonTrigger-efficiency.py -e [2016, 2016APV, 2017, 2018]
```


## Galapago details

### More about reading the Ntuples: The .dat file structure

This framework structures the samples by using 3 python classes: ```Sample```, ```Block``` and ```Tree```. They are defined in ```include/Sample.py```. They are used both in filling and harvesting steps.

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


