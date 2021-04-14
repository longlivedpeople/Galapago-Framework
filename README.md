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

## How this framework works: filling and harvesting steps

The starting point are the NTuples. The NTuples that are used for making the plots and the information required is given to the framework in the form of a .dat file.

<p>The plots creation process is divided in two steps: Filling and harvesting:</p>
<ul> 
  <li> In the <strong>filling step</strong> the framework runs over the events of the NTuples, applying the cuts and filling a set of predefined histograms. This step works with instances of the Sample.py, CutManager.py and processHandler.py classes, located in the include/ directory. The filled histograms are stored in several .root files that are given as an input to the harvesting step.</li>
  <li> The <strong>harvesting step</strong> works with instances of both the Sample.py and Canvas.py classes. It accesses the histograms stored in the previously created .root files and make the plots. For this step to succeed the .root files must exist, so the previous step needs to be run at least once time before.</li>
</ul>


## Quick start

To run the <strong>filling step</strong> with the default configuration just run the following command:
```
python fillPlots.py -o dir_name 
```
this will fill the set of predefined histograms declared in ```include/processHandler.py``` and store them in several ```.root``` files (one output ```.root``` file per ```.root``` sample file given as an input). By default this command runs on CONDOR, so one job will be launched per ```.root``` file. The output files will be saved in ```dir_name/``` that will be automatically created.

Once finished, the harvesting step can be run by running whatever ```harvesting_*.py``` file available, giving the previous output directory as the input directory e.g.:
```
python harvesting_BackgroundValidation.py -i dir_name 
```
that will create several png/pdf files with the harvested plots.

Before executing the harvesting step, we can <strong>check that all files run correctly in the filling step</strong> just by using the ```checkCompleted.py``` script:
```
python checkCompleted.py -d dir_name
```

### Galapago notebooks

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

## Reading the Ntuples: The .dat file structure

<p> Galapago works with Ntuples created with the IFCALongLived analyzer, that can be accesed here: https://github.com/longlivedpeople/IFCALongLivedAnalysis.</p>

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

### Using the .dat file to read samples

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


