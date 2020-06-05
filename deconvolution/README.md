# Deconvolution  
This directory contains code written for the cell type deconvolution of bulk RNA-seq data. This project is part of my graduation internship for the master Data Science for Life Sciences at the Hanze University of Applied Sciences. For more information, please see my thesis.


## Introduction
The goal of this project is to find interaction effects between eQTLS and covariates. Each step of the pipeline is called from the base directory. Each step also contains a command line argument to adjust to settings file. With this option, any setting can be altered by simply copying the default settings and replacing the variables as you please. Make sure to call the program with the right settings using the **-s** / **--settings** command line option.

Some test scripts and job files are also part of this repository. These files consist code not worthy of a complete step in the pipeline but are still relevant to the project.

## Prerequisites  

This program is developed in Pycharm 2019.3 (Professional Edition), performance on other system is not guaranteed.

The program requires the following packages to be installed:  

 * pandas ([v1.0.1](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))  
 * numpy ([v1.13.3](https://pypi.org/project/numpy/#history); [BSD License](https://www.numpy.org/license.html))  
 * scipy ([v1.4.1](https://docs.scipy.org/doc/scipy/reference/release.html); [BSD License](https://www.scipy.org/scipylib/license.html))  
 * matplotlib ([v3.1.3](https://github.com/matplotlib/matplotlib/releases); [BSD License](https://matplotlib.org/3.1.3/users/license.html))  
 * seaborn ([v0.10.0](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE))  
 * scikit-learn ([v0.22.2.post1](https://scikit-learn.org/stable/whats_new.html); [BSD 3-Clause License](https://github.com/scikit-learn/scikit-learn/blob/master/COPYING))
 * colour ([v0.1.5](https://pypi.org/project/colour/#history); [BSD 3-Clause License](https://pypi.org/project/colour/))  
 * xlrd* ([v1.2.0](https://pypi.org/project/xlrd/#history); [BSD License](https://pypi.org/project/xlrd/)) 
 * venn [v0.1.3](https://pypi.org/project/venn/#history); [GPLv3](https://pypi.org/project/venn/))
  
See 'Installing' on how to install these packages.

\* xlrd is only used on one of the test-scripts and is not required for the main program. 

## Installing  

Installing and running the program can be done by executing the following steps:

**Step 1: acquire the source files**      
Either [clone](https://bitbucket.org/martijnvochteloo/programming1/commits/all) or [download](https://bitbucket.org/martijnvochteloo/programming1/downloads/) the source files.

**Step 2: creating the virtual environment [OPTIONAL]**    
1) Open JetBrains PyCharm  
2) Click on 'File' -> 'Settings' (Ctrl + Alt + S) -> 'Project: progr1' -> 'Project Interpreter'  
3) Add a new interperter by pressing the gear icon -> 'Add'  
4) Make sure 'New environment' is selected, click on 'OK' and again click on 'OK'  
5) Activate the virtual environment by either restarting PyCharm or executing the following command in the terminal:  
```console
source venv/bin/activate
```  

**Step 3: installing the required packages**   
  
1) Make sure the virtual environment is active and run the command  
```console  
pip install -r requirements.txt
```  
2) Wait for the command to finish and Pycharm has updated the indices  


## Usage  
  
### Step 1: Matrix Preparation  

This step creates ordered matrices based on an eQTL file. For each eQTL, the corresponding genotype dosage values and gene expression values are matched. Furthermore, sample identifiers are made uniform.
  
This step also performs marker gene extractions, marker gene factorization, and partial deconvolution; all being part of the covariate matrix creation. 
   
 Other decommissioned steps involve masking the sample identifiers of the matrices, as well as creating groupings of non-NA containing eQTLs, and a correlation matrix for z-score comparison. 

Settings: [default_settings.json](matrix_preparation/settings/default_settings.json)  
Syntax:
```console  
python3 ./matrix_preparation.py
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-f** / **--force_steps**: The steps to force the program to redo, default: None.
 
  
### Step 2: Analyse Interactions 

This step performs multiple regression (MLR) analysis on the ordered expression/genotype matrix for each covariate in the covariate matrix (data from step 1). Each covariate is tested for having an interaction effect with the eQTL. Two distinct methods are implemented: (1) a wrapper around the java bases [eQTLInteractionAnalyser](https://github.com/molgenis/systemsgenetics/wiki/Discovery-of-hidden-confounders-of-QTLs) implementation developed by Patrick Deelen, and (2) my own custom implementation.  

My own implementation also performs permutation based FDR in the MLR analysis. Results are separated on technical covariates and covariates of interest. The former of which is used to validate the model.

#### Method A: [eQTLInteractionAnalyser](https://github.com/molgenis/systemsgenetics/wiki/Discovery-of-hidden-confounders-of-QTLs)
##### Step A1: Analyse Interaction per Group
Settings: [default_settings.json](analyse_interactions/settings/default_settings.json)  
Syntax:
```console  
python3 ./analyse_interactions.py
```  
Options:

 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-g** / **--groups**: The name of the group to analyse, default: 'all'.
 * **-force**: Force the program to redo all steps, default: False.
 * **-verbose**: Include steps and command prints, default: False. 
 
#### Step A2: Merge Groups  
Settings: [default_settings.json](merge_groups/settings/default_settings.json)  
Syntax:
```console  
python3 ./merge_groups.py
```  
Options:

 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-g** / **--groups**: The name of the group to analyse, default: 'all'.
 * **-force**: Force the program to redo all steps, default: False.
    

### Method B: Custom Interaction Analyser
Settings: [default_settings.json](custom_interaction_analyser/settings/default_settings.json)  

##### Step B1: Analyse All Interaction 
This step performs the interaction analyses on a partition of the complete data frame and saves the result as pickled files.
  
Syntax:
```console  
python3 ./custom_interaction_analyser.py
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-sr** / **--skip_rows**: The number of rows to skip in the input files, default: 0
 * **-ne** / **--n_eqtls**: The number of eQTLs in the input files, default: None (determine manually).
 * **-ns** / **--n_samples**: The number of samples in the input files, default: None (determine manually).
 * **-verbose**: Include steps and command prints, default: False.  

Example:
 * Job1, analyzes eQTLs 0-50: 
 ```console  
python3 ./custom_interaction_analyser.py -ne 50 
```  
 * Job2, analyzed eQTLs 50-100: 
 ```console  
python3 ./custom_interaction_analyser.py -sr 50 -ne 100 
``` 
 * Job3:   
   ...  
   
##### Step B2: Combine the Resuts
This step loads the pickled data and combines them into a complete interaction matrix. Also multiple-testing corrections are performed and the resulting FDR are compared to the original p-values. This code also creates a few visualizations of the interaction data.
  
Syntax:
```console  
python3 ./custom_interaction_analyser.py -combine
```  
Options:
 * **-combine**: Combine the created files, alternative functionality. Default: False.    
  
### Step 3: visualiser  

This step takes the results from the previous two steps and creates visualizations of it. Each plot can be called separately.

Settings: [default_settings.json](visualiser/settings/default_settings.json)  
Syntax:
```console  
python3 ./visualiser.py
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-a** / **--alpha**: The significance cut-off, default: 0.05.
 * **-p** / **--plots**: The name of the figures to be created, default: 'all'.
 * **-t** / **--top**: The number of top eQTLs to visualise, default: 1.
 * **-validate**: Validate that the input matrices match with each other and then quit, default: 'False'.  
  
## Author  

Martijn Vochteloo *(1)*

1. Institute for Life Science & Technology, Hanze University of Applied Science, Groningen, The Netherlands.

## License  

This project is licensed under the GNU GPL v3 license - see the [LICENSE.md](LICENSE.md) file for details
