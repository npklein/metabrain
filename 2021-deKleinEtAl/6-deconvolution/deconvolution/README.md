# Brain eQTL Deconvolution  
This directory contains code written for the cell type deconvolution of bulk RNA-seq data. This project is part of my graduation internship for the master Data Science for Life Sciences at the Hanze University of Applied Sciences. For more information, please see my thesis.


## Introduction
The goal of this project is to find interaction effects between eQTLS and covariates. Each step of the pipeline is called from the base directory. Each step also contains a command line argument to adjust the input / output directory name and the settings file. The **-n** / **--name** command line option defines the name of the input, as well as the output, directory. This name will be appended on the settings ***input_dir** variable(s). The reason for defining this setting as a command line argument is so subsets of data can be analyzed using the same settings file. With the settings file, any global variable can be altered in the program. Simply copy the default settings and replacing the variables as you please. Make sure to call the program with the right settings using the **-s** / **--settings** command line option.

## Prerequisites  

This program is developed in Pycharm 2019.3 (Professional Edition), performance on other system is not guaranteed.

The program requires the following packages to be installed:  

 * pandas ([v1.0.1](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))  
 * numpy ([v1.18.3](https://pypi.org/project/numpy/#history); [BSD License](https://www.numpy.org/license.html))  
 * scipy ([v1.4.1](https://docs.scipy.org/doc/scipy/reference/release.html); [BSD License](https://www.scipy.org/scipylib/license.html))
 * statsmodels ([v0.11.1](https://www.statsmodels.org/stable/index.html); [Modified BSD (3-clause) license](https://www.statsmodels.org/stable/index.html))    
 * matplotlib ([v3.2.1](https://github.com/matplotlib/matplotlib/releases); [BSD License](https://matplotlib.org/3.1.3/users/license.html))  
 * seaborn ([v0.10.0](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE))  
 * scikit-learn ([v0.22.2.post1](https://scikit-learn.org/stable/whats_new.html); [BSD 3-Clause License](https://github.com/scikit-learn/scikit-learn/blob/master/COPYING))
 * colour ([v0.1.5](https://pypi.org/project/colour/#history); [BSD 3-Clause License](https://pypi.org/project/colour/))
 * upsetplot [v0.4.0](https://pypi.org/project/UpSetPlot/#history); [BSD 3-Clause License](https://pypi.org/project/UpSetPlot//))   
 * xlrd* ([v1.2.0](https://pypi.org/project/xlrd/#history); [BSD License](https://pypi.org/project/xlrd/)) 

See 'Installing' on how to install these packages.

\* xlrd is only used on one of the test-scripts and is not required for the main program. 

## Installing  

Installing and running the program can be done by executing the following steps:

**Step 1: acquire the source files**      
Either [clone](https://github.com/mvochteloo/brain_eQTL.git) or [download](https://github.com/mvochteloo/brain_eQTL/archive/master.zip) the source files.

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
   
 Other decommissioned steps involve masking the sample identifiers of the matrices, as well as creating groupings of non-NA containing eQTLs, and a correlation matrix for z-score comparison. This is no longer used.

Settings: [default_settings.json](matrix_preparation/settings/default_settings.json)  
Syntax:
```console  
python3 ./matrix_preparation.py -n example_output
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-d** / **--disease**: The name of the disease to filter on, default: '' (i.e. no filter).
 * **-f** / **--force_steps**: The steps to force the program to redo, default: None.
 
  
### Step 2: Analyse Interactions 

This step performs multiple regression (MLR) analysis on the ordered expression/genotype matrix for each covariate in the covariate matrix (data from step 1). Each covariate is tested for having an interaction effect with the eQTL. Two distinct methods are implemented: **(1)** a wrapper around the java based [eQTLInteractionAnalyser](https://github.com/molgenis/systemsgenetics/wiki/Discovery-of-hidden-confounders-of-QTLs) implementation developed by Patrick Deelen, and **(2)** my own custom implementation. The first of which is decommissioned but is still included in the repository as it is part of process and therefore relevant for grading. 
  
My own implementation is what is actually used. This implementation of MLR also performs permutation based FDR in the MLR analysis. Results are separated on technical covariates and covariates of interest. The former of which are used to validate that the model works as intended.

Settings: [default_settings.json](custom_interaction_analyser/settings/default_settings.json)  

##### Step 2A: Multiple Linear Regression Analysis 
This step performs the interaction analyses on a partition of the complete data frame and saves the result as pickled files.
  
Syntax:
```console  
python3 ./custom_interaction_analyser.py -n example_output
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-sr** / **--skip_rows**: The number of rows to skip in the input files, default: 0
 * **-ne** / **--n_eqtls**: The number of eQTLs in the input files, default: None (determine automatically).
 * **-ns** / **--n_samples**: The number of samples in the input files, default: None (determine automatically).
 * **-verbose**: Include steps and command prints, default: False.  

Example:
 * Job1, analyzes eQTLs 0-50: 
 ```console  
python3 ./custom_interaction_analyser.py -n example_output -ne 50 
```  
 * Job2, analyzed eQTLs 50-100: 
 ```console  
python3 ./custom_interaction_analyser.py -n example_output -sr 50 -ne 100 
``` 
 * Job3:   
   ...  
   
**TIP**: to automate this process I developed three bash scripts that can help with this step. First, using the [create_CIA_jobs](jobs/create_CIA_jobs.py) script you can create Custom Interaction Analyser job files to submit to the cluster. 

Options:

 * **-j** / **--job**: The name of the job.
 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-f** / **--first**: The eQTL index of the first one to analyse, default: 0.
 * **-l** / **--last**: The eQTL index of the last one to analyse.
 * **-b** / **--batch**: The number of eQTLs per job. Default: 50.
 * **-ns** / **--n_samples**: The number of samples in the expression file.
 * **-e** / **--exclude**: The name of node to exclude, default: None.  
 
**Important note**: the maximum runtime of these job files is '05:59:00'. If the number of permutations in performed or the **-b** / **--batch** gets too big, the process won't be finished in time. I recommend using <75 for 10 permutations.

This program then creates N files named <job><n>.sh. After creating the files, you can start submitting them. IMPORTANT that you first submit <job>_0.sh and wait for it to start. You know it has started with 'custom_interaction_analyser/<output_directory>/permutation_order.pkl' exists. This file has to exist before continuing! 
You can then submit the rest of the job files using [start](jobs/start.sh):
 ```console  
./start.sh <job_prefix> <start_index> <stop_index>  
```  
Where job_prefix = **-j** / **--job**, start_index = 1, and stop_index = the highest job file number.
  
During the analysis you can check on the progress by using [status](jobs/status.sh):
 ```console  
./status.sh <job_prefix> <start_index> <stop_index>  
```  

If something went wrong, you can stop jobs easily using [stop](jobs/stop.sh):
 ```console  
./stop.sh
```  
This program stops all jobs that are written to the start.txt file.
      
##### Step 2B: Combine the Resuts
This step loads the pickled data and combines them into a complete interaction matrix. Also multiple-testing corrections are performed and the resulting FDR values are compared to the original p-values. This code also creates a few visualizations of the p-value distributions and fdr - pvalue comparisons.
  
Syntax:
```console  
python3 ./custom_interaction_analyser.py -n example_output -combine
```  
Options:
 * **-combine**: Combine the created files, alternative functionality. Default: False.    
  
### Step 3: Identify Cell Type Mediated eQTLs  

This step takes the results from the previous two steps and groups eQTLs based on cell type mediated effects. Furthermore, upsetplot(s) is/are created for overlapping interaction effects.

Settings: [default_settings.json](identify_ct_mediated_eqtls/settings/default_settings.json)  
Syntax:
```console  
python3 ./identify_ct_mediated_eqtls.py -n example_output
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-a** / **--alpha**: The significance cut-off, default: 0.05.
 * **-e** / **--extensions**: The output file formats, default: ['png'].
 * **-i** / **--interest**: The HGNC names to print the info of, default: None.  

### Step 4: Visualiser  

This step takes the results from the first two steps and creates visualizations of it. Each plot can be called separately.

Settings: [default_settings.json](visualiser/settings/default_settings.json)  
Syntax:
```console  
python3 ./visualiser.py -n example_output
```  
Options:

 * **-n** / **--name**: The name of the input/output directory.
 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-a** / **--alpha**: The significance cut-off, default: 0.05.
 * **-p** / **--plots**: The name of the figures to be created, default: 'all'.
 * **-t** / **--top**: The number of top eQTLs to visualise, default: 1.
 * **-i** / **--interest**: The indices of the eQTLS to visualise, default: None. If set, -t / --top is discarded.
 * **-e** / **--extension**: The output file format, default: 'png'.
 * **-validate**: Validate that the input matrices match with each other and then quit, default: 'False'.  
  
## Questions and Answers
**Q**: You are refering to your thesis; can I view it?  
**A**: The thesis of this project is under embargo until the 1st of July 2021.  
  
**Q**: Why are there more files in the repository than described in this README?  
**A**: This README only describes the steps for performing cell type deconvolution. However, some code contributes to the project but not directly to the cell type deconvolution pipeline. Some scripts are written to test a certain aspect of the results [test_scripts](test_scripts), others are made to create and manage slurm jobs [jobs](jobs), and other code such as [analyse_interactions](analyse_interactions) en [merge_groups](merge_groups) is no longer part of the pipeline. All this code is not strictly required for the project, but is included for the grading of the repository. One exception is the [general](general) folder. This code contains classes that are used for more than one step is  part of the pipeline although not explicitly mentioned in the steps.

**Q**: How to cite? (NOTE: only for citing the code in this repository; cite the results of the project as it is published)  
**A**: *'Vochteloo, M. (2020). Brain eQTL Deconvolution. Hanze University of Applied Sciences, Groningen, The Netherlands.'*  

## Author  

Martijn Vochteloo *(1)*

1. Institute for Life Science & Technology, Hanze University of Applied Science, Groningen, The Netherlands.

## License  

This project is licensed under the GNU GPL v3 license - see the [LICENSE.md](LICENSE.md) file for details
