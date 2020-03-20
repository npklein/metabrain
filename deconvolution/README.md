# Deconvolution  

## Introduction


## Prerequisites  

This program is developed in Pycharm 2019.3 (Professional Edition), performance on other system is not guaranteed.

The program requires the following packages to be installed:  

 * pandas ([v0.24.2](https://github.com/pandas-dev/pandas); [BSD 3-Clause License](https://github.com/pandas-dev/pandas/blob/master/LICENSE))  
 * numpy ([v1.18.1](https://github.com/numpy/numpy/releases); [BSD License](https://www.numpy.org/license.html))  
 * scipy ([v1.4.1](https://docs.scipy.org/doc/scipy/reference/release.html); [BSD License](https://www.scipy.org/scipylib/license.html))  
 * matplotlib ([v3.0.3](https://github.com/matplotlib/matplotlib/releases); [BSD License](https://matplotlib.org/3.1.0/users/license.html))  
 * seaborn ([v0.9.1](https://github.com/mwaskom/seaborn); [BSD 3-Clause License](https://github.com/mwaskom/seaborn/blob/master/LICENSE))  
 

See 'Installing' on how to install these packages.

## Installing  

Installing and running the program can be done by executing the following steps:

**Step 1: acquire the source files**      
Either [clone](https://bitbucket.org/martijnvochteloo/programming1/commits/all) or [download](https://bitbucket.org/martijnvochteloo/programming1/downloads/) the source files.

**Step 2: creating the virtual environment [OPTINAL]**    
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
  
### Step 1: matrix preparation  
Settings: [default_settings.json]('matrix_preparation/settings/default_settings.json')  
Syntax:
```console  
python3 ./matrix_preparation.py
```  
Options:

 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-f** / **--force_steps**: The steps to force the program to redo, default: None.
 
  
### Step 2: analyse interactions  
Settings: [default_settings.json]('analyse_interactions/settings/default_settings.json')  
Syntax:
```console  
python3 ./analyse_interactions.py
```  
Options:

 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-g** / **--groups**: The name of the group to analyse, default: 'all'.
 * **-force**: Force the program to redo all steps, default: False.
 * **-verbose**: Include steps and command prints, default: False. 
 
### Step 3: merge groups  
Settings: [default_settings.json]('merge_groups/settings/default_settings.json')  
Syntax:
```console  
python3 ./merge_groups.py
```  
Options:

 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-g** / **--groups**: The name of the group to analyse, default: 'all'.
 * **-force**: Force the program to redo all steps, default: False.
    
### Step 4: visualiser  
Settings: [default_settings.json]('visualiser/settings/default_settings.json')  
Syntax:
```console  
python3 ./visualiser.py
```  
Options:

 * **-s** / **--settings**: The settings input file (without '.json'), default: 'default_settings'.
 * **-p** / **--plots**: The name of the figures to be created, default: 'all'.
 * **-n** / **--n_eqtls**: The number of eQTLs to visualise, default: 1.
 * **-validate**: Validate that the input matrices match with each other and then quit, default: 'False'.  
  
## Author  

Martijn Vochteloo *(1)*

1. Institute for Life Science & Technology, Hanze University of Applied Science, Groningen, The Netherlands.

## License  

This project is licensed under the GNU GPL v3 license - see the [LICENSE.md](LICENSE.md) file for details
