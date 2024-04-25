# ANHALIZE

This software supports ocean modeling research at the Centre for Earth Observation Science (CEOS). 

The focus is on data manipulation, analysis and visualization of the 
[Arctic and Northern Hemisphere Atlantic (ANHA)]((https://canadian-nemo-ocean-modelling-forum-commuity-of-practice.readthedocs.io/en/latest/Institutions/UofA/Configurations/ANHA4/index.html)) configuration of the NEMO model. 

Currently, the development is concentrated in ANHA data. 





## Set up:

Install with

`python -m pip install .`





-------------------- Deprecated below this line. --------------------


Need to add the following environmental variables to your .bash_profile (or .bashrc, etc..), 
and edith paths to your needs:
``` 
#------------------------------------------------------------- 
#ANHALIZE setup
#-------------------------------------------------------------
export MASK_PATH='/root_path/user/ANALYSES/MASKS/'
export DATA_PATH='/root_path/user/NEMO/ANHA4/ANHA4-Wxx000-S/'
#-------------------------------------------------------------
```

