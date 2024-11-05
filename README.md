# ANHALIZE  (beta)

`ANHALIZE` is an analysis tool for the ANHA configuration of the NEMO model.

## Description

`ANHALIZE` is an analysis tool for the 
[Arctic and Northern Hemisphere Atlantic (ANHA)](https://canadian-nemo-ocean-modelling-forum-commuity-of-practice.readthedocs.io/en/latest/Institutions/UofA/Configurations/ANHA4/index.html) 
configuration of the [NEMO](https://www.nemo-ocean.eu/) model. The focus is on data manipulation, analysis, and visualization. 

This tool has been originally developed to support ocean modeling research at the 
Centre for Earth Observation Science (CEOS), at the University of Manitoba. 

NOTE: This code is stable, but new features are currently under development.



-----
## Getting Started

### Requirements

* Python (> 3.10.0)


### Installation

1. Clone this [repository](https://github.com/PORTAL-CEOS/ANHALIZE) with: 

    ```
    git clone git@github.com:PORTAL-CEOS/ANHALYZE.git
    ```

    **Additional Notes**

    To learn about what this means, and how to use Git, see 
    this [w3 tutorial](https://www.w3schools.com/git/default.asp?remote=github),    
    this [datacamp tutorial](https://www.datacamp.com/blog/how-to-learn-git),
    this [NHS Git guide](https://nhsdigital.github.io/rap-community-of-practice/training_resources/git/using-git-collaboratively/),
    or this [git tutorial](https://git-scm.com/docs/gittutorial).


2. Install the dependencies in your [python environment](https://docs.python.org/3/library/venv.html) with:
    ```
    cd ANHALYZE/
    python3 -m pip install -r requirements.txt
    ```
   
    Then install the package with:
    ```
    python3 -m pip install .
    ```
    
    If you want to install the package in development mode, you can do:
    ```
    python3 -m pip install  --editable .
    ```

    **Additional Notes**

    - Note1: If you are managing multiple venvs, you could use 
    [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/) as an organizing tool.
    - Note2: The minimum libraries required for this project are listed in `requirements.in` which was created
    automatically with
       ```
       pip install pipreqs
       pipreqs <project-directory> --savepath requirements.in --scan-notebooks    
       ```   
       Additionally, the `jupyter` library was added manually. 
   
       The required libraries including their dependencies are listed in `requirements.txt`.
       This was created automatically from `requirements.in` by using
       ```
       pip install pip-tools
       pip-compile    
       ```

    **Additional Help**

    Set up your environment using [pip](https://pypi.org/project/pip/).
    For more information on how to use virtual environments and why they are important, 
    see this [Real Python tutorial](https://realpython.com/python-virtual-environments-a-primer/), or 
    see the [NHS virtual environments guide](https://nhsdigital.github.io/rap-community-of-practice/training_resources/python/virtual-environments/why-use-virtual-environments/).

    If you are using Linux/MacOS you can use this simple example, 
    or for more information about setting your venv look [here](https://nhsdigital.github.io/rap-community-of-practice/training_resources/python/virtual-environments/venv/).

       ```
       python3 -m venv <venv-directory>
       source <venv-directory>/bin/activate
       python3 -m pip install -r requirements.txt 
       ```


## Usage

Once installed you can import the library like this:

```
import anhalyze as ah
```

### Example 

This example requires an `ANHA*_gridT.nc` file.

```
import anhalyze as ah

filename = 'ANHA*gridT.nc'

# Open file
aa = ah.AnhaDataset(filename)

# Do selection of lat_range, lon_range, and depth_range.
aaa = aa.sel(lat_range=[50,65],lon_range=[-93,-75],depth_range=[0,300])

# Plot region for a selected variable.
aaa.show_var_data_map(var='votemper')
``` 

## Additional Advanced Notes

### Masking

Masking is done internally. The default mask file will be automatically 
downloaded here `<package root path>/anhalyze/package_data/`,
the default mask name is `ANHA4_mask.nc`.

There are two alternate options for providing your own mask file:
 
#### Environment variable option:

You can add the following environmental variable(s) to your `.bash_profile` (or `.bashrc`, etc..), 
and edith paths to your needs:
``` 
#------------------------------------------------------------- 
#ANHALIZE setup
#-------------------------------------------------------------
export MASK_FILENAME='/root_path/user/../your_mask_location/mask_name.nc'
#-------------------------------------------------------------
```

#### Dynamic variable option:

You can provide your mask filename with full path when opening
an `ANHA*.nc` file, as follows:

```
import anhalyze as ah

# Open file
aa = ah.AnhaDataset(filename, mask_filename='mask_full_filename')
```

#### No-autodownload option:

`Anhalize` will first try to find your own mask, if none is provided, 
it will attempt to download the default mask and use that.
If you don't want this behaviour then edit the configuration file  `package_data.toml`,
located in `<package root path>/anhalyze/config/`. In the `[mask]` section, change
variable `autodownload_file` from 'true' to 'false'. This will prevent the file to be downloaded.
If you don't provide a valid mask alternative, the code will return an error message.
This assumes the default mask has not been downloaded already. 
If that is the case, you will need to delete it manually. 

-----
## Version History

* 0.7 (planned)
    * To release `AnhalyzeLocation` class.
    * pip install functionality
* 0.6 (planned)
    * To release `AnhalyzeProject` class.
* 0.5 (upcoming)
    * First main release, includes:
      * `AnhaDataset` class
      * Basic functionality of single ANHA `nc` files. 
        * data reading/writing
        * region selection 
        * map plotting
    * See [commit change]() or See [release history]()
* 0.0.1 (current)
    * Beta version, initial development.  

## License

This is a placeholder.



-------------------- Deprecated below this line. --------------------


Need to add the following environmental variables to your .bash_profile (or .bashrc, etc..), 
and edith paths to your needs:
``` 
#------------------------------------------------------------- 
#ANHALIZE setup
#-------------------------------------------------------------
export DATA_PATH='/root_path/user/NEMO/ANHA4/ANHA4-Wxx000-S/'
#-------------------------------------------------------------
```

