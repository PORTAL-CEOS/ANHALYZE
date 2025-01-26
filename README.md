# ANHALYZE

`ANHALYZE` is an analysis tool for the ANHA configuration of the NEMO model.

## Description

`ANHALYZE` is an analysis tool for the 
[Arctic and Northern Hemisphere Atlantic (ANHA)](https://canadian-nemo-ocean-modelling-forum-commuity-of-practice.readthedocs.io/en/latest/Institutions/UofA/Configurations/ANHA4/index.html) 
configuration of the [NEMO](https://www.nemo-ocean.eu/) model. The focus is on data manipulation, analysis, and visualization. 

This tool has been originally developed to support ocean modelling research at the 
Centre for Earth Observation Science (CEOS), at the University of Manitoba. 
Current efforts are ongoing to support the wider range of users of ANHA data.  

If you want to request a new feature, or let us know of a bug, please do so in the [issues](https://github.com/PORTAL-CEOS/ANHALYZE/issues).

The latest progress can be found in [`CHANGELOG`](https://github.com/PORTAL-CEOS/ANHALYZE/blob/main/CHANGELOG.md).

-----
## Getting Started

### Requirements

* Python (> 3.10.0)


### Installation

0. New to git/python world? Some notes for beginners.

    **About using git:**

    - For some tutorials on how to use Git, see 
    this [w3 tutorial](https://www.w3schools.com/git/default.asp?remote=github),    
    this [datacamp tutorial](https://www.datacamp.com/blog/how-to-learn-git),
    this [NHS Git guide](https://nhsdigital.github.io/rap-community-of-practice/training_resources/git/using-git-collaboratively/),
    or this [git tutorial](https://git-scm.com/docs/gittutorial).

    **About setting up your virtual environment before installing this package:**

    - Set up your environment using [pip](https://pypi.org/project/pip/).
    
    - For more information on how to use virtual environments and why they are important, 
    see this [Real Python tutorial](https://realpython.com/python-virtual-environments-a-primer/), or 
    see the [NHS virtual environments guide](https://nhsdigital.github.io/rap-community-of-practice/training_resources/python/virtual-environments/why-use-virtual-environments/).

    - If you are using Linux/MacOS you can use this simple example, 
    or for more information about setting your venv look [here](https://nhsdigital.github.io/rap-community-of-practice/training_resources/python/virtual-environments/venv/).

       ```
       python3 -m venv <venv-directory>
       source <venv-directory>/bin/activate
       python3 -m pip install -r requirements.txt 
       ```
   
    - If you are managing multiple venvs, you could use 
    [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/) as an organizing tool.



1. Clone this [repository](https://github.com/PORTAL-CEOS/ANHALYZE) with: 

    ```
    git clone git@github.com:PORTAL-CEOS/ANHALYZE.git
    ```

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

    - The minimum libraries required for this project are listed in `requirements.in` which was created
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

    - If after following the instructions above you are unable to import the anhalyze module 
      from your jupyter notebook. It maybe that you need to add the kernel of your python virtual environment. 
      You can follow [this](https://medium.com/@WamiqRaza/how-to-create-virtual-environment-jupyter-kernel-python-6836b50f4bf4), 
      [this](https://janakiev.com/blog/jupyter-virtual-envs/), 
      [this](https://www.hophr.com/tutorial-page/getting-import-error-jupyter-notebook-but-not-python-step-by-step-guide),
      or [this](https://cloudbytes.dev/snippets/run-jupyter-notebooks-with-python-virtual-environments) instructions.
      This may happen particularly if you are using multiple virtual environments. 
      
## Usage

Once installed you can import the library like this:

```
import anhalyze as ah
```

### Example 

For first users, we recomend you have a look at the [`anhalyze_tutorial`](https://github.com/PORTAL-CEOS/ANHALYZE/blob/plotting_dev/anhalyze/tutorials/anhalyze_tutorial.ipynb).

This is a short example and requires an `ANHA?-??????_y????m??d??_gridT.nc` file:

```
import anhalyze as ah

filename = 'ANHA?-??????_y????m??d??_gridT.nc'

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
and edit paths to your needs:
``` 
#------------------------------------------------------------- 
#ANHALYZE setup
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

`ANHALYZE` will first try to find your own mask. If none is provided, 
it will attempt to download the default mask and use that.
If you don't want this behaviour then edit the configuration file  `package_data.toml`,
located in `<package root path>/anhalyze/config/`. In the `[mask]` section, change
variable `autodownload_file` from 'true' to 'false'. This will prevent the file from being downloaded.
Then, if you don't provide a valid mask alternative, the code will return an error message.
This assumes the default mask has not been downloaded already. 
If that is the case, you will need to delete it manually. 

-----


## Future Work
    
* To release `AnhalyzeProject` class.
* To release `AnhalyzeLocation` class.
* PyPi packaging, for `pip install` functionality


## License

GNU-GPL license, see [`LICENSE`](https://github.com/PORTAL-CEOS/ANHALYZE/blob/main/LICENSE).