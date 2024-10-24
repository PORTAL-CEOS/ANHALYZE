#
# Setting up config file based on
# https://realpython.com/python-toml/#use-configuration-files-in-your-projects
#

import pathlib
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

path = pathlib.Path(__file__).parent / "package_data.toml"
with path.open(mode="rb") as fp:
    package_data = tomllib.load(fp)
