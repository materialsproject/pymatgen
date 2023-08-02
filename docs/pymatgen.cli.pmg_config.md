---
layout: default
title: pymatgen.cli.pmg_config.md
nav_exclude: true
---

# pymatgen.cli.pmg_config module

Implementation for pmg config CLI.


### pymatgen.cli.pmg_config.add_config_var(tokens: list[str], backup_suffix: str)
Add/update keys in .pmgrc.yaml config file.


### pymatgen.cli.pmg_config.build_bader(fortran_command='gfortran')
Build bader package.


* **Parameters**

    **fortran_command** –



### pymatgen.cli.pmg_config.build_enum(fortran_command: str = 'gfortran')
Build enum.


* **Parameters**

    **fortran_command** –



### pymatgen.cli.pmg_config.configure_pmg(args: Namespace)
Handle configure command.


### pymatgen.cli.pmg_config.install_software(install: Literal['enumlib', 'bader'])
Install all optional external software.


### pymatgen.cli.pmg_config.setup_cp2k_data(cp2k_data_dirs: list[str])
Setup CP2K basis and potential data directory.


### pymatgen.cli.pmg_config.setup_potcars(potcar_dirs: list[str])
Setup POTCAR directories.