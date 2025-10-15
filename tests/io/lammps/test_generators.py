from __future__ import annotations

import unittest

import pytest

from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData, LammpsForceField
from pymatgen.io.lammps.generators import (
    _BASE_LAMMPS_SETTINGS,
    BaseLammpsSetGenerator,
    LammpsMinimization,
    LammpsSettings,
)
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/io/lammps"


class TestLammpsSettings(unittest.TestCase):
    @classmethod
    def setup_class(cls):
        cls.cif = f"{TEST_DIR}/lgps.cif"
        cls.structure = Structure.from_file(cls.cif)

    def test_initialization_with_dict(self):
        """Test LammpsSettings initialization with dictionary."""
        settings_dict = {
            "units": "metal",
            "atom_style": "full",
            "dimension": 3,
            "boundary": ("p", "p", "p"),
            "ensemble": "nvt",
            "thermostat": "nose-hoover",
            "nsteps": 1000,
        }
        settings = LammpsSettings(**settings_dict)

        assert settings.units == "metal"
        assert settings.atom_style == "full"
        assert settings.dimension == 3
        assert settings.boundary == ("p", "p", "p")
        assert settings.ensemble == "nvt"
        assert settings.thermostat == "nose-hoover"
        assert settings.nsteps == 1000

    def test_initialization_with_validation(self):
        """Test LammpsSettings initialization with parameter validation."""
        # Valid settings should work
        valid_settings = {
            "units": "metal",
            "atom_style": "full",
            "boundary": ("p", "p", "p"),
            "ensemble": "nvt",
            "thermostat": "nose-hoover",
            "barostat": "nose-hoover",
            "min_style": "cg",
            "restart": "restart.file",
        }
        settings = LammpsSettings(validate_params=True, **valid_settings)
        assert settings.units == "metal"

    def test_validation_invalid_units(self):
        """Test validation with invalid units."""
        with pytest.raises(ValueError, match="Error validating key units: set to invalid_unit"):
            LammpsSettings(validate_params=True, units="invalid_unit")

    def test_validation_invalid_atom_style(self):
        """Test validation with invalid atom style."""
        with pytest.raises(ValueError, match="Error validating key atom_style: set to invalid_style"):
            LammpsSettings(validate_params=True, atom_style="invalid_style")

    def test_validation_invalid_boundary(self):
        """Test validation with invalid boundary conditions."""
        with pytest.raises(ValueError, match="Error validating key boundary: set to"):
            LammpsSettings(validate_params=True, boundary=("invalid", "p", "p"))

    def test_validation_invalid_ensemble(self):
        """Test validation with invalid ensemble."""
        with pytest.raises(ValueError, match="Error validating key ensemble: set to invalid_ensemble"):
            LammpsSettings(validate_params=True, ensemble="invalid_ensemble")

    def test_validation_invalid_thermostat(self):
        """Test validation with invalid thermostat."""
        with pytest.raises(ValueError, match="Error validating key thermostat: set to invalid_thermostat"):
            LammpsSettings(validate_params=True, thermostat="invalid_thermostat")

    def test_validation_invalid_barostat(self):
        """Test validation with invalid barostat."""
        with pytest.raises(ValueError, match="Error validating key barostat: set to invalid_barostat"):
            LammpsSettings(validate_params=True, barostat="invalid_barostat")

    def test_validation_invalid_min_style(self):
        """Test validation with invalid minimization style."""
        with pytest.raises(ValueError, match="Error validating key min_style: set to invalid_min_style"):
            LammpsSettings(validate_params=True, min_style="invalid_min_style")

    def test_validation_list_boundary(self):
        """Test validation with list boundary conditions."""
        # Valid list
        settings = LammpsSettings(validate_params=True, boundary=["p", "f", "s"])
        assert settings.boundary == ["p", "f", "s"]

        # Invalid list
        with pytest.raises(ValueError, match="Error validating key boundary: set to"):
            LammpsSettings(validate_params=True, boundary=["invalid", "p", "p"])

    def test_restart_validation(self):
        """Test restart parameter validation."""
        # Valid restart (string)
        settings = LammpsSettings(restart="restart.file")
        assert settings.restart == "restart.file"

        # Invalid restart (not string)
        with pytest.raises(ValueError, match="restart should be the path to the restart file"):
            LammpsSettings(restart=123)

    def test_pressure_validation(self):
        """Test pressure parameter validation."""
        # Valid 3-element pressure list
        settings = LammpsSettings(
            validate_params=True, ensemble="npt", start_pressure=[1.0, 1.0, 1.0], end_pressure=[2.0, 2.0, 2.0]
        )
        assert settings.start_pressure == [1.0, 1.0, 1.0]

        # Invalid pressure list length
        with pytest.raises(ValueError, match="start_pressure should be a list of 3 values"):
            LammpsSettings(
                validate_params=True,
                ensemble="npt",
                start_pressure=[1.0, 1.0],  # Only 2 values
            )

    def test_friction_timestep_warning(self):
        """Test warning when friction is smaller than timestep."""
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            LammpsSettings(
                validate_params=True,
                ensemble="nvt",
                thermostat="nose-hoover",
                friction=0.0001,  # Very small friction
                timestep=0.001,  # Larger timestep
            )
            assert len(w) > 0
            assert "Friction" in str(w[0].message)

    def test_minimization_validation(self):
        """Test minimization-specific validation."""
        # Valid minimization
        settings = LammpsSettings(validate_params=True, ensemble="minimize", nsteps=1000, tol=1e-6)
        assert settings.ensemble == "minimize"

        # Invalid nsteps for minimization
        with pytest.raises(ValueError, match="nsteps should be greater than 0 for minimization simulations"):
            LammpsSettings(validate_params=True, ensemble="minimize", nsteps=0)

    def test_minimization_tolerance_warning(self):
        """Test warning for large minimization tolerance."""
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            LammpsSettings(
                validate_params=True,
                ensemble="minimize",
                nsteps=1000,
                tol=1e-3,  # Large tolerance
            )
            assert len(w) > 0
            assert "Tolerance for minimization" in str(w[0].message)

    def test_as_dict_method(self):
        """Test the as_dict method."""
        settings = LammpsSettings(units="metal", atom_style="full", dimension=3)
        settings_dict = settings.as_dict()

        assert isinstance(settings_dict, dict)
        assert settings_dict["units"] == "metal"
        assert settings_dict["atom_style"] == "full"
        assert settings_dict["dimension"] == 3

    def test_update_method(self):
        """Test the update method."""
        settings = LammpsSettings(units="metal", atom_style="full")
        updates = {"dimension": 3, "ensemble": "nvt"}
        settings.update(updates)

        assert settings.dimension == 3
        assert settings.ensemble == "nvt"
        assert settings.units == "metal"  # Should remain unchanged

    def test_msonable_functionality(self):
        """Test MSONable serialization/deserialization."""
        settings = LammpsSettings(
            units="metal", atom_style="full", dimension=3, boundary=("p", "p", "p"), ensemble="nvt"
        )

        # Test as_dict
        settings_dict = settings.as_dict()
        assert "units" in settings_dict
        assert settings_dict["units"] == "metal"

        # Test from_dict
        new_settings = LammpsSettings.from_dict(settings_dict)
        assert new_settings.units == "metal"
        assert new_settings.atom_style == "full"
        assert new_settings.dimension == 3

    def test_default_values(self):
        """Test that default values are properly set."""
        settings = LammpsSettings(**_BASE_LAMMPS_SETTINGS["periodic"])

        # Check some default values from _COMMON_LAMMPS_SETTINGS
        assert settings.dimension == 3
        assert settings.pair_style == "lj/cut 10.0"
        assert settings.thermo == 100
        assert settings.start_temp == 300.0
        assert settings.end_temp == 300.0
        assert settings.start_pressure == 0.0
        assert settings.end_pressure == 0.0
        assert settings.log_interval == 100
        assert settings.traj_interval == 100
        assert settings.ensemble == "nvt"
        assert settings.thermostat == "nose-hoover"
        assert settings.barostat is None
        assert settings.nsteps == 1000
        assert settings.restart == ""
        assert settings.tol == 1e-6
        assert settings.min_style == "cg"


class TestBaseLammpsSetGenerator(unittest.TestCase):
    @classmethod
    def setup_class(cls):
        cls.cif = f"{TEST_DIR}/lgps.cif"
        cls.structure = Structure.from_file(cls.cif)

    def test_initialization_with_dict_settings(self):
        """Test BaseLammpsSetGenerator initialization with dictionary settings."""
        settings_dict = {
            "units": "metal",
            "atom_style": "full",
            "ensemble": "nvt",
            "nsteps": 1000,
        }
        generator = BaseLammpsSetGenerator(settings=settings_dict)

        assert isinstance(generator.settings, LammpsSettings)
        assert generator.settings.units == "metal"
        assert generator.settings.atom_style == "full"
        assert generator.settings.ensemble == "nvt"
        assert generator.settings.nsteps == 1000

    def test_initialization_with_lammps_settings(self):
        """Test BaseLammpsSetGenerator initialization with LammpsSettings object."""
        settings = LammpsSettings(units="real", atom_style="molecular", ensemble="npt")
        generator = BaseLammpsSetGenerator(settings=settings)

        assert generator.settings is settings
        assert generator.settings.units == "real"
        assert generator.settings.atom_style == "molecular"
        assert generator.settings.ensemble == "npt"

    def test_initialization_with_data_type(self):
        """Test BaseLammpsSetGenerator initialization with different data types."""
        # Periodic data type
        generator_periodic = BaseLammpsSetGenerator(data_type="periodic")
        assert generator_periodic.data_type == "periodic"
        assert generator_periodic.settings.units == "metal"
        assert generator_periodic.settings.atom_style == "atomic"

        # Molecular data type
        generator_molecular = BaseLammpsSetGenerator(data_type="molecular")
        assert generator_molecular.data_type == "molecular"
        assert generator_molecular.settings.units == "real"
        assert generator_molecular.settings.atom_style == "full"

    def test_initialization_with_calc_type(self):
        """Test BaseLammpsSetGenerator initialization with different calc types."""
        generator = BaseLammpsSetGenerator(calc_type="custom_calc")
        assert generator.calc_type == "custom_calc"

    def test_initialization_with_force_field_dict(self):
        """Test BaseLammpsSetGenerator initialization with force field dictionary."""
        force_field_dict = {"pair_style": "lj/cut 10.0", "pair_coeff": ["* * 1.0 1.0"]}
        generator = BaseLammpsSetGenerator(force_field=force_field_dict)

        assert isinstance(generator.force_field, LammpsForceField)
        assert generator.force_field.pair_style == "lj/cut 10.0"

    def test_initialization_with_force_field_object(self):
        """Test BaseLammpsSetGenerator initialization with LammpsForceField object."""
        force_field = LammpsForceField(pair_style="lj/cut 10.0", pair_coeff="* * 1.0 1.0")
        generator = BaseLammpsSetGenerator(force_field=force_field)

        assert generator.force_field is force_field

    def test_template_selection_nvt(self):
        """Test template selection for NVT ensemble."""
        generator = BaseLammpsSetGenerator(settings={"ensemble": "nvt"})
        assert generator.inputfile is not None
        # Should use md.template for NVT

    def test_template_selection_npt(self):
        """Test template selection for NPT ensemble."""
        generator = BaseLammpsSetGenerator(settings={"ensemble": "npt"})
        assert generator.inputfile is not None
        # Should use md.template for NPT

    def test_template_selection_minimize(self):
        """Test template selection for minimization ensemble."""
        generator = BaseLammpsSetGenerator(settings={"ensemble": "minimize"})
        assert generator.inputfile is not None
        # Should use minimization.template for minimize

    def test_template_selection_invalid_ensemble(self):
        """Test error for invalid ensemble."""
        with pytest.raises(ValueError, match="Unknown ensemble='invalid'"):
            BaseLammpsSetGenerator(settings={"ensemble": "invalid"})

    def test_include_defaults_true(self):
        """Test include_defaults=True behavior."""
        generator = BaseLammpsSetGenerator(settings={"units": "real"}, include_defaults=True)

        # Should include base settings and override with user settings
        base_settings = _BASE_LAMMPS_SETTINGS["periodic"]
        assert generator.settings.units == "real"  # User override
        assert generator.settings.dimension == base_settings["dimension"]  # From defaults
        assert generator.settings.pair_style == base_settings["pair_style"]  # From defaults

    def test_include_defaults_false(self):
        """Test include_defaults=False behavior."""
        generator = BaseLammpsSetGenerator(settings={"units": "real", "dimension": 2}, include_defaults=False)

        # Should only include user settings
        assert generator.settings.units == "real"
        assert generator.settings.dimension == 2
        # Other defaults should not be set
        assert not hasattr(generator.settings, "pair_style")

    def test_validate_params_true(self):
        """Test validate_params=True behavior."""
        # Valid settings should work
        generator = BaseLammpsSetGenerator(settings={"units": "metal", "atom_style": "full"}, validate_params=True)
        assert generator.settings.units == "metal"

    def test_validate_params_false(self):
        """Test validate_params=False behavior."""
        # Invalid settings should work when validation is disabled
        generator = BaseLammpsSetGenerator(settings={"units": "invalid_unit"}, validate_params=False)
        assert generator.settings.units == "invalid_unit"

    def test_update_settings_with_lammps_settings(self):
        """Test update_settings method with LammpsSettings object."""
        generator = BaseLammpsSetGenerator(settings={"units": "metal", "atom_style": "full"})

        updates = {"dimension": 3, "ensemble": "nvt"}
        generator.update_settings(updates)

        assert generator.settings.dimension == 3
        assert generator.settings.ensemble == "nvt"
        assert generator.settings.units == "metal"  # Should remain unchanged

    def test_update_settings_with_dict(self):
        """Test update_settings method with dictionary settings."""
        generator = BaseLammpsSetGenerator(settings={"units": "metal"})

        updates = {"atom_style": "full", "dimension": 3}
        generator.update_settings(updates)

        assert generator.settings.atom_style == "full"
        assert generator.settings.dimension == 3

    def test_update_settings_include_defaults_true(self):
        """Test update_settings with include_defaults=True."""
        generator = BaseLammpsSetGenerator(settings={"units": "real"}, include_defaults=False)

        updates = {"dimension": 3}
        generator.update_settings(updates, include_defaults=True)

        # Should include base settings
        base_settings = _BASE_LAMMPS_SETTINGS["periodic"]
        assert generator.settings.units == "real"  # Original setting
        assert generator.settings.dimension == 3  # Updated setting
        assert generator.settings.pair_style == base_settings["pair_style"]  # From defaults

    def test_update_settings_validate_params_true(self):
        """Test update_settings with validate_params=True."""
        generator = BaseLammpsSetGenerator(settings={"units": "metal"})

        # Valid updates should work
        updates = {"atom_style": "full", "ensemble": "nvt"}
        generator.update_settings(updates, validate_params=True)
        assert generator.settings.atom_style == "full"

    def test_update_settings_validate_params_false(self):
        """Test update_settings with validate_params=False."""
        generator = BaseLammpsSetGenerator(settings={"units": "metal"})

        # Invalid updates should work when validation is disabled
        updates = {"units": "invalid_unit"}
        generator.update_settings(updates, validate_params=False)
        assert generator.settings.units == "invalid_unit"

    def test_get_input_set_with_structure(self):
        """Test get_input_set method with Structure."""
        generator = BaseLammpsSetGenerator(settings={"ensemble": "nvt"})
        input_set = generator.get_input_set(self.structure)

        assert input_set is not None
        assert hasattr(input_set, "data")
        assert hasattr(input_set, "inputfile")

    def test_get_input_set_with_molecule(self):
        """Test get_input_set method with Molecule."""
        from pymatgen.core import Molecule

        molecule = Molecule(["H", "H"], [[0, 0, 0], [0, 0, 0.74]])
        generator = BaseLammpsSetGenerator(data_type="molecular", settings={"ensemble": "nvt"})
        input_set = generator.get_input_set(molecule)

        assert input_set is not None
        assert hasattr(input_set, "data")
        assert hasattr(input_set, "inputfile")

    def test_get_input_set_with_lammps_data(self):
        """Test get_input_set method with LammpsData."""
        lammps_data = LammpsData.from_structure(self.structure)
        generator = BaseLammpsSetGenerator(settings={"ensemble": "nvt"})
        input_set = generator.get_input_set(lammps_data)

        assert input_set is not None
        assert hasattr(input_set, "data")
        assert hasattr(input_set, "inputfile")

    def test_get_input_set_with_additional_data(self):
        """Test get_input_set method with additional data."""
        lammps_data = LammpsData.from_structure(self.structure)
        additional_data = LammpsData.from_structure(self.structure)

        generator = BaseLammpsSetGenerator(settings={"ensemble": "nvt"})
        input_set = generator.get_input_set(lammps_data, additional_data=additional_data)

        assert input_set is not None
        assert hasattr(input_set, "data")
        assert hasattr(input_set, "inputfile")

    def test_keep_stages_true(self):
        """Test keep_stages=True behavior."""
        generator = BaseLammpsSetGenerator(keep_stages=True, settings={"ensemble": "nvt"})
        assert generator.keep_stages is True
        assert generator.inputfile is not None

    def test_keep_stages_false(self):
        """Test keep_stages=False behavior."""
        generator = BaseLammpsSetGenerator(keep_stages=False, settings={"ensemble": "nvt"})
        assert generator.keep_stages is False
        assert generator.inputfile is not None

    def test_override_updates_false(self):
        """Test override_updates=False behavior (default)."""
        generator = BaseLammpsSetGenerator(settings={"ensemble": "nvt"})
        # This is the default behavior, so just test that it's set correctly
        assert hasattr(generator, "override_updates")  # This attribute should exist

    def test_custom_inputfile(self):
        """Test initialization with custom inputfile."""
        from pymatgen.io.lammps.inputs import LammpsInputFile

        custom_input = LammpsInputFile()
        generator = BaseLammpsSetGenerator(inputfile=custom_input)

        assert generator.inputfile is custom_input

    def test_custom_inputfile_string(self):
        """Test initialization with custom inputfile as string path."""
        generator = BaseLammpsSetGenerator(inputfile="custom.in")

        assert generator.inputfile == "custom.in"

    def test_custom_inputfile_path(self):
        """Test initialization with custom inputfile as Path object."""
        from pathlib import Path

        custom_path = Path("custom.in")
        generator = BaseLammpsSetGenerator(inputfile=custom_path)

        assert generator.inputfile == custom_path


class TestLammpsMinimization(unittest.TestCase):
    @classmethod
    def setup_class(cls):
        cls.filename = f"{TEST_DIR}/lgps.in"
        cls.cif = f"{TEST_DIR}/lgps.cif"
        cls.structure = Structure.from_file(cls.cif)

    def test_get_input_set(self):
        lmp_min = LammpsMinimization(keep_stages=False).get_input_set(self.structure)
        assert list(lmp_min.data.as_dict()) == list(
            LammpsData.from_structure(self.structure, atom_style="full").as_dict()
        )
        assert (
            lmp_min.data.as_dict()["atoms"].to_numpy()
            == LammpsData.from_structure(self.structure, atom_style="full").as_dict()["atoms"].to_numpy()
        ).all()
        """assert lmp_min.inputfile.stages == [
            {
                "stage_name": "Stage 1",
                "commands": [
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("boundary", "p p p"),
                    ("#", "2) System definition"),
                    ("read_data", "input.data"),
                    ("neigh_modify", "every 1 delay 0 check yes"),
                    ("#", "3) Simulation settings"),
                    ("pair_style", "$pair_style"),
                    ("bond_style", "$bond_style"),
                    ("angle_style", "$angle_style"),
                    ("dihedral_style", "$dihedral_style"),
                    ("improper_style", "$improper_style"),
                    ("include", "forcefield.lammps"),
                    ("extra_data", "$extra_data"),
                    ("#", "4) Energy minimization"),
                    ("thermo", "5"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                    ("min_style", "$min_style"),
                    ("fix", "1 all box/relax $psymm $start_pressure vmax 0.001"),
                    ("minimize", "$tol $tol $nsteps 10000000"),
                    ("#", "5) Write output data"),
                    ("write_data", "output.data"),
                    ("write_restart", "run.restart"),
                ],
            }
        ]

        lmp_min = LammpsMinimization(units="atomic", dimension=2, keep_stages=False).get_input_set(self.structure)
        assert list(lmp_min.data.as_dict()) == list(LammpsData.from_structure(self.structure).as_dict())
        assert (
            lmp_min.data.as_dict()["atoms"].to_numpy()
            == LammpsData.from_structure(self.structure).as_dict()["atoms"].to_numpy()
        ).all()
        assert lmp_min.inputfile.stages == [
            {
                "stage_name": "Stage 1",
                "commands": [
                    ("units", "atomic"),
                    ("atom_style", "full"),
                    ("dimension", "2"),
                    ("boundary", "p p p"),
                    ("#", "2) System definition"),
                    ("read_data", "system.data"),
                    ("neigh_modify", "every 1 delay 5 check yes"),
                    ("#", "3) Simulation settings"),
                    ("Unspecified", "force field!"),
                    ("#", "4) Energy minimization"),
                    ("thermo", "5"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 100000"),
                    ("#", "5) Write output data"),
                    ("write_data", "run.data"),
                    ("write_restart", "run.restart"),
                ],
            }
        ]

        lmp_min = LammpsMinimization(keep_stages=True).get_input_set(self.structure)
        assert lmp_min.inputfile.stages == [
            {
                "stage_name": "1) Initialization",
                "commands": [
                    ("units", "metal"),
                    ("atom_style", "full"),
                    ("dimension", "3"),
                    ("boundary", "p p p"),
                ],
            },
            {
                "stage_name": "2) System definition",
                "commands": [
                    ("read_data", "system.data"),
                    ("neigh_modify", "every 1 delay 5 check yes"),
                ],
            },
            {
                "stage_name": "3) Simulation settings",
                "commands": [("Unspecified", "force field!")],
            },
            {
                "stage_name": "4) Energy minimization",
                "commands": [
                    ("thermo", "5"),
                    ("thermo_style", "custom step lx ly lz press pxx pyy pzz pe"),
                    ("dump", "dmp all atom 5 run.dump"),
                ],
            },
            {
                "stage_name": "Stage 5",
                "commands": [
                    ("min_style", "cg"),
                    ("fix", "1 all box/relax iso 0.0 vmax 0.001"),
                    ("minimize", "1.0e-16 1.0e-16 5000 100000"),
                ],
            },
            {
                "stage_name": "5) Write output data",
                "commands": [
                    ("write_data", "run.data"),
                    ("write_restart", "run.restart"),
                ],
            },
        ]

        assert lmp_min.inputfile.nstages == 6
        assert lmp_min.inputfile.ncomments == 0
"""
