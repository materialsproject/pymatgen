Change log
==========

v2019.3.13
----------
* Streamlined Site, PeriodicSite, Molecule and Structure code by abandoning
  immutability for Site and PeriodicSite.
* VaspInput class now supports a run_vasp method, which can be used to code
  runnable python scripts for running simple calculations (custodian still
  recommended for more complex calculations.). For example, the following is a
  kpoint convergence script that can be submitted in a queue.

```python
from pymatgen import MPRester
from pymatgen.io.vasp.sets import MPRelaxSet


VASP_CMD = ["mpirun", "-machinefile", "$PBS_NODEFILE", "-np", "16", "vasp"]


def main():
    mpr = MPRester()
    structure = mpr.get_structure("Li2O")
    for k_dens in [100, 200, 400, 800]:
        vis = MPRelaxSet(structure, 
            user_kpoints_settings={"reciprocal_density": k_dens})
        kpoints = vi["KPOINTS"].kpts[0][0]
        d = "Li2O_kpoints_%d" % kpoints
        vi.run_vasp(d, vasp_cmd=VASP_CMD)
        # Use the final structure as the new initial structure to speed up calculations.
        structure = Vasprun("%s/vasprun.xml" % d).final_structure


if __name__ == "__main__":
    main()
```

* Many pymatgen from_file methods now support pathlib.Path as well as strings.
* Misc bug fixes.

