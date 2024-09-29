# JDFTx io

This folder contains all class objects and methods used for io operations with JDFTX.

This includes:

- JDFTXInfile:
  - Initializable from pre-existing JDFTx in files or dictionaries
  - Modifiable (acts as dictionary)
  - Writes new in files
  - Extracts Structure objects
- JDFTXOutfile:
  - Parses JDFTx out file
  - Contains all typically relevant output variables from a JDFTx geometric optimization or single-point calculation.
  - Contains hierarchy of class objects (each contained by the former) for storing data at the following call frequencies.
    - JDFTXOutfileSlice: Per call of JDFTx executable
      - One "slice" contains all data output in a single call of JDFTx,
          as broken up in the out file by the "**** JDFTx" flag
      - JOutStructures: Per call of JDFTx executable (same frequency as parent)
        - List of JOutStructure, used in building Trajectory objects by JDFTXOutfile (does not inherit Trajectory or Structure)
      - JOutStructure: Per geometric optimization update (only one for single-point)
        - Inheritor of Structure object (contains structural data), also contains electronic minimization data (see below) and convergence data relevant to the geometric optimization (forces and Wolfeminimization variables)
      - JEiters: Per geometric optimization update (same frequency as parent)
        - List of JEiter as well as convergence data relevant to electronic optimization.
      - JEiter: Per SCF update
            - Contains all electronic data logged in out file at SCF update frequency.

This folder is currently missing:

- Broader output parsing
  - Future versions will still have JDFTXOutfile and all subclasses. However, this class object will likely not be called directly, and will be a class property of a broader JDFTXOutput class capable of handling dumped output files (ie electronic density arrays or DOS text files) in tandem to the out file.
- File organization
  - Future versions will streamline the user-end interface by having all output-related methods imported through the outputs.py module and input-related methods imported through the inputs.py module. Class objects may be moved to these modules from their standalone modules to keep the directory tidy.
- Sets
  - Common input sets are currently missing from the sets module.
