Change log
==========

v2017.9.1
---------
* Massive refactoring of PhaseDiagram. Now, PDAnalyzer is completely defunct
  and all analysis is carried out within PhaseDiagram itself, e.g.,
  pd.get_e_above_hull as opposed to PDAnalyzer(pd).get_e_above_hull.
* Refactoring of structure prediction. Now in
  pymatgen.analysis.structure_prediction.
* New core Spectrum object and associated pymatgen.vis.plotters.SpectrumPlotter.
* Parsing energies from gen_scfman module in Qchem 5 (Brandon Wood)
* Improvements to LAMMPSData, vasp IO.
