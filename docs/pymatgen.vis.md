---
layout: default
title: pymatgen.vis.md
nav_exclude: true
---

# pymatgen.vis package

The vis package implements various visualization tools. For example, a VTK
Structure viewer.



* [pymatgen.vis.plotters module](pymatgen.vis.plotters.md)


    * [`SpectrumPlotter`](pymatgen.vis.plotters.md#pymatgen.vis.plotters.SpectrumPlotter)


        * [`SpectrumPlotter.add_spectra()`](pymatgen.vis.plotters.md#pymatgen.vis.plotters.SpectrumPlotter.add_spectra)


        * [`SpectrumPlotter.add_spectrum()`](pymatgen.vis.plotters.md#pymatgen.vis.plotters.SpectrumPlotter.add_spectrum)


        * [`SpectrumPlotter.get_plot()`](pymatgen.vis.plotters.md#pymatgen.vis.plotters.SpectrumPlotter.get_plot)


        * [`SpectrumPlotter.save_plot()`](pymatgen.vis.plotters.md#pymatgen.vis.plotters.SpectrumPlotter.save_plot)


        * [`SpectrumPlotter.show()`](pymatgen.vis.plotters.md#pymatgen.vis.plotters.SpectrumPlotter.show)


* [pymatgen.vis.structure_chemview module](pymatgen.vis.structure_chemview.md)


    * [`quick_view()`](pymatgen.vis.structure_chemview.md#pymatgen.vis.structure_chemview.quick_view)


* [pymatgen.vis.structure_vtk module](pymatgen.vis.structure_vtk.md)


    * [`MultiStructuresInteractorStyle`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresInteractorStyle)


        * [`MultiStructuresInteractorStyle.keyPressEvent()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresInteractorStyle.keyPressEvent)


    * [`MultiStructuresVis`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis)


        * [`MultiStructuresVis.DEFAULT_ANIMATED_MOVIE_OPTIONS`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.DEFAULT_ANIMATED_MOVIE_OPTIONS)


        * [`MultiStructuresVis.apply_tags()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.apply_tags)


        * [`MultiStructuresVis.display_help()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.display_help)


        * [`MultiStructuresVis.display_info()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.display_info)


        * [`MultiStructuresVis.display_warning()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.display_warning)


        * [`MultiStructuresVis.erase_info()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.erase_info)


        * [`MultiStructuresVis.erase_warning()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.erase_warning)


        * [`MultiStructuresVis.set_animated_movie_options()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.set_animated_movie_options)


        * [`MultiStructuresVis.set_structure()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.set_structure)


        * [`MultiStructuresVis.set_structures()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.MultiStructuresVis.set_structures)


    * [`StructureInteractorStyle`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureInteractorStyle)


        * [`StructureInteractorStyle.keyPressEvent()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureInteractorStyle.keyPressEvent)


        * [`StructureInteractorStyle.leftButtonPressEvent()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureInteractorStyle.leftButtonPressEvent)


        * [`StructureInteractorStyle.leftButtonReleaseEvent()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureInteractorStyle.leftButtonReleaseEvent)


        * [`StructureInteractorStyle.mouseMoveEvent()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureInteractorStyle.mouseMoveEvent)


    * [`StructureVis`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis)


        * [`StructureVis.add_bonds()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_bonds)


        * [`StructureVis.add_edges()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_edges)


        * [`StructureVis.add_faces()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_faces)


        * [`StructureVis.add_line()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_line)


        * [`StructureVis.add_partial_sphere()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_partial_sphere)


        * [`StructureVis.add_picker()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_picker)


        * [`StructureVis.add_picker_fixed()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_picker_fixed)


        * [`StructureVis.add_polyhedron()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_polyhedron)


        * [`StructureVis.add_site()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_site)


        * [`StructureVis.add_text()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_text)


        * [`StructureVis.add_triangle()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.add_triangle)


        * [`StructureVis.display_help()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.display_help)


        * [`StructureVis.orthongonalize_structure()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.orthongonalize_structure)


        * [`StructureVis.redraw()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.redraw)


        * [`StructureVis.rotate_view()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.rotate_view)


        * [`StructureVis.set_structure()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.set_structure)


        * [`StructureVis.show()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.show)


        * [`StructureVis.write_image()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.write_image)


        * [`StructureVis.zoom()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.StructureVis.zoom)


    * [`make_movie()`](pymatgen.vis.structure_vtk.md#pymatgen.vis.structure_vtk.make_movie)