MKDIR tmpconda
cd tmpconda
FOR %%A IN (latexcodec tabulate monty pybtex palettable spglib pydispatcher pymatgen) DO conda skeleton pypi %%A
FOR %%A IN (latexcodec tabulate monty pybtex palettable spglib pydispatcher pymatgen) DO conda build %%A
FOR %%A IN (latexcodec tabulate monty pybtex palettable spglib pydispatcher pymatgen) DO anaconda upload %HOMEPATH%\Miniconda3\win-64\%AA-*py35*.tar.bz2
cd ..
RMDIR /s /q tmpconda
