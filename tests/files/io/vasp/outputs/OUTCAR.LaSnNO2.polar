 vasp.5.4.4.18Apr17-6-g9f103f2a35 (build Feb 12 2018 00:57:43) complex          
  
 executed on             LinuxIFC date 2019.05.15  18:05:34
 running on   64 total cores
 distrk:  each k-point on   64 cores,    1 groups
 distr:  one band on NCORES_PER_BAND=   1 cores,   64 groups


--------------------------------------------------------------------------------------------------------


 INCAR:
 POTCAR:   PAW_PBE La 06Sep2000                   
 POTCAR:   PAW_PBE N 08Apr2002                    
 POTCAR:   PAW_PBE O 08Apr2002                    
 POTCAR:   PAW_PBE Sn_d 06Sep2000                 

 ----------------------------------------------------------------------------- 
|                                                                             |
|           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
|           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
|           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
|           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
|           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
|           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
|                                                                             |
|      For optimal performance we recommend to set                            |
|        NCORE= 4 - approx SQRT( number of cores)                             |
|      NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE).     |
|      This setting can  greatly improve the performance of VASP for DFT.     |
|      The default,   NCORE=1            might be grossly inefficient         |
|      on modern multi-core architectures or massively parallel machines.     |
|      Do your own testing !!!!                                               |
|      Unfortunately you need to use the default for GW and RPA calculations. |
|      (for HF NCORE is supported but not extensively tested yet)             |
|                                                                             |
 ----------------------------------------------------------------------------- 

 POTCAR:   PAW_PBE La 06Sep2000                   
   VRHFIN =La : [core=Kr4d]                                                     
   LEXCH  = PE                                                                  
   EATOM  =   865.2204 eV,   63.5919 Ry                                         
                                                                                
   TITEL  = PAW_PBE La 06Sep2000                                                
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =    2.300    partial core radius                                     
   POMASS =  138.900; ZVAL   =   11.000    mass and valenz                      
   RCORE  =    2.800    outmost cutoff radius                                   
   RWIGS  =    2.900; RWIGS  =    1.535    wigner-seitz radius (au A)           
   ENMAX  =  219.313; ENMIN  =  164.485 eV                                      
   RCLOC  =    1.601    cutoff for local pot                                    
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  583.575                                                            
   DEXC   =     .000                                                            
   RMAX   =    3.113    core radius for proj-oper                               
   RAUG   =    1.300    factor for augmentation sphere                          
   RDEP   =    2.802    radius for radial grids                                 
   QCUT   =   -4.015; QGAM   =    8.030    optimization parameters              
                                                                                
   Description                                                                  
     l     E      TYP  RCUT    TYP  RCUT                                        
     0   .000     23  2.800                                                     
     0   .000     23  2.800                                                     
     1   .000     23  2.500                                                     
     1  2.000     23  2.500                                                     
     2   .000     23  2.500                                                     
     2   .000     23  2.500                                                     
     3   .000     23  2.800                                                     
     3   .000     23  2.800                                                     
  local pseudopotential read in
  partial core-charges read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           2  read in
    real space projection operators read in
  non local Contribution for L=           2  read in
    real space projection operators read in
  non local Contribution for L=           3  read in
    real space projection operators read in
  non local Contribution for L=           3  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           8
   number of lm-projection operators is LMMAX =          32
 
 POTCAR:   PAW_PBE N 08Apr2002                    
   VRHFIN =N: s2p3                                                              
   LEXCH  = PE                                                                  
   EATOM  =   264.5486 eV,   19.4438 Ry                                         
                                                                                
   TITEL  = PAW_PBE N 08Apr2002                                                 
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =     .000    partial core radius                                     
   POMASS =   14.001; ZVAL   =    5.000    mass and valenz                      
   RCORE  =    1.500    outmost cutoff radius                                   
   RWIGS  =    1.400; RWIGS  =     .741    wigner-seitz radius (au A)           
   ENMAX  =  400.000; ENMIN  =  300.000 eV                                      
   ICORE  =        2    local potential                                         
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  627.112                                                            
   DEXC   =     .000                                                            
   RMAX   =    2.247    core radius for proj-oper                               
   RAUG   =    1.300    factor for augmentation sphere                          
   RDEP   =    1.514    radius for radial grids                                 
   QCUT   =   -5.562; QGAM   =   11.124    optimization parameters              
                                                                                
   Description                                                                  
     l     E      TYP  RCUT    TYP  RCUT                                        
     0   .000     23  1.200                                                     
     0   .000     23  1.200                                                     
     1   .000     23  1.500                                                     
     1   .700     23  1.500                                                     
     2   .000      7  1.500                                                     
  local pseudopotential read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           4
   number of lm-projection operators is LMMAX =           8
 
 POTCAR:   PAW_PBE O 08Apr2002                    
   VRHFIN =O: s2p4                                                              
   LEXCH  = PE                                                                  
   EATOM  =   432.3788 eV,   31.7789 Ry                                         
                                                                                
   TITEL  = PAW_PBE O 08Apr2002                                                 
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        0    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =     .000    partial core radius                                     
   POMASS =   16.000; ZVAL   =    6.000    mass and valenz                      
   RCORE  =    1.520    outmost cutoff radius                                   
   RWIGS  =    1.550; RWIGS  =     .820    wigner-seitz radius (au A)           
   ENMAX  =  400.000; ENMIN  =  300.000 eV                                      
   ICORE  =        2    local potential                                         
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  605.392                                                            
   DEXC   =     .000                                                            
   RMAX   =    2.264    core radius for proj-oper                               
   RAUG   =    1.300    factor for augmentation sphere                          
   RDEP   =    1.550    radius for radial grids                                 
   QCUT   =   -5.520; QGAM   =   11.041    optimization parameters              
                                                                                
   Description                                                                  
     l     E      TYP  RCUT    TYP  RCUT                                        
     0   .000     23  1.200                                                     
     0  -.700     23  1.200                                                     
     1   .000     23  1.520                                                     
     1   .600     23  1.520                                                     
     2   .000      7  1.500                                                     
  local pseudopotential read in
  atomic valenz-charges read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           4
   number of lm-projection operators is LMMAX =           8
 
 POTCAR:   PAW_PBE Sn_d 06Sep2000                 
   VRHFIN =Sn: s2p2                                                             
   LEXCH  = PE                                                                  
   EATOM  =  1893.0782 eV,  139.1373 Ry                                         
                                                                                
   TITEL  = PAW_PBE Sn_d 06Sep2000                                              
   LULTRA =        F    use ultrasoft PP ?                                      
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no                           
   RPACOR =    2.300    partial core radius                                     
   POMASS =  118.710; ZVAL   =   14.000    mass and valenz                      
   RCORE  =    2.500    outmost cutoff radius                                   
   RWIGS  =    2.960; RWIGS  =    1.566    wigner-seitz radius (au A)           
   ENMAX  =  241.090; ENMIN  =  180.817 eV                                      
   ICORE  =        3    local potential                                         
   LCOR   =        T    correct aug charges                                     
   LPAW   =        T    paw PP                                                  
   EAUG   =  439.825                                                            
   DEXC   =    -.001                                                            
   RMAX   =    2.969    core radius for proj-oper                               
   RAUG   =    1.300    factor for augmentation sphere                          
   RDEP   =    2.659    radius for radial grids                                 
   QCUT   =   -4.209; QGAM   =    8.419    optimization parameters              
                                                                                
   Description                                                                  
     l     E      TYP  RCUT    TYP  RCUT                                        
     2   .000     23  2.500                                                     
     2  -.400     23  2.500                                                     
     0   .000     23  2.500                                                     
     0  1.000     23  2.500                                                     
     1   .000     23  2.500                                                     
     1   .000     23  2.500                                                     
     3  -.100      7  2.500                                                     
  local pseudopotential read in
  partial core-charges read in
  atomic valenz-charges read in
  non local Contribution for L=           2  read in
    real space projection operators read in
  non local Contribution for L=           2  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           0  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
  non local Contribution for L=           1  read in
    real space projection operators read in
    PAW grid and wavefunctions read in
 
   number of l-projection  operators is LMAX  =           6
   number of lm-projection operators is LMMAX =          18
 
 Optimization of the real space projectors (new method)

 maximal supplied QI-value         = 13.42
 optimisation between [QCUT,QGAM] = [ 11.67, 23.35] = [ 38.16,152.62] Ry 
 Optimized for a Real-space Cutoff    1.63 Angstroem

   l    n(q)    QCUT    max X(q) W(low)/X(q) W(high)/X(q)  e(spline) 
   0     11    11.673     5.065    0.51E-05    0.90E-06    0.19E-07
   0     11    11.673    15.079    0.39E-04    0.13E-04    0.19E-06
   1     11    11.673     3.052    0.28E-05    0.12E-05    0.45E-07
   1     11    11.673     5.332    0.88E-05    0.13E-04    0.19E-06
   2     10    11.673     4.118    0.13E-04    0.64E-05    0.91E-07
   2     10    11.673    34.547    0.61E-04    0.33E-04    0.17E-06
   3     10    11.673     6.436    0.16E-04    0.45E-04    0.11E-06
   3     10    11.673    18.530    0.17E-03    0.14E-03    0.37E-06
 Optimization of the real space projectors (new method)

 maximal supplied QI-value         = 25.13
 optimisation between [QCUT,QGAM] = [ 11.56, 23.12] = [ 37.43,149.72] Ry 
 Optimized for a Real-space Cutoff    1.40 Angstroem

   l    n(q)    QCUT    max X(q) W(low)/X(q) W(high)/X(q)  e(spline) 
   0     10    11.561     4.511    0.22E-04    0.14E-04    0.74E-07
   0     10    11.561    37.799    0.45E-04    0.77E-04    0.61E-06
   1      9    11.561     2.478    0.97E-05    0.16E-04    0.10E-06
   1      9    11.561    11.620    0.23E-03    0.13E-03    0.12E-05
 Optimization of the real space projectors (new method)

 maximal supplied QI-value         = 24.76
 optimisation between [QCUT,QGAM] = [ 11.64, 23.27] = [ 37.91,151.63] Ry 
 Optimized for a Real-space Cutoff    1.08 Angstroem

   l    n(q)    QCUT    max X(q) W(low)/X(q) W(high)/X(q)  e(spline) 
   0      7    11.635     4.192    0.54E-04    0.80E-04    0.26E-07
   0      7    11.635     8.473    0.70E-04    0.40E-03    0.21E-06
   1      7    11.635     2.474    0.39E-04    0.14E-03    0.23E-07
   1      7    11.635     3.912    0.25E-03    0.37E-03    0.22E-06
 Optimization of the real space projectors (new method)

 maximal supplied QI-value         = 15.12
 optimisation between [QCUT,QGAM] = [ 11.64, 23.28] = [ 37.95,151.78] Ry 
 Optimized for a Real-space Cutoff    1.36 Angstroem

   l    n(q)    QCUT    max X(q) W(low)/X(q) W(high)/X(q)  e(spline) 
   2      9    11.641     6.369    0.11E-03    0.18E-03    0.58E-07
   2      9    11.641    10.755    0.17E-03    0.18E-03    0.15E-06
   0      9    11.641     8.594    0.12E-04    0.73E-06    0.48E-08
   0      9    11.641    13.234    0.13E-03    0.43E-04    0.95E-07
   1      9    11.641     7.249    0.80E-05    0.19E-04    0.11E-07
   1      9    11.641    72.638    0.24E-03    0.18E-03    0.20E-06
 PAW_PBE La 06Sep2000                   :
 energy of atom  1       EATOM= -865.2204
 kinetic energy error for atom=    0.0009 (will be added to EATOM!!)
 PAW_PBE N 08Apr2002                    :
 energy of atom  2       EATOM= -264.5486
 kinetic energy error for atom=    0.0125 (will be added to EATOM!!)
 PAW_PBE O 08Apr2002                    :
 energy of atom  3       EATOM= -432.3788
 kinetic energy error for atom=    0.0208 (will be added to EATOM!!)
 PAW_PBE Sn_d 06Sep2000                 :
 energy of atom  4       EATOM=-1893.0782
 kinetic energy error for atom=    0.0026 (will be added to EATOM!!)
 
 
 POSCAR: La6 Sn6 N6 O12                          
  positions in direct lattice
  velocities in cartesian coordinates
 exchange correlation table for  LEXCH =        8
   RHO(1)=    0.500       N(1)  =     2000
   RHO(2)=  100.500       N(2)  =     4000
 


--------------------------------------------------------------------------------------------------------


 ion  position               nearest neighbor table
   1  0.000  0.000  0.269-  14 2.42  15 2.42  13 2.42  20 2.46  21 2.46  19 2.46   7 2.72  29 3.68
                            30 3.68  28 3.68
   2  0.000  0.000  0.769-  17 2.42  18 2.42  16 2.42  23 2.46  24 2.46  22 2.46   8 2.72  26 3.68
                            27 3.68  25 3.68
   3  0.333  0.667  0.238-  13 2.42  14 2.42  15 2.42  20 2.45  19 2.45  21 2.45
   4  0.667  0.333  0.738-  16 2.42  17 2.42  18 2.42  23 2.45  22 2.45  24 2.45
   5  0.333  0.667  0.738-  17 2.42  16 2.42  18 2.42  22 2.45  23 2.45  24 2.45
   6  0.667  0.333  0.238-  14 2.42  13 2.42  15 2.42  19 2.45  20 2.45  21 2.45
   7  0.000  0.000  0.481-  29 2.18  30 2.18  28 2.18   1 2.72
   8  0.000  0.000  0.981-  26 2.18  27 2.18  25 2.18   2 2.72
   9  0.333  0.667  0.018-  25 2.18  26 2.18  27 2.18
  10  0.667  0.333  0.518-  28 2.18  29 2.18  30 2.18
  11  0.333  0.667  0.518-  29 2.18  28 2.18  30 2.18
  12  0.667  0.333  0.018-  26 2.18  25 2.18  27 2.18
  13  0.313  0.000  0.166-  25 2.13   1 2.42   3 2.42   6 2.42
  14  0.000  0.313  0.166-  26 2.13   1 2.42   6 2.42   3 2.42
  15  0.687  0.687  0.166-  27 2.13   1 2.42   6 2.42   3 2.42
  16  0.687  0.000  0.666-  28 2.13   2 2.42   4 2.42   5 2.42
  17  0.000  0.687  0.666-  29 2.13   2 2.42   5 2.42   4 2.42
  18  0.313  0.313  0.666-  30 2.13   2 2.42   5 2.42   4 2.42
  19  0.646  0.000  0.336-  28 2.12   6 2.45   3 2.45   1 2.46
  20  0.000  0.646  0.336-  29 2.12   3 2.45   6 2.45   1 2.46
  21  0.354  0.354  0.336-  30 2.12   3 2.45   6 2.45   1 2.46
  22  0.354  0.000  0.836-  25 2.12   5 2.45   4 2.45   2 2.46
  23  0.000  0.354  0.836-  26 2.12   4 2.45   5 2.45   2 2.46
  24  0.646  0.646  0.836-  27 2.12   4 2.45   5 2.45   2 2.46
  25  0.333  0.000  0.001-  22 2.12  13 2.13   9 2.18  12 2.18   8 2.18   2 3.68  26 3.76  27 3.76
                            26 3.76  27 3.76  26 3.76  27 3.76
  26  0.000  0.333  0.001-  23 2.12  14 2.13  12 2.18   9 2.18   8 2.18   2 3.68  27 3.76  25 3.76
                            25 3.76  27 3.76  27 3.76  25 3.76
  27  0.667  0.667  0.001-  24 2.12  15 2.13  12 2.18   9 2.18   8 2.18   2 3.68  26 3.76  25 3.76
                            25 3.76  26 3.76  26 3.76  25 3.76
  28  0.667  0.000  0.501-  19 2.12  16 2.13  10 2.18  11 2.18   7 2.18   1 3.68  29 3.76  30 3.76
                            29 3.76  30 3.76  30 3.76  29 3.76
  29  0.000  0.667  0.501-  20 2.12  17 2.13  11 2.18  10 2.18   7 2.18   1 3.68  30 3.76  28 3.76
                            28 3.76  30 3.76  30 3.76  28 3.76
  30  0.333  0.333  0.501-  21 2.12  18 2.13  11 2.18  10 2.18   7 2.18   1 3.68  29 3.76  28 3.76
                            28 3.76  29 3.76  29 3.76  28 3.76
 
  LATTYP: Found a hexagonal cell.
 ALAT       =     6.5064997491
 C/A-ratio  =     1.9711520010
  
  Lattice vectors:
  
 A1 = (   6.5065000000,   0.0000000000,   0.0000000000)
 A2 = (  -3.2532500000,   5.6347940000,   0.0000000000)
 A3 = (   0.0000000000,   0.0000000000,  12.8253000000)


Analysis of symmetry for initial positions (statically):
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 hexagonal supercell.


 Subroutine GETGRP returns: Found 12 space group operations
 (whereof  6 operations were pure point group operations)
 out of a pool of 24 trial point group operations.


The static configuration has the point symmetry C_3v.
 The point group associated with its full space group is C_6v.


Analysis of symmetry for dynamics (positions and initial velocities):
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 hexagonal supercell.


 Subroutine GETGRP returns: Found 12 space group operations
 (whereof  6 operations were pure point group operations)
 out of a pool of 24 trial point group operations.


The dynamic configuration has the point symmetry C_3v.
 The point group associated with its full space group is C_6v.


Analysis of structural, dynamic, and magnetic symmetry:
=====================================================================
 Subroutine PRICEL returns:
 Original cell was already a primitive cell.
 

 Routine SETGRP: Setting up the symmetry group for a 
 hexagonal supercell.


 Subroutine GETGRP returns: Found 12 space group operations
 (whereof  6 operations were pure point group operations)
 out of a pool of 24 trial point group operations.


The magnetic configuration has the point symmetry C_3v.
 The point group associated with its full space group is C_6v.


 Subroutine INISYM returns: Found 12 space group operations
 (whereof  6 operations are pure point group operations),
 and found     1 'primitive' translations

 
 
 KPOINTS: pymatgen v1.1 with grid density = 1583 /

Automatic generation of k-mesh.
Space group operators:
 irot       det(A)        alpha          n_x          n_y          n_z        tau_x        tau_y        tau_z
    1     1.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000
    2     1.000000   120.000000     0.000000     0.000000    -1.000000     0.000000     0.000000     0.000000
    3     1.000000   120.000000     0.000000     0.000000     1.000000     0.000000     0.000000     0.000000
    4    -1.000000   180.000000     0.866025    -0.500000     0.000000     0.000000     0.000000     0.000000
    5    -1.000000   180.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000
    6    -1.000000   180.000000     0.866025     0.500000     0.000000     0.000000     0.000000     0.000000
    7     1.000000    60.000000     0.000000     0.000000    -1.000000     0.000001     0.000000     0.500000
    8     1.000000   180.000000     0.000000     0.000000     1.000000     0.000000     0.000000     0.500000
    9     1.000000    60.000000     0.000000     0.000000     1.000000     0.000000     0.000000     0.500000
   10    -1.000000   180.000000     1.000000     0.000000     0.000000     0.000001     0.000000     0.500000
   11    -1.000000   180.000000     0.500000    -0.866025     0.000000     0.000000     0.000000     0.500000
   12    -1.000000   180.000000    -0.500000    -0.866025     0.000000     0.000000     0.000000     0.500000
 
 Subroutine IBZKPT returns following result:
 ===========================================
 
 Found      8 irreducible k-points:
 
 Following reciprocal coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
  0.250000  0.000000  0.000000      6.000000
  0.500000  0.000000  0.000000      3.000000
  0.250000  0.250000  0.000000      6.000000
  0.000000  0.000000  0.500000      1.000000
  0.250000  0.000000  0.500000      6.000000
  0.500000  0.000000  0.500000      3.000000
  0.250000  0.250000  0.500000      6.000000
 
 Following cartesian coordinates:
            Coordinates               Weight
  0.000000  0.000000  0.000000      1.000000
  0.038423  0.022184  0.000000      6.000000
  0.076846  0.044367  0.000000      3.000000
  0.038423  0.066551  0.000000      6.000000
  0.000000  0.000000  0.038985      1.000000
  0.038423  0.022184  0.038985      6.000000
  0.076846  0.044367  0.038985      3.000000
  0.038423  0.066551  0.038985      6.000000
 
 TETIRR: Found     36 inequivalent tetrahedra from      192
 
 Subroutine IBZKPT_HF returns following result:
 ==============================================
 
 Found     32 k-points in 1st BZ
 the following     32 k-points will be used (e.g. in the exchange kernel)
 Following reciprocal coordinates:   # in IRBZ
  0.000000  0.000000  0.000000    0.03125000   1 t-inv F
  0.250000  0.000000  0.000000    0.03125000   2 t-inv F
  0.500000  0.000000  0.000000    0.03125000   3 t-inv F
  0.250000  0.250000  0.000000    0.03125000   4 t-inv F
  0.000000  0.000000  0.500000    0.03125000   5 t-inv F
  0.250000  0.000000  0.500000    0.03125000   6 t-inv F
  0.500000  0.000000  0.500000    0.03125000   7 t-inv F
  0.250000  0.250000  0.500000    0.03125000   8 t-inv F
 -0.250000  0.250000  0.000000    0.03125000   2 t-inv F
  0.000000 -0.250000  0.000000    0.03125000   2 t-inv F
  0.000000  0.250000  0.000000    0.03125000   2 t-inv F
  0.250000 -0.250000  0.000000    0.03125000   2 t-inv F
 -0.250000  0.000000  0.000000    0.03125000   2 t-inv F
 -0.500000  0.500000  0.000000    0.03125000   3 t-inv F
  0.000000 -0.500000  0.000000    0.03125000   3 t-inv F
 -0.500000  0.250000  0.000000    0.03125000   4 t-inv F
  0.250000 -0.500000  0.000000    0.03125000   4 t-inv F
 -0.250000  0.500000  0.000000    0.03125000   4 t-inv F
 -0.250000 -0.250000  0.000000    0.03125000   4 t-inv F
  0.500000 -0.250000  0.000000    0.03125000   4 t-inv F
 -0.250000  0.250000  0.500000    0.03125000   6 t-inv F
  0.000000 -0.250000  0.500000    0.03125000   6 t-inv F
  0.000000  0.250000  0.500000    0.03125000   6 t-inv F
  0.250000 -0.250000  0.500000    0.03125000   6 t-inv F
 -0.250000  0.000000  0.500000    0.03125000   6 t-inv F
 -0.500000  0.500000  0.500000    0.03125000   7 t-inv F
  0.000000 -0.500000  0.500000    0.03125000   7 t-inv F
 -0.500000  0.250000  0.500000    0.03125000   8 t-inv F
  0.250000 -0.500000  0.500000    0.03125000   8 t-inv F
 -0.250000  0.500000  0.500000    0.03125000   8 t-inv F
 -0.250000 -0.250000  0.500000    0.03125000   8 t-inv F
  0.500000 -0.250000  0.500000    0.03125000   8 t-inv F


--------------------------------------------------------------------------------------------------------




 Dimension of arrays:
   k-points           NKPTS =      8   k-points in BZ     NKDIM =      8   number of bands    NBANDS=    192
   number of dos      NEDOS =   3000   number of ions     NIONS =     30
   non local maximal  LDIM  =      8   non local SUM 2l+1 LMDIM =     32
   total plane-waves  NPLWV = 221184
   max r-space proj   IRMAX =   8539   max aug-charges    IRDMAX=  24776
   dimension x,y,z NGX =    48 NGY =   48 NGZ =   96
   dimension x,y,z NGXF=    96 NGYF=   96 NGZF=  192
   support grid    NGXF=    96 NGYF=   96 NGZF=  192
   ions per type =               6   6  12   6
   NGX,Y,Z   is equivalent  to a cutoff of  12.26, 12.26, 12.44 a.u.
   NGXF,Y,Z  is equivalent  to a cutoff of  24.53, 24.53, 24.89 a.u.

 SYSTEM =  unknown system                          
 POSCAR =  La6 Sn6 N6 O12                          

 Startparameter for this run:
   NWRITE =      2    write-flag & timer
   PREC   = accura    normal or accurate (medium, high low for compatibility)
   ISTART =      0    job   : 0-new  1-cont  2-samecut
   ICHARG =      2    charge: 1-file 2-atom 10-const
   ISPIN  =      2    spin polarized calculation?
   LNONCOLLINEAR =      F non collinear calculations
   LSORBIT =      F    spin-orbit coupling
   INIWAV =      1    electr: 0-lowe 1-rand  2-diag
   LASPH  =      F    aspherical Exc in radial PAW
   METAGGA=      F    non-selfconsistent MetaGGA calc.

 Electronic Relaxation 1
   ENCUT  =  520.0 eV  38.22 Ry    6.18 a.u.  12.10 12.10 23.85*2*pi/ulx,y,z
   ENINI  =  520.0     initial cutoff
   ENAUG  =  627.1 eV  augmentation charge cutoff
   NELM   =    100;   NELMIN=  2; NELMDL=  0     # of ELM steps 
   EDIFF  = 0.1E-04   stopping-criterion for ELM
   LREAL  =      T    real-space projection
   NLSPLINE    = F    spline interpolate recip. space projectors
   LCOMPAT=      F    compatible to vasp.4.4
   GGA_COMPAT  = T    GGA compatible to vasp.4.4-vasp.4.6
   LMAXPAW     = -100 max onsite density
   LMAXMIX     =    2 max onsite mixed and CHGCAR
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   ROPT   =   -0.00025  -0.00025  -0.00025  -0.00025
 Ionic relaxation
   EDIFFG = 0.1E-03   stopping-criterion for IOM
   NSW    =      0    number of steps for IOM
   NBLOCK =      1;   KBLOCK =      1    inner block; outer block 
   IBRION =     -1    ionic relax: 0-MD 1-quasi-New 2-CG
   NFREE  =      0    steps in history (QN), initial steepest desc. (CG)
   ISIF   =      3    stress and relaxation
   IWAVPR =     10    prediction:  0-non 1-charg 2-wave 3-comb
   ISYM   =      2    0-nonsym 1-usesym 2-fastsym
   LCORR  =      T    Harris-Foulkes like correction to forces

   POTIM  = 0.5000    time-step for ionic-motion
   TEIN   =    0.0    initial temperature
   TEBEG  =    0.0;   TEEND  =   0.0 temperature during run
   SMASS  =  -3.00    Nose mass-parameter (am)
   estimated Nose-frequenzy (Omega)   =  0.10E-29 period in steps =****** mass=  -0.967E-27a.u.
   SCALEE = 1.0000    scale energy and forces
   NPACO  =    256;   APACO  = 16.0  distance and # of slots for P.C.
   PSTRESS=    0.0 pullay stress

  Mass of Ions in am
   POMASS = 138.90 14.00 16.00118.71
  Ionic Valenz
   ZVAL   =  11.00  5.00  6.00 14.00
  Atomic Wigner-Seitz radii
   RWIGS  =  -1.00 -1.00 -1.00 -1.00
  virtual crystal weights 
   VCA    =   1.00  1.00  1.00  1.00
   NELECT =     252.0000    total number of electrons
   NUPDOWN=      -1.0000    fix difference up-down

 DOS related values:
   EMIN   =  10.00;   EMAX   =-10.00  energy-range for DOS
   EFERMI =   0.00
   ISMEAR =    -5;   SIGMA  =   0.05  broadening in eV -4-tet -1-fermi 0-gaus

 Electronic relaxation 2 (details)
   IALGO  =     38    algorithm
   LDIAG  =      T    sub-space diagonalisation (order eigenvalues)
   LSUBROT=      F    optimize rotation matrix (better conditioning)
   TURBO    =      0    0=normal 1=particle mesh
   IRESTART =      0    0=no restart 2=restart with 2 vectors
   NREBOOT  =      0    no. of reboots
   NMIN     =      0    reboot dimension
   EREF     =   0.00    reference energy to select bands
   IMIX   =      4    mixing-type and parameters
     AMIX     =   0.40;   BMIX     =  1.00
     AMIX_MAG =   1.60;   BMIX_MAG =  1.00
     AMIN     =   0.10
     WC   =   100.;   INIMIX=   1;  MIXPRE=   1;  MAXMIX= -45

 Intra band minimization:
   WEIMIN = 0.0000     energy-eigenvalue tresh-hold
   EBREAK =  0.13E-07  absolut break condition
   DEPER  =   0.30     relativ break condition  

   TIME   =   0.40     timestep for ELM

  volume/ion in A,a.u.               =      15.67       105.77
  Fermi-wavevector in a.u.,A,eV,Ry     =   1.329772  2.512906 24.059112  1.768295
  Thomas-Fermi vector in A             =   2.458910
 
 Write flags
   LWAVE        =      F    write WAVECAR
   LDOWNSAMPLE  =      F    k-point downsampling of WAVECAR
   LCHARG       =      T    write CHGCAR
   LVTOT        =      F    write LOCPOT, total local potential
   LVHAR        =      T    write LOCPOT, Hartree potential only
   LELF         =      F    write electronic localiz. function (ELF)
   LORBIT       =     11    0 simple, 1 ext, 2 COOP (PROOUT), +10 PAW based schemes


 Dipole corrections
   LMONO  =      F    monopole corrections only (constant potential shift)
   LDIPOL =      F    correct potential (dipole corrections)
   IDIPOL =      0    1-x, 2-y, 3-z, 4-all directions 
   EPSILON=  1.0000000 bulk dielectric constant

 Exchange correlation treatment:
   GGA     =    --    GGA type
   LEXCH   =     8    internal setting for exchange type
   VOSKOWN=      0    Vosko Wilk Nusair interpolation
   LHFCALC =     F    Hartree Fock is set to
   LHFONE  =     F    Hartree Fock one center treatment
   AEXX    =    0.0000 exact exchange contribution

 Linear response parameters
   LEPSILON=     F    determine dielectric tensor
   LRPA    =     F    only Hartree local field effects (RPA)
   LNABLA  =     F    use nabla operator in PAW spheres
   LVEL    =     F    velocity operator in full k-point grid
   LINTERFAST=   F  fast interpolation
   KINTER  =     0    interpolate to denser k-point grid
   CSHIFT  =0.1000    complex shift for real part using Kramers Kronig
   OMEGAMAX=  -1.0    maximum frequency
   DEG_THRESHOLD= 0.2000000E-02 threshold for treating states as degnerate
   RTIME   =   -0.100 relaxation time in fs
  (WPLASMAI=    0.000 imaginary part of plasma frequency in eV, 0.658/RTIME)
   DFIELD  = 0.0000000 0.0000000 0.0000000 field for delta impulse in time
 
 Orbital magnetization related:
   ORBITALMAG=     F  switch on orbital magnetization
   LCHIMAG   =     F  perturbation theory with respect to B field
   DQ        =  0.001000  dq finite difference perturbation B field
   LLRAUG    =     F  two centre corrections for induced B field

 PEAD related settings:
   LPEAD      =     T    switch on PEAD
   IPEAD      =     4    finite difference order for dpsi/dk
   LCALCPOL   =     T    calculate macroscopic polarization
   LCALCEPS   =     F    calculate dielectric tensor
   EFIELD_PEAD=    0.0000    0.0000    0.0000
   SKIP_EDOTP =     F
   LRPA       =     F
   SKIP_SCF   =     F


--------------------------------------------------------------------------------------------------------


 Static calculation
 charge density and potential will be updated during run
 spin polarized calculation
 Variant of blocked Davidson
 Davidson routine will perform the subspace rotation
 perform sub-space diagonalisation
    after iterative eigenvector-optimisation
 modified Broyden-mixing scheme, WC =      100.0
 initial mixing is a Kerker type mixing with AMIX =  0.4000 and BMIX =      1.0000
 Hartree-type preconditioning will be used
 using additional bands           66
 real space projection scheme for non local part
 use partial core corrections
 calculate Harris-corrections to forces 
   (improved forces if not selfconsistent)
 use gradient corrections 
 use of overlap-Matrix (Vanderbilt PP)
 Fermi weights with tetrahedron method with Bloechl corrections


--------------------------------------------------------------------------------------------------------


  energy-cutoff  :      520.00
  volume of cell :      470.21
      direct lattice vectors                 reciprocal lattice vectors
     6.506500000  0.000000000  0.000000000     0.153692461  0.088734389  0.000000000
    -3.253250000  5.634794000  0.000000000     0.000000000  0.177468777  0.000000000
     0.000000000  0.000000000 12.825300000     0.000000000  0.000000000  0.077970886

  length of vectors
     6.506500000  6.506499749 12.825300000     0.177468770  0.177468777  0.077970886


 
 k-points in units of 2pi/SCALE and weight: pymatgen v1.1 with grid density = 1583 /
   0.00000000  0.00000000  0.00000000       0.031
   0.03842312  0.02218360  0.00000000       0.188
   0.07684623  0.04436719  0.00000000       0.094
   0.03842312  0.06655079  0.00000000       0.188
   0.00000000  0.00000000  0.03898544       0.031
   0.03842312  0.02218360  0.03898544       0.188
   0.07684623  0.04436719  0.03898544       0.094
   0.03842312  0.06655079  0.03898544       0.188
 
 k-points in reciprocal lattice and weights: pymatgen v1.1 with grid density = 1583 /
   0.00000000  0.00000000  0.00000000       0.031
   0.25000000  0.00000000  0.00000000       0.188
   0.50000000  0.00000000  0.00000000       0.094
   0.25000000  0.25000000  0.00000000       0.188
   0.00000000  0.00000000  0.50000000       0.031
   0.25000000  0.00000000  0.50000000       0.188
   0.50000000  0.00000000  0.50000000       0.094
   0.25000000  0.25000000  0.50000000       0.188
 
 position of ions in fractional coordinates (direct lattice) 
   0.00000000  0.00000000  0.26865800
   0.00000000  0.00000000  0.76865800
   0.33333300  0.66666700  0.23843600
   0.66666700  0.33333300  0.73843600
   0.33333300  0.66666700  0.73843600
   0.66666700  0.33333300  0.23843600
   0.00000000  0.00000000  0.48088900
   0.00000000  0.00000000  0.98088900
   0.33333300  0.66666700  0.01818700
   0.66666700  0.33333300  0.51818700
   0.33333300  0.66666700  0.51818700
   0.66666700  0.33333300  0.01818700
   0.31260100  0.00000000  0.16616400
   0.00000000  0.31260100  0.16616400
   0.68739900  0.68739900  0.16616400
   0.68739900  0.00000000  0.66616400
   0.00000000  0.68739900  0.66616400
   0.31260100  0.31260100  0.66616400
   0.64635000  0.00000000  0.33600400
   0.00000000  0.64635000  0.33600400
   0.35365000  0.35365000  0.33600400
   0.35365000  0.00000000  0.83600400
   0.00000000  0.35365000  0.83600400
   0.64635000  0.64635000  0.83600400
   0.33333300  0.00000000  0.00064900
   0.00000000  0.33333300  0.00064900
   0.66666700  0.66666700  0.00064900
   0.66666700  0.00000000  0.50064900
   0.00000000  0.66666700  0.50064900
   0.33333300  0.33333300  0.50064900
 
 position of ions in cartesian coordinates  (Angst):
   0.00000000  0.00000000  3.44561945
   0.00000000  0.00000000  9.85826945
  -0.00000325  3.75653121  3.05801323
   3.25325325  1.87826279  9.47066323
  -0.00000325  3.75653121  9.47066323
   3.25325325  1.87826279  3.05801323
   0.00000000  0.00000000  6.16754569
   0.00000000  0.00000000 12.58019569
  -0.00000325  3.75653121  0.23325373
   3.25325325  1.87826279  6.64590373
  -0.00000325  3.75653121  6.64590373
   3.25325325  1.87826279  0.23325373
   2.03393841  0.00000000  2.13110315
  -1.01696920  1.76144224  2.13110315
   2.23628080  3.87335176  2.13110315
   4.47256159  0.00000000  8.54375315
  -2.23628080  3.87335176  8.54375315
   1.01696920  1.76144224  8.54375315
   4.20547627  0.00000000  4.30935210
  -2.10273814  3.64204910  4.30935210
   1.15051186  1.99274490  4.30935210
   2.30102373  0.00000000 10.72200210
  -1.15051186  1.99274490 10.72200210
   2.10273814  3.64204910 10.72200210
   2.16883116  0.00000000  0.00832362
  -1.08441558  1.87826279  0.00832362
   2.16883442  3.75653121  0.00832362
   4.33766884  0.00000000  6.42097362
  -2.16883442  3.75653121  6.42097362
   1.08441558  1.87826279  6.42097362
 


--------------------------------------------------------------------------------------------------------


 k-point  1 :   0.0000 0.0000 0.0000  plane waves:   12587
 k-point  2 :   0.2500 0.0000 0.0000  plane waves:   12663
 k-point  3 :   0.5000 0.0000 0.0000  plane waves:   12644
 k-point  4 :   0.2500 0.2500 0.0000  plane waves:   12648
 k-point  5 :   0.0000 0.0000 0.5000  plane waves:   12660
 k-point  6 :   0.2500 0.0000 0.5000  plane waves:   12656
 k-point  7 :   0.5000 0.0000 0.5000  plane waves:   12652
 k-point  8 :   0.2500 0.2500 0.5000  plane waves:   12674

 maximum and minimum number of plane-waves per node :     12674    12587

 maximum number of plane-waves:     12674
 maximum index in each direction: 
   IXMAX=   12   IYMAX=   12   IZMAX=   23
   IXMIN=  -12   IYMIN=  -12   IZMIN=  -24


 serial   3D FFT for wavefunctions
 parallel 3D FFT for charge:
    minimum data exchange during FFTs selected (reduces bandwidth)


 total amount of memory used by VASP MPI-rank0    93604. kBytes
=======================================================================

   base      :      30000. kBytes
   nonlr-proj:      27111. kBytes
   fftplans  :       2081. kBytes
   grid      :      21338. kBytes
   one-center:       2949. kBytes
   wavefun   :      10125. kBytes
 
     INWAV:  cpu time    0.0001: real time    0.0001
 Broyden mixing: mesh for mixing (old mesh)
   NGX = 25   NGY = 25   NGZ = 47
  (NGX  = 96   NGY  = 96   NGZ  =192)
  gives a total of  29375 points

 initial charge density was supplied:
 charge density of overlapping atoms calculated
 number of electron     252.0000000 magnetization      18.0000000
 keeping initial charge density in first step


--------------------------------------------------------------------------------------------------------


 Maximum index for non-local projection operator         8244
 Maximum index for augmentation-charges          914 (set IRDMAX)


--------------------------------------------------------------------------------------------------------


 First call to EWALD:  gamma=   0.228
 Maximum number of real-space cells 4x 4x 2
 Maximum number of reciprocal cells 2x 2x 4

    FEWALD:  cpu time    0.0074: real time    0.0078


--------------------------------------- Iteration      1(   1)  ---------------------------------------


    POTLOK:  cpu time    0.2424: real time    0.2653
    SETDIJ:  cpu time    0.0603: real time    0.0615
     EDDAV:  cpu time   17.0256: real time   17.0249
 BZINTS: Fermi energy: 16.105343;252.000000 electrons
         Band energy: 917.363771;  BLOECHL correction: -0.091658
       DOS:  cpu time    0.0985: real time    0.0986
    --------------------------------------------
      LOOP:  cpu time   17.4269: real time   17.4504

 eigenvalue-minimisations  :  6976
 total energy-change (2. order) : 0.2260840E+04  (-0.1152797E+05)
 number of electron     252.0000000 magnetization      18.0000000
 augmentation part      252.0000000 magnetization      18.0000000

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5415.87925739
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       981.13668492
  PAW double counting   =     14689.62286600   -14939.65885362
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =       917.36377142
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      2260.84035560 eV

  energy without entropy =     2260.84035560  energy(sigma->0) =     2260.84035560


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   2)  ---------------------------------------


     EDDAV:  cpu time   14.7536: real time   14.7707
 BZINTS: Fermi energy:  9.828508;252.000000 electrons
         Band energy:***********;  BLOECHL correction: -0.039971
       DOS:  cpu time    0.0754: real time    0.0753
    --------------------------------------------
      LOOP:  cpu time   14.8290: real time   14.8460

 eigenvalue-minimisations  :  6144
 total energy-change (2. order) :-0.2079033E+04  (-0.1717551E+04)
 number of electron     252.0000000 magnetization      18.0000000
 augmentation part      252.0000000 magnetization      18.0000000

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5415.87925739
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       981.13668492
  PAW double counting   =     14689.62286600   -14939.65885362
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1161.66920347
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =       181.80738071 eV

  energy without entropy =      181.80738071  energy(sigma->0) =      181.80738071


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   3)  ---------------------------------------


     EDDAV:  cpu time   27.9755: real time   27.9896
 BZINTS: Fermi energy:  4.901786;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0246: real time    0.0245
    --------------------------------------------
      LOOP:  cpu time   28.0000: real time   28.0140

 eigenvalue-minimisations  :  9920
 total energy-change (2. order) :-0.4155889E+03  (-0.3956262E+03)
 number of electron     252.0000000 magnetization      18.0000000
 augmentation part      252.0000000 magnetization      18.0000000

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5415.87925739
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       981.13668492
  PAW double counting   =     14689.62286600   -14939.65885362
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1577.25810284
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -233.78151866 eV

  energy without entropy =     -233.78151866  energy(sigma->0) =     -233.78151866


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   4)  ---------------------------------------


     EDDAV:  cpu time   25.3850: real time   25.3992
 BZINTS: Fermi energy:  4.642652;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0243: real time    0.0243
    --------------------------------------------
      LOOP:  cpu time   25.4093: real time   25.4235

 eigenvalue-minimisations  :  9280
 total energy-change (2. order) :-0.1906418E+02  (-0.1904032E+02)
 number of electron     252.0000000 magnetization      18.0000000
 augmentation part      252.0000000 magnetization      18.0000000

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5415.87925739
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       981.13668492
  PAW double counting   =     14689.62286600   -14939.65885362
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1596.32228375
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -252.84569957 eV

  energy without entropy =     -252.84569957  energy(sigma->0) =     -252.84569957


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   5)  ---------------------------------------


     EDDAV:  cpu time   27.4413: real time   27.4574
 BZINTS: Fermi energy:  4.633160;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0224: real time    0.0224
    CHARGE:  cpu time    0.5842: real time    0.5881
    MIXING:  cpu time    0.0106: real time    0.0133
    --------------------------------------------
      LOOP:  cpu time   28.0585: real time   28.0812

 eigenvalue-minimisations  :  9792
 total energy-change (2. order) :-0.8192214E+00  (-0.8188358E+00)
 number of electron     252.0000092 magnetization      12.4256675
 augmentation part       46.9850255 magnetization      12.7879402

 Broyden mixing:
  rms(total) = 0.36954E+01    rms(broyden)= 0.36946E+01
  rms(prec ) = 0.45079E+01
  weight for this iteration     100.00

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5415.87925739
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       981.13668492
  PAW double counting   =     14689.62286600   -14939.65885362
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1597.14150519
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -253.66492101 eV

  energy without entropy =     -253.66492101  energy(sigma->0) =     -253.66492101


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   6)  ---------------------------------------


    POTLOK:  cpu time    0.2142: real time    0.2401
    SETDIJ:  cpu time    0.0553: real time    0.0553
     EDDAV:  cpu time   26.1148: real time   26.0980
 BZINTS: Fermi energy:  6.420677;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0225: real time    0.0225
    CHARGE:  cpu time    0.5805: real time    0.5800
    MIXING:  cpu time    0.0068: real time    0.0073
    --------------------------------------------
      LOOP:  cpu time   26.9942: real time   27.0034

 eigenvalue-minimisations  :  9792
 total energy-change (2. order) : 0.3483920E+02  (-0.2037754E+02)
 number of electron     252.0000066 magnetization       9.2101349
 augmentation part       43.7724630 magnetization       8.8952650

 Broyden mixing:
  rms(total) = 0.13961E+01    rms(broyden)= 0.13948E+01
  rms(prec ) = 0.15878E+01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   0.8252
  0.8252

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5660.57516159
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =      1001.10927125
  PAW double counting   =     15283.47242916   -15521.31053675
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1349.77686492
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -218.82571858 eV

  energy without entropy =     -218.82571858  energy(sigma->0) =     -218.82571858


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   7)  ---------------------------------------


    POTLOK:  cpu time    0.2141: real time    0.2408
    SETDIJ:  cpu time    0.0559: real time    0.0559
     EDDAV:  cpu time   25.6109: real time   25.5943
 BZINTS: Fermi energy:  6.098004;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0227: real time    0.0227
    CHARGE:  cpu time    0.5814: real time    0.5809
    MIXING:  cpu time    0.0040: real time    0.0040
    --------------------------------------------
      LOOP:  cpu time   26.4889: real time   26.4985

 eigenvalue-minimisations  :  9536
 total energy-change (2. order) : 0.2648859E+00  (-0.1327770E+01)
 number of electron     252.0000071 magnetization       3.8684764
 augmentation part       43.4460786 magnetization       3.5833340

 Broyden mixing:
  rms(total) = 0.77076E+00    rms(broyden)= 0.77073E+00
  rms(prec ) = 0.86933E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.1007
  0.9463  1.2552

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5629.05392533
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       998.56057940
  PAW double counting   =     15451.10298641   -15694.14983777
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1373.27577969
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -218.56083271 eV

  energy without entropy =     -218.56083271  energy(sigma->0) =     -218.56083271


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   8)  ---------------------------------------


    POTLOK:  cpu time    0.2359: real time    0.2416
    SETDIJ:  cpu time    0.0557: real time    0.0557
     EDDAV:  cpu time   33.8296: real time   33.8070
 BZINTS: Fermi energy:  5.967880;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0225: real time    0.0225
    CHARGE:  cpu time    0.5770: real time    0.5767
    MIXING:  cpu time    0.0044: real time    0.0044
    --------------------------------------------
      LOOP:  cpu time   34.7252: real time   34.7078

 eigenvalue-minimisations  : 11648
 total energy-change (2. order) :-0.6766266E+00  (-0.4027102E+00)
 number of electron     252.0000072 magnetization       1.9263380
 augmentation part       43.4853786 magnetization       1.7145081

 Broyden mixing:
  rms(total) = 0.30418E+00    rms(broyden)= 0.30405E+00
  rms(prec ) = 0.34346E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.2834
  2.1834  0.9165  0.7503

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5636.57329752
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       998.72579323
  PAW double counting   =     15748.15040270   -16001.21134213
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1356.58415986
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.23745932 eV

  energy without entropy =     -219.23745932  energy(sigma->0) =     -219.23745932


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(   9)  ---------------------------------------


    POTLOK:  cpu time    0.2118: real time    0.2394
    SETDIJ:  cpu time    0.0556: real time    0.0555
     EDDAV:  cpu time   20.2563: real time   20.2449
 BZINTS: Fermi energy:  6.102638;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0226: real time    0.0226
    CHARGE:  cpu time    0.5771: real time    0.5767
    MIXING:  cpu time    0.0044: real time    0.0044
    --------------------------------------------
      LOOP:  cpu time   21.1278: real time   21.1436

 eigenvalue-minimisations  :  8064
 total energy-change (2. order) :-0.2013029E+00  (-0.7221805E-01)
 number of electron     252.0000072 magnetization       0.8918462
 augmentation part       43.4094370 magnetization       0.7833926

 Broyden mixing:
  rms(total) = 0.12585E+00    rms(broyden)= 0.12584E+00
  rms(prec ) = 0.13848E+00
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.3465
  2.5244  1.2764  0.8896  0.6954

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5644.42003909
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.31811652
  PAW double counting   =     15829.30053127   -16085.34365360
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1346.54886159
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.43876222 eV

  energy without entropy =     -219.43876222  energy(sigma->0) =     -219.43876222


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  10)  ---------------------------------------


    POTLOK:  cpu time    0.2444: real time    0.2444
    SETDIJ:  cpu time    0.0551: real time    0.0551
     EDDAV:  cpu time   26.0584: real time   26.0414
 BZINTS: Fermi energy:  6.174461;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0250: real time    0.0250
    CHARGE:  cpu time    0.5757: real time    0.5752
    MIXING:  cpu time    0.0046: real time    0.0046
    --------------------------------------------
      LOOP:  cpu time   26.9632: real time   26.9457

 eigenvalue-minimisations  :  9536
 total energy-change (2. order) :-0.8085802E-01  (-0.2095675E-01)
 number of electron     252.0000071 magnetization       0.5659987
 augmentation part       43.3759251 magnetization       0.5220272

 Broyden mixing:
  rms(total) = 0.70861E-01    rms(broyden)= 0.70828E-01
  rms(prec ) = 0.79399E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.2603
  2.6320  1.4012  0.7957  0.7957  0.6769

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5648.41240798
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.51975044
  PAW double counting   =     15844.82945841   -16102.93460333
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1340.77696206
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.51962024 eV

  energy without entropy =     -219.51962024  energy(sigma->0) =     -219.51962024


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  11)  ---------------------------------------


    POTLOK:  cpu time    0.2114: real time    0.2391
    SETDIJ:  cpu time    0.0556: real time    0.0555
     EDDAV:  cpu time   27.1752: real time   27.1586
 BZINTS: Fermi energy:  6.144785;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0225: real time    0.0225
    CHARGE:  cpu time    0.5780: real time    0.5773
    MIXING:  cpu time    0.0048: real time    0.0048
    --------------------------------------------
      LOOP:  cpu time   28.0474: real time   28.0578

 eigenvalue-minimisations  :  9984
 total energy-change (2. order) :-0.7538343E-02  (-0.3286447E-02)
 number of electron     252.0000072 magnetization       0.1476091
 augmentation part       43.3723567 magnetization       0.1249153

 Broyden mixing:
  rms(total) = 0.35447E-01    rms(broyden)= 0.35444E-01
  rms(prec ) = 0.37650E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.3386
  2.8549  1.7875  1.1062  0.9241  0.7425  0.6165

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5648.20336035
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.46330286
  PAW double counting   =     15846.63452842   -16105.32290234
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1340.35387144
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.52715859 eV

  energy without entropy =     -219.52715859  energy(sigma->0) =     -219.52715859


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  12)  ---------------------------------------


    POTLOK:  cpu time    0.2119: real time    0.9534
    SETDIJ:  cpu time    0.0552: real time    0.0552
     EDDAV:  cpu time   25.5488: real time   25.5327
 BZINTS: Fermi energy:  6.122868;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0228: real time    0.0228
    CHARGE:  cpu time    0.5794: real time    0.5789
    MIXING:  cpu time    0.0053: real time    0.0053
    --------------------------------------------
      LOOP:  cpu time   26.4233: real time   27.1483

 eigenvalue-minimisations  :  9408
 total energy-change (2. order) :-0.8440246E-02  (-0.1770773E-02)
 number of electron     252.0000072 magnetization       0.0757784
 augmentation part       43.3875066 magnetization       0.0670223

 Broyden mixing:
  rms(total) = 0.19301E-01    rms(broyden)= 0.19294E-01
  rms(prec ) = 0.22315E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.2884
  2.9436  2.0087  1.1668  0.8152  0.7410  0.7410  0.6025

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5647.95710969
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.40760561
  PAW double counting   =     15844.48494531   -16103.58264848
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1340.14353585
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.53559883 eV

  energy without entropy =     -219.53559883  energy(sigma->0) =     -219.53559883


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  13)  ---------------------------------------


    POTLOK:  cpu time    0.2120: real time    0.2422
    SETDIJ:  cpu time    0.0557: real time    0.0557
     EDDAV:  cpu time   26.3148: real time   26.2989
 BZINTS: Fermi energy:  6.140910;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0248: real time    0.0248
    CHARGE:  cpu time    0.5790: real time    0.5786
    MIXING:  cpu time    0.0053: real time    0.0052
    --------------------------------------------
      LOOP:  cpu time   27.1916: real time   27.2055

 eigenvalue-minimisations  :  9792
 total energy-change (2. order) :-0.8919530E-04  (-0.2532614E-03)
 number of electron     252.0000072 magnetization       0.0262319
 augmentation part       43.3887126 magnetization       0.0226821

 Broyden mixing:
  rms(total) = 0.10437E-01    rms(broyden)= 0.10436E-01
  rms(prec ) = 0.11557E-01
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.3189
  3.0072  2.3987  1.0888  1.0888  0.9676  0.7512  0.6684  0.5806

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5648.41605518
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.43496463
  PAW double counting   =     15841.43841342   -16100.42455190
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1339.82360327
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.53568803 eV

  energy without entropy =     -219.53568803  energy(sigma->0) =     -219.53568803


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  14)  ---------------------------------------


    POTLOK:  cpu time    0.2122: real time    0.2414
    SETDIJ:  cpu time    0.0557: real time    0.0556
     EDDAV:  cpu time   25.9184: real time   25.9072
 BZINTS: Fermi energy:  6.140440;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0209: real time    0.0209
    CHARGE:  cpu time    0.5738: real time    0.5733
    MIXING:  cpu time    0.0051: real time    0.0051
    --------------------------------------------
      LOOP:  cpu time   26.7861: real time   26.8036

 eigenvalue-minimisations  :  9728
 total energy-change (2. order) :-0.4468014E-04  (-0.1267342E-03)
 number of electron     252.0000072 magnetization       0.0084596
 augmentation part       43.3847977 magnetization       0.0076206

 Broyden mixing:
  rms(total) = 0.29476E-02    rms(broyden)= 0.29432E-02
  rms(prec ) = 0.34075E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.3007
  3.1223  2.3613  1.3782  1.0937  0.8488  0.8488  0.8272  0.6570  0.5687

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5648.95602286
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.46648534
  PAW double counting   =     15838.66861018   -16097.53118039
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1339.43876924
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.53573271 eV

  energy without entropy =     -219.53573271  energy(sigma->0) =     -219.53573271


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  15)  ---------------------------------------


    POTLOK:  cpu time    0.2112: real time    0.2404
    SETDIJ:  cpu time    0.0553: real time    0.0552
     EDDAV:  cpu time   22.3598: real time   22.3701
 BZINTS: Fermi energy:  6.140586;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0226: real time    0.0226
    CHARGE:  cpu time    0.5783: real time    0.5779
    MIXING:  cpu time    0.0057: real time    0.0057
    --------------------------------------------
      LOOP:  cpu time   23.2328: real time   23.2718

 eigenvalue-minimisations  :  9024
 total energy-change (2. order) :-0.1587059E-04  (-0.1010815E-04)
 number of electron     252.0000072 magnetization       0.0039171
 augmentation part       43.3846457 magnetization       0.0038011

 Broyden mixing:
  rms(total) = 0.14263E-02    rms(broyden)= 0.14258E-02
  rms(prec ) = 0.16017E-02
  weight for this iteration     100.00

 eigenvalues of (default mixing * dielectric matrix)
  average eigenvalue GAMMA=   1.2826
  3.1631  2.2905  1.7385  1.1050  0.8905  0.8905  0.7787  0.7787  0.6321  0.5585

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5649.04927181
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.46875355
  PAW double counting   =     15837.86152341   -16096.68343074
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1339.38846726
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.53574858 eV

  energy without entropy =     -219.53574858  energy(sigma->0) =     -219.53574858


--------------------------------------------------------------------------------------------------------




--------------------------------------- Iteration      1(  16)  ---------------------------------------


    POTLOK:  cpu time    0.2131: real time    0.2404
    SETDIJ:  cpu time    0.0554: real time    0.0553
     EDDAV:  cpu time   20.4884: real time   20.4771
 BZINTS: Fermi energy:  6.140520;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000
       DOS:  cpu time    0.0220: real time    0.0220
    --------------------------------------------
      LOOP:  cpu time   20.7789: real time   20.7949

 eigenvalue-minimisations  :  7616
 total energy-change (2. order) :-0.1392287E-05  (-0.1563490E-05)
 number of electron     252.0000072 magnetization       0.0039171
 augmentation part       43.3846457 magnetization       0.0038011

 Free energy of the ion-electron system (eV)
  ---------------------------------------------------
  alpha Z        PSCENC =       865.42118643
  Ewald energy   TEWEN  =    -18162.44880550
  -Hartree energ DENC   =     -5649.06024777
  -exchange      EXHF   =         0.00000000
  -V(xc)+E(xc)   XCENC  =       999.46772438
  PAW double counting   =     15837.57004162   -16096.37502326
  entropy T*S    EENTRO =         0.00000000
  eigenvalues    EBANDS =     -1339.39338920
  atomic energy  EATOM  =     23325.28276333
  Solvation  Ediel_sol  =         0.00000000
  ---------------------------------------------------
  free energy    TOTEN  =      -219.53574997 eV

  energy without entropy =     -219.53574997  energy(sigma->0) =     -219.53574997


--------------------------------------------------------------------------------------------------------




 average (electrostatic) potential at core
  the test charge radii are     1.1480  0.7089  0.7215  1.0894
  (the norm of the test charge is              1.0000)
       1 -95.1033       2 -95.1033       3 -95.1487       4 -95.1487       5 -95.1487
       6 -95.1487       7 -60.8822       8 -60.8823       9 -60.8664      10 -60.8664
      11 -60.8664      12 -60.8664      13 -67.6299      14 -67.6299      15 -67.6297
      16 -67.6300      17 -67.6300      18 -67.6297      19 -67.6119      20 -67.6119
      21 -67.6117      22 -67.6119      23 -67.6119      24 -67.6117      25 -99.9948
      26 -99.9948      27 -99.9946      28 -99.9948      29 -99.9948      30 -99.9947
 
 
 
 E-fermi :   6.1405     XC(G=0): -10.4637     alpha+bet :-14.8976


 spin component 1

 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -23.8512      1.00000
      2     -23.8505      1.00000
      3     -23.6292      1.00000
      4     -23.6292      1.00000
      5     -23.5997      1.00000
      6     -23.5997      1.00000
      7     -14.5957      1.00000
      8     -14.5760      1.00000
      9     -14.5302      1.00000
     10     -14.5301      1.00000
     11     -14.5298      1.00000
     12     -14.5297      1.00000
     13     -14.4194      1.00000
     14     -14.4190      1.00000
     15     -14.4145      1.00000
     16     -14.4140      1.00000
     17     -14.3737      1.00000
     18     -14.3737      1.00000
     19     -14.3737      1.00000
     20     -14.3737      1.00000
     21     -14.3071      1.00000
     22     -14.3071      1.00000
     23     -14.3070      1.00000
     24     -14.3070      1.00000
     25     -14.2464      1.00000
     26     -14.2450      1.00000
     27     -14.2432      1.00000
     28     -14.2415      1.00000
     29     -14.2377      1.00000
     30     -14.2376      1.00000
     31     -14.2376      1.00000
     32     -14.2376      1.00000
     33     -14.2211      1.00000
     34     -14.2211      1.00000
     35     -14.2210      1.00000
     36     -14.2210      1.00000
     37     -11.3955      1.00000
     38     -11.3954      1.00000
     39     -11.3812      1.00000
     40     -11.3810      1.00000
     41     -11.2561      1.00000
     42     -11.2208      1.00000
     43     -11.1303      1.00000
     44     -11.1302      1.00000
     45     -11.1156      1.00000
     46     -11.1154      1.00000
     47     -10.4841      1.00000
     48     -10.4433      1.00000
     49      -7.9195      1.00000
     50      -7.8405      1.00000
     51      -7.8124      1.00000
     52      -7.7146      1.00000
     53      -7.2881      1.00000
     54      -7.2881      1.00000
     55      -7.2867      1.00000
     56      -7.2867      1.00000
     57      -7.0363      1.00000
     58      -7.0362      1.00000
     59      -7.0316      1.00000
     60      -7.0315      1.00000
     61      -6.9079      1.00000
     62      -6.8973      1.00000
     63      -6.8972      1.00000
     64      -6.8924      1.00000
     65      -6.8924      1.00000
     66      -6.7583      1.00000
     67      -6.5890      1.00000
     68      -6.4637      1.00000
     69      -5.9546      1.00000
     70      -5.8301      1.00000
     71      -5.7919      1.00000
     72      -5.7708      1.00000
     73      -0.3940      1.00000
     74      -0.3939      1.00000
     75      -0.3788      1.00000
     76      -0.3787      1.00000
     77       1.6666      1.00000
     78       1.8571      1.00000
     79       1.8571      1.00000
     80       1.8584      1.00000
     81       1.8585      1.00000
     82       1.9791      1.00000
     83       2.4264      1.00000
     84       2.5322      1.00000
     85       2.7460      1.00000
     86       2.7460      1.00000
     87       2.9728      1.00000
     88       2.9730      1.00000
     89       3.1820      1.00000
     90       3.1821      1.00000
     91       3.1844      1.00000
     92       3.1844      1.00000
     93       3.6463      1.00000
     94       3.6659      1.00000
     95       3.9251      1.00000
     96       3.9436      1.00000
     97       4.0914      1.00000
     98       4.0914      1.00000
     99       4.1461      1.00000
    100       4.1461      1.00000
    101       4.3700      1.00000
    102       4.3700      1.00000
    103       4.6398      1.00000
    104       4.6398      1.00000
    105       4.6638      1.00000
    106       4.6638      1.00000
    107       4.8017      1.00000
    108       4.8076      1.00000
    109       4.8083      1.00000
    110       4.8345      1.00000
    111       4.8528      1.00000
    112       4.8529      1.00000
    113       5.0278      1.00000
    114       5.0278      1.00000
    115       5.0750      1.00000
    116       5.0751      1.00000
    117       5.0906      1.00000
    118       5.0907      1.00000
    119       5.3531      1.00000
    120       5.3531      1.00000
    121       5.6056      1.00000
    122       5.6547      1.00000
    123       5.8476      1.00000
    124       5.9938      1.00000
    125       6.0305      1.00000
    126       6.0596      1.00000
    127       7.2548      0.00000
    128       7.8553      0.00000
    129      10.4235      0.00000
    130      10.4235      0.00000
    131      10.4441      0.00000
    132      10.4442      0.00000
    133      10.6100      0.00000
    134      10.6100      0.00000
    135      10.6710      0.00000
    136      10.6710      0.00000
    137      10.7007      0.00000
    138      10.7231      0.00000
    139      10.7898      0.00000
    140      10.7925      0.00000
    141      10.9225      0.00000
    142      10.9225      0.00000
    143      10.9366      0.00000
    144      10.9366      0.00000
    145      10.9496      0.00000
    146      10.9507      0.00000
    147      11.0312      0.00000
    148      11.0312      0.00000
    149      11.0776      0.00000
    150      11.0776      0.00000
    151      11.1274      0.00000
    152      11.1274      0.00000
    153      11.1303      0.00000
    154      11.1644      0.00000
    155      11.1644      0.00000
    156      11.1721      0.00000
    157      11.1798      0.00000
    158      11.1942      0.00000
    159      11.3081      0.00000
    160      11.3188      0.00000
    161      11.3581      0.00000
    162      11.3813      0.00000
    163      11.3878      0.00000
    164      11.4310      0.00000
    165      11.4310      0.00000
    166      11.4440      0.00000
    167      11.4440      0.00000
    168      11.4517      0.00000
    169      11.5398      0.00000
    170      11.5398      0.00000
    171      11.5642      0.00000
    172      11.5751      0.00000
    173      11.5751      0.00000
    174      11.6250      0.00000
    175      11.6250      0.00000
    176      11.6259      0.00000
    177      11.6817      0.00000
    178      11.6817      0.00000
    179      11.9870      0.00000
    180      12.0086      0.00000
    181      12.2054      0.00000
    182      12.2054      0.00000
    183      12.3511      0.00000
    184      12.3512      0.00000
    185      12.3785      0.00000
    186      12.6656      0.00000
    187      12.6656      0.00000
    188      12.6931      0.00000
    189      12.6931      0.00000
    190      12.7409      0.00000
    191      12.9028      0.00000
    192      12.9028      0.00000

 k-point     2 :       0.2500    0.0000    0.0000
  band No.  band energies     occupation 
      1     -23.8134      1.00000
      2     -23.8128      1.00000
      3     -23.6516      1.00000
      4     -23.6515      1.00000
      5     -23.6179      1.00000
      6     -23.6178      1.00000
      7     -14.5844      1.00000
      8     -14.5715      1.00000
      9     -14.5491      1.00000
     10     -14.5456      1.00000
     11     -14.5354      1.00000
     12     -14.5321      1.00000
     13     -14.4190      1.00000
     14     -14.4166      1.00000
     15     -14.4133      1.00000
     16     -14.4121      1.00000
     17     -14.3826      1.00000
     18     -14.3817      1.00000
     19     -14.3586      1.00000
     20     -14.3586      1.00000
     21     -14.3080      1.00000
     22     -14.3079      1.00000
     23     -14.3009      1.00000
     24     -14.3008      1.00000
     25     -14.2487      1.00000
     26     -14.2472      1.00000
     27     -14.2452      1.00000
     28     -14.2438      1.00000
     29     -14.2328      1.00000
     30     -14.2327      1.00000
     31     -14.2314      1.00000
     32     -14.2313      1.00000
     33     -14.2267      1.00000
     34     -14.2263      1.00000
     35     -14.2223      1.00000
     36     -14.2222      1.00000
     37     -11.6970      1.00000
     38     -11.5488      1.00000
     39     -11.4852      1.00000
     40     -11.3537      1.00000
     41     -11.1704      1.00000
     42     -11.1159      1.00000
     43     -11.0411      1.00000
     44     -10.8955      1.00000
     45     -10.8572      1.00000
     46     -10.7518      1.00000
     47     -10.6792      1.00000
     48     -10.6715      1.00000
     49      -7.9335      1.00000
     50      -7.9101      1.00000
     51      -7.8806      1.00000
     52      -7.7835      1.00000
     53      -7.4189      1.00000
     54      -7.4167      1.00000
     55      -7.3220      1.00000
     56      -7.3124      1.00000
     57      -7.2687      1.00000
     58      -7.1527      1.00000
     59      -7.1453      1.00000
     60      -7.0397      1.00000
     61      -6.9026      1.00000
     62      -6.8251      1.00000
     63      -6.7482      1.00000
     64      -6.7435      1.00000
     65      -6.6142      1.00000
     66      -6.5899      1.00000
     67      -6.5190      1.00000
     68      -6.3322      1.00000
     69      -6.0059      1.00000
     70      -5.9373      1.00000
     71      -5.8917      1.00000
     72      -5.8537      1.00000
     73      -0.4893      1.00000
     74      -0.3741      1.00000
     75      -0.2162      1.00000
     76      -0.2017      1.00000
     77       1.5214      1.00000
     78       1.7915      1.00000
     79       1.8944      1.00000
     80       2.1921      1.00000
     81       2.3180      1.00000
     82       2.3342      1.00000
     83       2.3607      1.00000
     84       2.4128      1.00000
     85       2.6323      1.00000
     86       2.7207      1.00000
     87       2.8450      1.00000
     88       2.9120      1.00000
     89       3.2077      1.00000
     90       3.2678      1.00000
     91       3.3475      1.00000
     92       3.4457      1.00000
     93       3.5247      1.00000
     94       3.6413      1.00000
     95       3.7398      1.00000
     96       3.8924      1.00000
     97       3.9759      1.00000
     98       4.0335      1.00000
     99       4.1502      1.00000
    100       4.2731      1.00000
    101       4.2944      1.00000
    102       4.3540      1.00000
    103       4.4035      1.00000
    104       4.4368      1.00000
    105       4.4448      1.00000
    106       4.5656      1.00000
    107       4.5987      1.00000
    108       4.6543      1.00000
    109       4.7188      1.00000
    110       4.7851      1.00000
    111       4.8047      1.00000
    112       4.8662      1.00000
    113       4.8840      1.00000
    114       4.9517      1.00000
    115       5.0720      1.00000
    116       5.1530      1.00000
    117       5.1608      1.00000
    118       5.2224      1.00000
    119       5.2576      1.00000
    120       5.2764      1.00000
    121       5.8206      1.00000
    122       5.9004      1.00000
    123       5.9531      1.00000
    124       5.9817      1.00000
    125       5.9983      1.00000
    126       6.0444      1.00000
    127       8.0495      0.00000
    128       8.3006      0.00000
    129       9.8404      0.00000
    130      10.0588      0.00000
    131      10.1066      0.00000
    132      10.3474      0.00000
    133      10.4942      0.00000
    134      10.5063      0.00000
    135      10.6200      0.00000
    136      10.6769      0.00000
    137      10.7373      0.00000
    138      10.7508      0.00000
    139      10.7822      0.00000
    140      10.8359      0.00000
    141      10.8562      0.00000
    142      10.8866      0.00000
    143      10.9125      0.00000
    144      10.9270      0.00000
    145      10.9363      0.00000
    146      10.9722      0.00000
    147      11.0402      0.00000
    148      11.0416      0.00000
    149      11.0467      0.00000
    150      11.0572      0.00000
    151      11.0993      0.00000
    152      11.1211      0.00000
    153      11.1251      0.00000
    154      11.1379      0.00000
    155      11.1618      0.00000
    156      11.1626      0.00000
    157      11.1734      0.00000
    158      11.1845      0.00000
    159      11.2057      0.00000
    160      11.2136      0.00000
    161      11.2273      0.00000
    162      11.2381      0.00000
    163      11.2566      0.00000
    164      11.2867      0.00000
    165      11.3107      0.00000
    166      11.3315      0.00000
    167      11.3565      0.00000
    168      11.3567      0.00000
    169      11.3893      0.00000
    170      11.4686      0.00000
    171      11.4704      0.00000
    172      11.5089      0.00000
    173      11.5196      0.00000
    174      11.5364      0.00000
    175      11.6989      0.00000
    176      11.7306      0.00000
    177      11.8037      0.00000
    178      11.8440      0.00000
    179      11.8836      0.00000
    180      11.9927      0.00000
    181      12.0184      0.00000
    182      12.0742      0.00000
    183      12.1491      0.00000
    184      12.3003      0.00000
    185      12.3123      0.00000
    186      12.3239      0.00000
    187      12.4454      0.00000
    188      12.5194      0.00000
    189      12.5317      0.00000
    190      12.6185      0.00000
    191      12.7598      0.00000
    192      12.9060      0.00000

 k-point     3 :       0.5000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -23.7383      1.00000
      2     -23.7379      1.00000
      3     -23.7120      1.00000
      4     -23.7117      1.00000
      5     -23.6353      1.00000
      6     -23.6351      1.00000
      7     -14.5629      1.00000
      8     -14.5619      1.00000
      9     -14.5604      1.00000
     10     -14.5595      1.00000
     11     -14.5589      1.00000
     12     -14.5399      1.00000
     13     -14.4207      1.00000
     14     -14.4146      1.00000
     15     -14.4008      1.00000
     16     -14.4007      1.00000
     17     -14.3993      1.00000
     18     -14.3989      1.00000
     19     -14.3280      1.00000
     20     -14.3279      1.00000
     21     -14.3265      1.00000
     22     -14.3264      1.00000
     23     -14.2937      1.00000
     24     -14.2936      1.00000
     25     -14.2507      1.00000
     26     -14.2491      1.00000
     27     -14.2389      1.00000
     28     -14.2379      1.00000
     29     -14.2378      1.00000
     30     -14.2368      1.00000
     31     -14.2277      1.00000
     32     -14.2274      1.00000
     33     -14.2269      1.00000
     34     -14.2269      1.00000
     35     -14.2252      1.00000
     36     -14.2251      1.00000
     37     -11.8251      1.00000
     38     -11.6143      1.00000
     39     -11.3192      1.00000
     40     -11.2691      1.00000
     41     -11.1917      1.00000
     42     -11.1038      1.00000
     43     -11.0179      1.00000
     44     -10.9543      1.00000
     45     -10.8836      1.00000
     46     -10.8674      1.00000
     47     -10.6147      1.00000
     48     -10.4002      1.00000
     49      -7.9593      1.00000
     50      -7.8944      1.00000
     51      -7.6911      1.00000
     52      -7.6899      1.00000
     53      -7.6620      1.00000
     54      -7.6449      1.00000
     55      -7.6061      1.00000
     56      -7.6010      1.00000
     57      -7.2276      1.00000
     58      -7.1119      1.00000
     59      -7.0271      1.00000
     60      -6.9821      1.00000
     61      -6.9600      1.00000
     62      -6.8637      1.00000
     63      -6.7983      1.00000
     64      -6.7967      1.00000
     65      -6.4008      1.00000
     66      -6.3521      1.00000
     67      -6.3398      1.00000
     68      -6.2831      1.00000
     69      -6.0921      1.00000
     70      -6.0593      1.00000
     71      -5.9470      1.00000
     72      -5.9064      1.00000
     73      -0.6088      1.00000
     74      -0.3943      1.00000
     75       0.5291      1.00000
     76       0.5452      1.00000
     77       0.6052      1.00000
     78       0.6235      1.00000
     79       2.2824      1.00000
     80       2.3259      1.00000
     81       2.4208      1.00000
     82       2.4234      1.00000
     83       2.4379      1.00000
     84       2.6558      1.00000
     85       2.8618      1.00000
     86       2.9735      1.00000
     87       3.2629      1.00000
     88       3.2751      1.00000
     89       3.3249      1.00000
     90       3.3519      1.00000
     91       3.3963      1.00000
     92       3.4708      1.00000
     93       3.5349      1.00000
     94       3.5587      1.00000
     95       3.7088      1.00000
     96       3.7351      1.00000
     97       3.7411      1.00000
     98       3.9743      1.00000
     99       3.9918      1.00000
    100       4.0563      1.00000
    101       4.0650      1.00000
    102       4.2116      1.00000
    103       4.3273      1.00000
    104       4.3483      1.00000
    105       4.4743      1.00000
    106       4.4802      1.00000
    107       4.4954      1.00000
    108       4.5529      1.00000
    109       4.6028      1.00000
    110       4.6621      1.00000
    111       4.7548      1.00000
    112       4.8158      1.00000
    113       4.8481      1.00000
    114       4.9322      1.00000
    115       5.0111      1.00000
    116       5.1638      1.00000
    117       5.1669      1.00000
    118       5.3073      1.00000
    119       5.3310      1.00000
    120       5.3532      1.00000
    121       5.8214      1.00000
    122       5.8869      1.00000
    123       5.9597      1.00000
    124       5.9714      1.00000
    125       6.0970      1.00000
    126       6.1331      1.00000
    127       9.2207      0.00000
    128       9.2620      0.00000
    129       9.3051      0.00000
    130       9.3233      0.00000
    131       9.3401      0.00000
    132       9.3643      0.00000
    133      10.2916      0.00000
    134      10.3342      0.00000
    135      10.4507      0.00000
    136      10.4722      0.00000
    137      10.7852      0.00000
    138      10.8121      0.00000
    139      10.8449      0.00000
    140      10.8649      0.00000
    141      10.8728      0.00000
    142      10.9293      0.00000
    143      10.9475      0.00000
    144      10.9492      0.00000
    145      10.9661      0.00000
    146      10.9746      0.00000
    147      11.0063      0.00000
    148      11.0232      0.00000
    149      11.0418      0.00000
    150      11.0669      0.00000
    151      11.0676      0.00000
    152      11.0850      0.00000
    153      11.1009      0.00000
    154      11.1323      0.00000
    155      11.1327      0.00000
    156      11.1443      0.00000
    157      11.1468      0.00000
    158      11.1508      0.00000
    159      11.1518      0.00000
    160      11.1618      0.00000
    161      11.1683      0.00000
    162      11.1698      0.00000
    163      11.1923      0.00000
    164      11.2491      0.00000
    165      11.2514      0.00000
    166      11.2882      0.00000
    167      11.2924      0.00000
    168      11.3212      0.00000
    169      11.4310      0.00000
    170      11.4448      0.00000
    171      11.4618      0.00000
    172      11.4857      0.00000
    173      11.4985      0.00000
    174      11.5007      0.00000
    175      11.5145      0.00000
    176      11.5830      0.00000
    177      11.6816      0.00000
    178      11.7220      0.00000
    179      11.8274      0.00000
    180      11.8317      0.00000
    181      11.9855      0.00000
    182      12.0270      0.00000
    183      12.2348      0.00000
    184      12.2395      0.00000
    185      12.2583      0.00000
    186      12.3151      0.00000
    187      12.5260      0.00000
    188      12.5744      0.00000
    189      12.6483      0.00000
    190      12.6553      0.00000
    191      12.9262      0.00000
    192      13.0033      0.00000

 k-point     4 :       0.2500    0.2500    0.0000
  band No.  band energies     occupation 
      1     -23.7529      1.00000
      2     -23.7525      1.00000
      3     -23.6841      1.00000
      4     -23.6839      1.00000
      5     -23.6486      1.00000
      6     -23.6484      1.00000
      7     -14.5695      1.00000
      8     -14.5635      1.00000
      9     -14.5606      1.00000
     10     -14.5604      1.00000
     11     -14.5447      1.00000
     12     -14.5447      1.00000
     13     -14.4141      1.00000
     14     -14.4119      1.00000
     15     -14.4108      1.00000
     16     -14.4086      1.00000
     17     -14.3969      1.00000
     18     -14.3948      1.00000
     19     -14.3339      1.00000
     20     -14.3338      1.00000
     21     -14.3069      1.00000
     22     -14.3068      1.00000
     23     -14.3061      1.00000
     24     -14.3060      1.00000
     25     -14.2469      1.00000
     26     -14.2455      1.00000
     27     -14.2454      1.00000
     28     -14.2440      1.00000
     29     -14.2360      1.00000
     30     -14.2351      1.00000
     31     -14.2280      1.00000
     32     -14.2278      1.00000
     33     -14.2271      1.00000
     34     -14.2269      1.00000
     35     -14.2239      1.00000
     36     -14.2237      1.00000
     37     -11.6523      1.00000
     38     -11.6316      1.00000
     39     -11.4609      1.00000
     40     -11.4350      1.00000
     41     -11.2722      1.00000
     42     -11.1383      1.00000
     43     -10.9188      1.00000
     44     -10.8092      1.00000
     45     -10.7924      1.00000
     46     -10.7714      1.00000
     47     -10.6082      1.00000
     48     -10.5964      1.00000
     49      -7.8380      1.00000
     50      -7.8162      1.00000
     51      -7.8095      1.00000
     52      -7.7371      1.00000
     53      -7.4735      1.00000
     54      -7.4721      1.00000
     55      -7.4585      1.00000
     56      -7.4362      1.00000
     57      -7.4106      1.00000
     58      -7.3952      1.00000
     59      -7.3667      1.00000
     60      -7.2975      1.00000
     61      -6.8744      1.00000
     62      -6.7798      1.00000
     63      -6.6424      1.00000
     64      -6.5397      1.00000
     65      -6.5126      1.00000
     66      -6.4645      1.00000
     67      -6.4327      1.00000
     68      -6.1050      1.00000
     69      -6.0887      1.00000
     70      -6.0542      1.00000
     71      -5.9746      1.00000
     72      -5.8974      1.00000
     73      -0.1966      1.00000
     74      -0.1804      1.00000
     75      -0.0366      1.00000
     76      -0.0347      1.00000
     77       0.7325      1.00000
     78       0.9104      1.00000
     79       2.1149      1.00000
     80       2.3385      1.00000
     81       2.3744      1.00000
     82       2.4355      1.00000
     83       2.6167      1.00000
     84       2.7526      1.00000
     85       2.9062      1.00000
     86       3.1174      1.00000
     87       3.1429      1.00000
     88       3.1570      1.00000
     89       3.2813      1.00000
     90       3.3656      1.00000
     91       3.4566      1.00000
     92       3.5507      1.00000
     93       3.5677      1.00000
     94       3.7152      1.00000
     95       3.7249      1.00000
     96       3.7426      1.00000
     97       3.7610      1.00000
     98       3.8259      1.00000
     99       3.8479      1.00000
    100       4.0106      1.00000
    101       4.0485      1.00000
    102       4.0687      1.00000
    103       4.2483      1.00000
    104       4.4149      1.00000
    105       4.4636      1.00000
    106       4.5010      1.00000
    107       4.5635      1.00000
    108       4.6711      1.00000
    109       4.6990      1.00000
    110       4.7958      1.00000
    111       4.8002      1.00000
    112       4.8408      1.00000
    113       4.8751      1.00000
    114       4.8830      1.00000
    115       5.1352      1.00000
    116       5.1551      1.00000
    117       5.1628      1.00000
    118       5.2655      1.00000
    119       5.2682      1.00000
    120       5.2737      1.00000
    121       5.8160      1.00000
    122       5.9078      1.00000
    123       5.9176      1.00000
    124       6.0180      1.00000
    125       6.0678      1.00000
    126       6.0872      1.00000
    127       8.8010      0.00000
    128       8.9093      0.00000
    129       9.3636      0.00000
    130       9.4110      0.00000
    131       9.4367      0.00000
    132       9.4973      0.00000
    133      10.3490      0.00000
    134      10.4441      0.00000
    135      10.5256      0.00000
    136      10.5414      0.00000
    137      10.5795      0.00000
    138      10.6222      0.00000
    139      10.8849      0.00000
    140      10.8881      0.00000
    141      10.8975      0.00000
    142      10.9073      0.00000
    143      10.9528      0.00000
    144      10.9652      0.00000
    145      10.9782      0.00000
    146      10.9970      0.00000
    147      11.0230      0.00000
    148      11.0423      0.00000
    149      11.0579      0.00000
    150      11.0859      0.00000
    151      11.0993      0.00000
    152      11.1101      0.00000
    153      11.1127      0.00000
    154      11.1195      0.00000
    155      11.1310      0.00000
    156      11.1412      0.00000
    157      11.1503      0.00000
    158      11.1602      0.00000
    159      11.1616      0.00000
    160      11.1693      0.00000
    161      11.1733      0.00000
    162      11.1822      0.00000
    163      11.2414      0.00000
    164      11.2597      0.00000
    165      11.2761      0.00000
    166      11.2820      0.00000
    167      11.3273      0.00000
    168      11.3768      0.00000
    169      11.4490      0.00000
    170      11.4560      0.00000
    171      11.4720      0.00000
    172      11.4970      0.00000
    173      11.5038      0.00000
    174      11.5130      0.00000
    175      11.5740      0.00000
    176      11.6419      0.00000
    177      11.6425      0.00000
    178      11.7017      0.00000
    179      11.7070      0.00000
    180      11.8188      0.00000
    181      11.9345      0.00000
    182      12.0109      0.00000
    183      12.0116      0.00000
    184      12.2327      0.00000
    185      12.3915      0.00000
    186      12.3982      0.00000
    187      12.5064      0.00000
    188      12.5272      0.00000
    189      12.7868      0.00000
    190      12.7901      0.00000
    191      12.8627      0.00000
    192      12.9785      0.00000

 k-point     5 :       0.0000    0.0000    0.5000
  band No.  band energies     occupation 
      1     -23.8509      1.00000
      2     -23.8509      1.00000
      3     -23.6292      1.00000
      4     -23.6292      1.00000
      5     -23.5997      1.00000
      6     -23.5997      1.00000
      7     -14.5859      1.00000
      8     -14.5858      1.00000
      9     -14.5300      1.00000
     10     -14.5300      1.00000
     11     -14.5300      1.00000
     12     -14.5299      1.00000
     13     -14.4193      1.00000
     14     -14.4192      1.00000
     15     -14.4143      1.00000
     16     -14.4142      1.00000
     17     -14.3737      1.00000
     18     -14.3737      1.00000
     19     -14.3737      1.00000
     20     -14.3737      1.00000
     21     -14.3071      1.00000
     22     -14.3071      1.00000
     23     -14.3070      1.00000
     24     -14.3070      1.00000
     25     -14.2457      1.00000
     26     -14.2457      1.00000
     27     -14.2424      1.00000
     28     -14.2423      1.00000
     29     -14.2376      1.00000
     30     -14.2376      1.00000
     31     -14.2376      1.00000
     32     -14.2376      1.00000
     33     -14.2211      1.00000
     34     -14.2211      1.00000
     35     -14.2211      1.00000
     36     -14.2210      1.00000
     37     -11.3888      1.00000
     38     -11.3888      1.00000
     39     -11.3887      1.00000
     40     -11.3887      1.00000
     41     -11.2392      1.00000
     42     -11.2392      1.00000
     43     -11.1225      1.00000
     44     -11.1225      1.00000
     45     -11.1224      1.00000
     46     -11.1224      1.00000
     47     -10.4626      1.00000
     48     -10.4626      1.00000
     49      -7.8860      1.00000
     50      -7.8860      1.00000
     51      -7.7653      1.00000
     52      -7.7653      1.00000
     53      -7.2875      1.00000
     54      -7.2875      1.00000
     55      -7.2875      1.00000
     56      -7.2875      1.00000
     57      -7.0344      1.00000
     58      -7.0344      1.00000
     59      -7.0343      1.00000
     60      -7.0343      1.00000
     61      -6.8944      1.00000
     62      -6.8944      1.00000
     63      -6.8944      1.00000
     64      -6.8944      1.00000
     65      -6.8168      1.00000
     66      -6.8168      1.00000
     67      -6.5423      1.00000
     68      -6.5423      1.00000
     69      -5.8717      1.00000
     70      -5.8717      1.00000
     71      -5.7954      1.00000
     72      -5.7954      1.00000
     73      -0.3849      1.00000
     74      -0.3849      1.00000
     75      -0.3848      1.00000
     76      -0.3848      1.00000
     77       1.8130      1.00000
     78       1.8130      1.00000
     79       1.8441      1.00000
     80       1.8441      1.00000
     81       1.8441      1.00000
     82       1.8441      1.00000
     83       2.4833      1.00000
     84       2.4833      1.00000
     85       2.8739      1.00000
     86       2.8740      1.00000
     87       2.8740      1.00000
     88       2.8740      1.00000
     89       3.1838      1.00000
     90       3.1838      1.00000
     91       3.1838      1.00000
     92       3.1838      1.00000
     93       3.7886      1.00000
     94       3.7886      1.00000
     95       3.7938      1.00000
     96       3.7938      1.00000
     97       4.0996      1.00000
     98       4.0996      1.00000
     99       4.0996      1.00000
    100       4.0996      1.00000
    101       4.5448      1.00000
    102       4.5448      1.00000
    103       4.5448      1.00000
    104       4.5448      1.00000
    105       4.6888      1.00000
    106       4.6888      1.00000
    107       4.6888      1.00000
    108       4.6888      1.00000
    109       4.8011      1.00000
    110       4.8011      1.00000
    111       4.8214      1.00000
    112       4.8214      1.00000
    113       5.0640      1.00000
    114       5.0640      1.00000
    115       5.0641      1.00000
    116       5.0641      1.00000
    117       5.2551      1.00000
    118       5.2551      1.00000
    119       5.2551      1.00000
    120       5.2551      1.00000
    121       5.7418      1.00000
    122       5.7418      1.00000
    123       5.8436      1.00000
    124       5.8436      1.00000
    125       6.0140      1.00000
    126       6.0141      1.00000
    127       7.5376      0.00000
    128       7.5377      0.00000
    129      10.4389      0.00000
    130      10.4389      0.00000
    131      10.4389      0.00000
    132      10.4389      0.00000
    133      10.6464      0.00000
    134      10.6464      0.00000
    135      10.6465      0.00000
    136      10.6465      0.00000
    137      10.7096      0.00000
    138      10.7096      0.00000
    139      10.8754      0.00000
    140      10.8754      0.00000
    141      10.8772      0.00000
    142      10.8772      0.00000
    143      10.8983      0.00000
    144      10.8983      0.00000
    145      10.8983      0.00000
    146      10.8983      0.00000
    147      11.0749      0.00000
    148      11.0749      0.00000
    149      11.0749      0.00000
    150      11.0749      0.00000
    151      11.1448      0.00000
    152      11.1448      0.00000
    153      11.1448      0.00000
    154      11.1448      0.00000
    155      11.1532      0.00000
    156      11.1532      0.00000
    157      11.1865      0.00000
    158      11.1865      0.00000
    159      11.3432      0.00000
    160      11.3432      0.00000
    161      11.3644      0.00000
    162      11.3644      0.00000
    163      11.4013      0.00000
    164      11.4013      0.00000
    165      11.4342      0.00000
    166      11.4342      0.00000
    167      11.4342      0.00000
    168      11.4342      0.00000
    169      11.5581      0.00000
    170      11.5581      0.00000
    171      11.5582      0.00000
    172      11.5582      0.00000
    173      11.5685      0.00000
    174      11.5685      0.00000
    175      11.6573      0.00000
    176      11.6573      0.00000
    177      11.6573      0.00000
    178      11.6573      0.00000
    179      12.2294      0.00000
    180      12.2295      0.00000
    181      12.2295      0.00000
    182      12.2295      0.00000
    183      12.2768      0.00000
    184      12.2768      0.00000
    185      12.7492      0.00000
    186      12.7492      0.00000
    187      12.7492      0.00000
    188      12.7492      0.00000
    189      12.8637      0.00000
    190      12.8637      0.00000
    191      12.9739      0.00000
    192      12.9739      0.00000

 k-point     6 :       0.2500    0.0000    0.5000
  band No.  band energies     occupation 
      1     -23.8131      1.00000
      2     -23.8131      1.00000
      3     -23.6515      1.00000
      4     -23.6515      1.00000
      5     -23.6179      1.00000
      6     -23.6179      1.00000
      7     -14.5780      1.00000
      8     -14.5779      1.00000
      9     -14.5407      1.00000
     10     -14.5406      1.00000
     11     -14.5405      1.00000
     12     -14.5405      1.00000
     13     -14.4179      1.00000
     14     -14.4178      1.00000
     15     -14.4127      1.00000
     16     -14.4127      1.00000
     17     -14.3822      1.00000
     18     -14.3821      1.00000
     19     -14.3586      1.00000
     20     -14.3586      1.00000
     21     -14.3080      1.00000
     22     -14.3080      1.00000
     23     -14.3009      1.00000
     24     -14.3008      1.00000
     25     -14.2480      1.00000
     26     -14.2479      1.00000
     27     -14.2445      1.00000
     28     -14.2444      1.00000
     29     -14.2327      1.00000
     30     -14.2327      1.00000
     31     -14.2314      1.00000
     32     -14.2313      1.00000
     33     -14.2265      1.00000
     34     -14.2265      1.00000
     35     -14.2223      1.00000
     36     -14.2222      1.00000
     37     -11.6052      1.00000
     38     -11.6052      1.00000
     39     -11.4679      1.00000
     40     -11.4679      1.00000
     41     -11.1477      1.00000
     42     -11.1476      1.00000
     43     -10.9277      1.00000
     44     -10.9276      1.00000
     45     -10.7728      1.00000
     46     -10.7728      1.00000
     47     -10.7123      1.00000
     48     -10.7123      1.00000
     49      -7.9164      1.00000
     50      -7.9164      1.00000
     51      -7.8351      1.00000
     52      -7.8351      1.00000
     53      -7.4184      1.00000
     54      -7.4184      1.00000
     55      -7.3195      1.00000
     56      -7.3195      1.00000
     57      -7.2099      1.00000
     58      -7.2099      1.00000
     59      -7.1409      1.00000
     60      -7.1409      1.00000
     61      -6.8138      1.00000
     62      -6.8138      1.00000
     63      -6.7339      1.00000
     64      -6.7339      1.00000
     65      -6.6061      1.00000
     66      -6.6061      1.00000
     67      -6.4404      1.00000
     68      -6.4404      1.00000
     69      -5.9590      1.00000
     70      -5.9590      1.00000
     71      -5.8806      1.00000
     72      -5.8806      1.00000
     73      -0.4326      1.00000
     74      -0.4326      1.00000
     75      -0.2088      1.00000
     76      -0.2088      1.00000
     77       1.6114      1.00000
     78       1.6115      1.00000
     79       2.1275      1.00000
     80       2.1275      1.00000
     81       2.2653      1.00000
     82       2.2653      1.00000
     83       2.3690      1.00000
     84       2.3690      1.00000
     85       2.7845      1.00000
     86       2.7845      1.00000
     87       2.8475      1.00000
     88       2.8475      1.00000
     89       3.2237      1.00000
     90       3.2238      1.00000
     91       3.3603      1.00000
     92       3.3603      1.00000
     93       3.6095      1.00000
     94       3.6095      1.00000
     95       3.7188      1.00000
     96       3.7188      1.00000
     97       4.0668      1.00000
     98       4.0668      1.00000
     99       4.1680      1.00000
    100       4.1680      1.00000
    101       4.3215      1.00000
    102       4.3215      1.00000
    103       4.4303      1.00000
    104       4.4303      1.00000
    105       4.5799      1.00000
    106       4.5799      1.00000
    107       4.6570      1.00000
    108       4.6570      1.00000
    109       4.7273      1.00000
    110       4.7274      1.00000
    111       4.8778      1.00000
    112       4.8778      1.00000
    113       4.9534      1.00000
    114       4.9534      1.00000
    115       5.0933      1.00000
    116       5.0934      1.00000
    117       5.1632      1.00000
    118       5.1632      1.00000
    119       5.1854      1.00000
    120       5.1854      1.00000
    121       5.8611      1.00000
    122       5.8611      1.00000
    123       5.9919      1.00000
    124       5.9920      1.00000
    125       6.0218      1.00000
    126       6.0218      1.00000
    127       8.1600      0.00000
    128       8.1600      0.00000
    129       9.9488      0.00000
    130       9.9488      0.00000
    131      10.1823      0.00000
    132      10.1823      0.00000
    133      10.5319      0.00000
    134      10.5319      0.00000
    135      10.6610      0.00000
    136      10.6611      0.00000
    137      10.7626      0.00000
    138      10.7626      0.00000
    139      10.8142      0.00000
    140      10.8142      0.00000
    141      10.8885      0.00000
    142      10.8885      0.00000
    143      10.9000      0.00000
    144      10.9000      0.00000
    145      10.9338      0.00000
    146      10.9338      0.00000
    147      11.0441      0.00000
    148      11.0441      0.00000
    149      11.0572      0.00000
    150      11.0572      0.00000
    151      11.1097      0.00000
    152      11.1097      0.00000
    153      11.1354      0.00000
    154      11.1354      0.00000
    155      11.1571      0.00000
    156      11.1571      0.00000
    157      11.1763      0.00000
    158      11.1763      0.00000
    159      11.2076      0.00000
    160      11.2076      0.00000
    161      11.2419      0.00000
    162      11.2419      0.00000
    163      11.2723      0.00000
    164      11.2723      0.00000
    165      11.3099      0.00000
    166      11.3099      0.00000
    167      11.3729      0.00000
    168      11.3729      0.00000
    169      11.4199      0.00000
    170      11.4200      0.00000
    171      11.4796      0.00000
    172      11.4796      0.00000
    173      11.5222      0.00000
    174      11.5222      0.00000
    175      11.7459      0.00000
    176      11.7459      0.00000
    177      11.8351      0.00000
    178      11.8351      0.00000
    179      11.9274      0.00000
    180      11.9274      0.00000
    181      12.1503      0.00000
    182      12.1503      0.00000
    183      12.2502      0.00000
    184      12.2502      0.00000
    185      12.3731      0.00000
    186      12.3731      0.00000
    187      12.5266      0.00000
    188      12.5266      0.00000
    189      12.5870      0.00000
    190      12.5870      0.00000
    191      12.8459      0.00000
    192      12.8459      0.00000

 k-point     7 :       0.5000    0.0000    0.5000
  band No.  band energies     occupation 
      1     -23.7381      1.00000
      2     -23.7381      1.00000
      3     -23.7119      1.00000
      4     -23.7119      1.00000
      5     -23.6352      1.00000
      6     -23.6352      1.00000
      7     -14.5608      1.00000
      8     -14.5607      1.00000
      9     -14.5597      1.00000
     10     -14.5596      1.00000
     11     -14.5513      1.00000
     12     -14.5512      1.00000
     13     -14.4179      1.00000
     14     -14.4179      1.00000
     15     -14.4001      1.00000
     16     -14.4000      1.00000
     17     -14.3998      1.00000
     18     -14.3998      1.00000
     19     -14.3280      1.00000
     20     -14.3279      1.00000
     21     -14.3265      1.00000
     22     -14.3264      1.00000
     23     -14.2937      1.00000
     24     -14.2936      1.00000
     25     -14.2499      1.00000
     26     -14.2499      1.00000
     27     -14.2384      1.00000
     28     -14.2383      1.00000
     29     -14.2374      1.00000
     30     -14.2373      1.00000
     31     -14.2276      1.00000
     32     -14.2275      1.00000
     33     -14.2269      1.00000
     34     -14.2269      1.00000
     35     -14.2252      1.00000
     36     -14.2251      1.00000
     37     -11.7282      1.00000
     38     -11.7282      1.00000
     39     -11.2654      1.00000
     40     -11.2654      1.00000
     41     -11.2093      1.00000
     42     -11.2093      1.00000
     43     -10.9279      1.00000
     44     -10.9279      1.00000
     45     -10.8999      1.00000
     46     -10.8999      1.00000
     47     -10.4999      1.00000
     48     -10.4999      1.00000
     49      -7.9238      1.00000
     50      -7.9238      1.00000
     51      -7.6790      1.00000
     52      -7.6790      1.00000
     53      -7.6742      1.00000
     54      -7.6742      1.00000
     55      -7.6032      1.00000
     56      -7.6032      1.00000
     57      -7.0982      1.00000
     58      -7.0982      1.00000
     59      -7.0550      1.00000
     60      -7.0550      1.00000
     61      -6.9209      1.00000
     62      -6.9209      1.00000
     63      -6.8467      1.00000
     64      -6.8467      1.00000
     65      -6.3869      1.00000
     66      -6.3869      1.00000
     67      -6.1985      1.00000
     68      -6.1985      1.00000
     69      -6.1053      1.00000
     70      -6.1052      1.00000
     71      -5.9570      1.00000
     72      -5.9570      1.00000
     73      -0.5072      1.00000
     74      -0.5072      1.00000
     75       0.5746      1.00000
     76       0.5746      1.00000
     77       0.5750      1.00000
     78       0.5750      1.00000
     79       2.3668      1.00000
     80       2.3668      1.00000
     81       2.3831      1.00000
     82       2.3831      1.00000
     83       2.6161      1.00000
     84       2.6161      1.00000
     85       2.9791      1.00000
     86       2.9791      1.00000
     87       3.0927      1.00000
     88       3.0927      1.00000
     89       3.1004      1.00000
     90       3.1004      1.00000
     91       3.4352      1.00000
     92       3.4352      1.00000
     93       3.8398      1.00000
     94       3.8398      1.00000
     95       3.9001      1.00000
     96       3.9001      1.00000
     97       3.9504      1.00000
     98       3.9504      1.00000
     99       3.9952      1.00000
    100       3.9953      1.00000
    101       4.1102      1.00000
    102       4.1102      1.00000
    103       4.1287      1.00000
    104       4.1287      1.00000
    105       4.2575      1.00000
    106       4.2575      1.00000
    107       4.4710      1.00000
    108       4.4710      1.00000
    109       4.8151      1.00000
    110       4.8151      1.00000
    111       4.8232      1.00000
    112       4.8232      1.00000
    113       4.9678      1.00000
    114       4.9678      1.00000
    115       5.1517      1.00000
    116       5.1517      1.00000
    117       5.1578      1.00000
    118       5.1578      1.00000
    119       5.2724      1.00000
    120       5.2724      1.00000
    121       5.9109      1.00000
    122       5.9109      1.00000
    123       6.0005      1.00000
    124       6.0006      1.00000
    125       6.0414      1.00000
    126       6.0414      1.00000
    127       9.2593      0.00000
    128       9.2593      0.00000
    129       9.3171      0.00000
    130       9.3171      0.00000
    131       9.3380      0.00000
    132       9.3380      0.00000
    133      10.3679      0.00000
    134      10.3679      0.00000
    135      10.3991      0.00000
    136      10.3991      0.00000
    137      10.7938      0.00000
    138      10.7938      0.00000
    139      10.8682      0.00000
    140      10.8682      0.00000
    141      10.8698      0.00000
    142      10.8698      0.00000
    143      10.9451      0.00000
    144      10.9451      0.00000
    145      10.9622      0.00000
    146      10.9622      0.00000
    147      11.0480      0.00000
    148      11.0480      0.00000
    149      11.0524      0.00000
    150      11.0524      0.00000
    151      11.0748      0.00000
    152      11.0749      0.00000
    153      11.0964      0.00000
    154      11.0964      0.00000
    155      11.1321      0.00000
    156      11.1321      0.00000
    157      11.1365      0.00000
    158      11.1365      0.00000
    159      11.1579      0.00000
    160      11.1579      0.00000
    161      11.1747      0.00000
    162      11.1747      0.00000
    163      11.2604      0.00000
    164      11.2604      0.00000
    165      11.2814      0.00000
    166      11.2814      0.00000
    167      11.3197      0.00000
    168      11.3197      0.00000
    169      11.4439      0.00000
    170      11.4439      0.00000
    171      11.4822      0.00000
    172      11.4822      0.00000
    173      11.4945      0.00000
    174      11.4945      0.00000
    175      11.5313      0.00000
    176      11.5313      0.00000
    177      11.7804      0.00000
    178      11.7804      0.00000
    179      11.7917      0.00000
    180      11.7917      0.00000
    181      12.1489      0.00000
    182      12.1489      0.00000
    183      12.2106      0.00000
    184      12.2106      0.00000
    185      12.3096      0.00000
    186      12.3097      0.00000
    187      12.4747      0.00000
    188      12.4747      0.00000
    189      12.5963      0.00000
    190      12.5963      0.00000
    191      12.9535      0.00000
    192      12.9535      0.00000

 k-point     8 :       0.2500    0.2500    0.5000
  band No.  band energies     occupation 
      1     -23.7527      1.00000
      2     -23.7527      1.00000
      3     -23.6840      1.00000
      4     -23.6840      1.00000
      5     -23.6485      1.00000
      6     -23.6485      1.00000
      7     -14.5665      1.00000
      8     -14.5665      1.00000
      9     -14.5527      1.00000
     10     -14.5527      1.00000
     11     -14.5526      1.00000
     12     -14.5525      1.00000
     13     -14.4126      1.00000
     14     -14.4125      1.00000
     15     -14.4104      1.00000
     16     -14.4103      1.00000
     17     -14.3959      1.00000
     18     -14.3958      1.00000
     19     -14.3339      1.00000
     20     -14.3338      1.00000
     21     -14.3069      1.00000
     22     -14.3068      1.00000
     23     -14.3061      1.00000
     24     -14.3060      1.00000
     25     -14.2462      1.00000
     26     -14.2461      1.00000
     27     -14.2447      1.00000
     28     -14.2446      1.00000
     29     -14.2356      1.00000
     30     -14.2356      1.00000
     31     -14.2279      1.00000
     32     -14.2279      1.00000
     33     -14.2270      1.00000
     34     -14.2270      1.00000
     35     -14.2239      1.00000
     36     -14.2238      1.00000
     37     -11.5663      1.00000
     38     -11.5663      1.00000
     39     -11.5440      1.00000
     40     -11.5440      1.00000
     41     -11.2163      1.00000
     42     -11.2163      1.00000
     43     -10.8521      1.00000
     44     -10.8520      1.00000
     45     -10.6893      1.00000
     46     -10.6893      1.00000
     47     -10.6749      1.00000
     48     -10.6749      1.00000
     49      -7.8222      1.00000
     50      -7.8222      1.00000
     51      -7.7757      1.00000
     52      -7.7757      1.00000
     53      -7.4648      1.00000
     54      -7.4648      1.00000
     55      -7.4441      1.00000
     56      -7.4441      1.00000
     57      -7.4041      1.00000
     58      -7.4041      1.00000
     59      -7.3548      1.00000
     60      -7.3548      1.00000
     61      -6.8053      1.00000
     62      -6.8053      1.00000
     63      -6.6033      1.00000
     64      -6.6033      1.00000
     65      -6.5593      1.00000
     66      -6.5593      1.00000
     67      -6.2058      1.00000
     68      -6.2057      1.00000
     69      -6.0317      1.00000
     70      -6.0317      1.00000
     71      -5.9677      1.00000
     72      -5.9676      1.00000
     73      -0.1210      1.00000
     74      -0.1210      1.00000
     75      -0.1114      1.00000
     76      -0.1114      1.00000
     77       0.8142      1.00000
     78       0.8142      1.00000
     79       2.2811      1.00000
     80       2.2812      1.00000
     81       2.5190      1.00000
     82       2.5190      1.00000
     83       2.5971      1.00000
     84       2.5971      1.00000
     85       2.9670      1.00000
     86       2.9670      1.00000
     87       3.0182      1.00000
     88       3.0182      1.00000
     89       3.2860      1.00000
     90       3.2860      1.00000
     91       3.5838      1.00000
     92       3.5838      1.00000
     93       3.5981      1.00000
     94       3.5981      1.00000
     95       3.7050      1.00000
     96       3.7050      1.00000
     97       3.9771      1.00000
     98       3.9772      1.00000
     99       4.0546      1.00000
    100       4.0546      1.00000
    101       4.0874      1.00000
    102       4.0874      1.00000
    103       4.2091      1.00000
    104       4.2091      1.00000
    105       4.3874      1.00000
    106       4.3874      1.00000
    107       4.4348      1.00000
    108       4.4348      1.00000
    109       4.9146      1.00000
    110       4.9146      1.00000
    111       4.9492      1.00000
    112       4.9492      1.00000
    113       4.9697      1.00000
    114       4.9697      1.00000
    115       5.0809      1.00000
    116       5.0809      1.00000
    117       5.1395      1.00000
    118       5.1395      1.00000
    119       5.1930      1.00000
    120       5.1930      1.00000
    121       5.9116      1.00000
    122       5.9116      1.00000
    123       6.0030      1.00000
    124       6.0031      1.00000
    125       6.0133      1.00000
    126       6.0133      1.00000
    127       8.8591      0.00000
    128       8.8591      0.00000
    129       9.4082      0.00000
    130       9.4082      0.00000
    131       9.4620      0.00000
    132       9.4620      0.00000
    133      10.4005      0.00000
    134      10.4005      0.00000
    135      10.5450      0.00000
    136      10.5450      0.00000
    137      10.5656      0.00000
    138      10.5656      0.00000
    139      10.8932      0.00000
    140      10.8932      0.00000
    141      10.8980      0.00000
    142      10.8980      0.00000
    143      10.9636      0.00000
    144      10.9636      0.00000
    145      10.9955      0.00000
    146      10.9955      0.00000
    147      11.0172      0.00000
    148      11.0172      0.00000
    149      11.0648      0.00000
    150      11.0648      0.00000
    151      11.1128      0.00000
    152      11.1128      0.00000
    153      11.1154      0.00000
    154      11.1154      0.00000
    155      11.1325      0.00000
    156      11.1325      0.00000
    157      11.1451      0.00000
    158      11.1451      0.00000
    159      11.1617      0.00000
    160      11.1617      0.00000
    161      11.1686      0.00000
    162      11.1686      0.00000
    163      11.2759      0.00000
    164      11.2759      0.00000
    165      11.2831      0.00000
    166      11.2831      0.00000
    167      11.3619      0.00000
    168      11.3619      0.00000
    169      11.4641      0.00000
    170      11.4641      0.00000
    171      11.4676      0.00000
    172      11.4677      0.00000
    173      11.4968      0.00000
    174      11.4968      0.00000
    175      11.6202      0.00000
    176      11.6202      0.00000
    177      11.6863      0.00000
    178      11.6863      0.00000
    179      11.7210      0.00000
    180      11.7210      0.00000
    181      12.1568      0.00000
    182      12.1568      0.00000
    183      12.2263      0.00000
    184      12.2263      0.00000
    185      12.2275      0.00000
    186      12.2275      0.00000
    187      12.6266      0.00000
    188      12.6266      0.00000
    189      12.7126      0.00000
    190      12.7127      0.00000
    191      13.0342      0.00000
    192      13.0342      0.00000

 spin component 2

 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -23.8511      1.00000
      2     -23.8504      1.00000
      3     -23.6291      1.00000
      4     -23.6291      1.00000
      5     -23.5997      1.00000
      6     -23.5996      1.00000
      7     -14.5958      1.00000
      8     -14.5761      1.00000
      9     -14.5303      1.00000
     10     -14.5302      1.00000
     11     -14.5299      1.00000
     12     -14.5298      1.00000
     13     -14.4195      1.00000
     14     -14.4191      1.00000
     15     -14.4146      1.00000
     16     -14.4141      1.00000
     17     -14.3738      1.00000
     18     -14.3738      1.00000
     19     -14.3738      1.00000
     20     -14.3738      1.00000
     21     -14.3072      1.00000
     22     -14.3072      1.00000
     23     -14.3071      1.00000
     24     -14.3071      1.00000
     25     -14.2465      1.00000
     26     -14.2451      1.00000
     27     -14.2433      1.00000
     28     -14.2416      1.00000
     29     -14.2378      1.00000
     30     -14.2377      1.00000
     31     -14.2377      1.00000
     32     -14.2377      1.00000
     33     -14.2212      1.00000
     34     -14.2212      1.00000
     35     -14.2211      1.00000
     36     -14.2211      1.00000
     37     -11.3955      1.00000
     38     -11.3954      1.00000
     39     -11.3812      1.00000
     40     -11.3810      1.00000
     41     -11.2561      1.00000
     42     -11.2208      1.00000
     43     -11.1303      1.00000
     44     -11.1302      1.00000
     45     -11.1156      1.00000
     46     -11.1154      1.00000
     47     -10.4841      1.00000
     48     -10.4433      1.00000
     49      -7.9194      1.00000
     50      -7.8404      1.00000
     51      -7.8123      1.00000
     52      -7.7145      1.00000
     53      -7.2880      1.00000
     54      -7.2880      1.00000
     55      -7.2867      1.00000
     56      -7.2867      1.00000
     57      -7.0363      1.00000
     58      -7.0362      1.00000
     59      -7.0315      1.00000
     60      -7.0315      1.00000
     61      -6.9079      1.00000
     62      -6.8972      1.00000
     63      -6.8972      1.00000
     64      -6.8923      1.00000
     65      -6.8923      1.00000
     66      -6.7583      1.00000
     67      -6.5890      1.00000
     68      -6.4637      1.00000
     69      -5.9545      1.00000
     70      -5.8300      1.00000
     71      -5.7918      1.00000
     72      -5.7708      1.00000
     73      -0.3939      1.00000
     74      -0.3938      1.00000
     75      -0.3787      1.00000
     76      -0.3786      1.00000
     77       1.6667      1.00000
     78       1.8573      1.00000
     79       1.8573      1.00000
     80       1.8586      1.00000
     81       1.8586      1.00000
     82       1.9792      1.00000
     83       2.4265      1.00000
     84       2.5322      1.00000
     85       2.7460      1.00000
     86       2.7461      1.00000
     87       2.9729      1.00000
     88       2.9730      1.00000
     89       3.1821      1.00000
     90       3.1821      1.00000
     91       3.1845      1.00000
     92       3.1845      1.00000
     93       3.6463      1.00000
     94       3.6660      1.00000
     95       3.9251      1.00000
     96       3.9436      1.00000
     97       4.0914      1.00000
     98       4.0915      1.00000
     99       4.1461      1.00000
    100       4.1462      1.00000
    101       4.3700      1.00000
    102       4.3700      1.00000
    103       4.6398      1.00000
    104       4.6398      1.00000
    105       4.6638      1.00000
    106       4.6638      1.00000
    107       4.8017      1.00000
    108       4.8076      1.00000
    109       4.8084      1.00000
    110       4.8345      1.00000
    111       4.8528      1.00000
    112       4.8529      1.00000
    113       5.0278      1.00000
    114       5.0278      1.00000
    115       5.0751      1.00000
    116       5.0752      1.00000
    117       5.0906      1.00000
    118       5.0907      1.00000
    119       5.3531      1.00000
    120       5.3531      1.00000
    121       5.6056      1.00000
    122       5.6547      1.00000
    123       5.8477      1.00000
    124       5.9938      1.00000
    125       6.0305      1.00000
    126       6.0596      1.00000
    127       7.2550      0.00000
    128       7.8554      0.00000
    129      10.4236      0.00000
    130      10.4237      0.00000
    131      10.4443      0.00000
    132      10.4444      0.00000
    133      10.6102      0.00000
    134      10.6102      0.00000
    135      10.6712      0.00000
    136      10.6713      0.00000
    137      10.7009      0.00000
    138      10.7233      0.00000
    139      10.7899      0.00000
    140      10.7927      0.00000
    141      10.9227      0.00000
    142      10.9227      0.00000
    143      10.9367      0.00000
    144      10.9367      0.00000
    145      10.9498      0.00000
    146      10.9509      0.00000
    147      11.0314      0.00000
    148      11.0314      0.00000
    149      11.0778      0.00000
    150      11.0778      0.00000
    151      11.1274      0.00000
    152      11.1274      0.00000
    153      11.1304      0.00000
    154      11.1645      0.00000
    155      11.1645      0.00000
    156      11.1722      0.00000
    157      11.1799      0.00000
    158      11.1944      0.00000
    159      11.3082      0.00000
    160      11.3189      0.00000
    161      11.3581      0.00000
    162      11.3814      0.00000
    163      11.3879      0.00000
    164      11.4311      0.00000
    165      11.4311      0.00000
    166      11.4441      0.00000
    167      11.4441      0.00000
    168      11.4518      0.00000
    169      11.5399      0.00000
    170      11.5400      0.00000
    171      11.5643      0.00000
    172      11.5753      0.00000
    173      11.5753      0.00000
    174      11.6251      0.00000
    175      11.6251      0.00000
    176      11.6260      0.00000
    177      11.6818      0.00000
    178      11.6818      0.00000
    179      11.9873      0.00000
    180      12.0090      0.00000
    181      12.2057      0.00000
    182      12.2057      0.00000
    183      12.3514      0.00000
    184      12.3514      0.00000
    185      12.3786      0.00000
    186      12.6658      0.00000
    187      12.6658      0.00000
    188      12.6933      0.00000
    189      12.6933      0.00000
    190      12.7412      0.00000
    191      12.9030      0.00000
    192      12.9041      0.00000

 k-point     2 :       0.2500    0.0000    0.0000
  band No.  band energies     occupation 
      1     -23.8133      1.00000
      2     -23.8127      1.00000
      3     -23.6515      1.00000
      4     -23.6514      1.00000
      5     -23.6178      1.00000
      6     -23.6177      1.00000
      7     -14.5845      1.00000
      8     -14.5716      1.00000
      9     -14.5492      1.00000
     10     -14.5457      1.00000
     11     -14.5355      1.00000
     12     -14.5322      1.00000
     13     -14.4191      1.00000
     14     -14.4167      1.00000
     15     -14.4134      1.00000
     16     -14.4122      1.00000
     17     -14.3827      1.00000
     18     -14.3818      1.00000
     19     -14.3587      1.00000
     20     -14.3587      1.00000
     21     -14.3081      1.00000
     22     -14.3081      1.00000
     23     -14.3010      1.00000
     24     -14.3009      1.00000
     25     -14.2488      1.00000
     26     -14.2473      1.00000
     27     -14.2453      1.00000
     28     -14.2439      1.00000
     29     -14.2329      1.00000
     30     -14.2328      1.00000
     31     -14.2315      1.00000
     32     -14.2314      1.00000
     33     -14.2268      1.00000
     34     -14.2264      1.00000
     35     -14.2224      1.00000
     36     -14.2223      1.00000
     37     -11.6969      1.00000
     38     -11.5488      1.00000
     39     -11.4852      1.00000
     40     -11.3537      1.00000
     41     -11.1704      1.00000
     42     -11.1159      1.00000
     43     -11.0411      1.00000
     44     -10.8955      1.00000
     45     -10.8572      1.00000
     46     -10.7518      1.00000
     47     -10.6792      1.00000
     48     -10.6715      1.00000
     49      -7.9335      1.00000
     50      -7.9100      1.00000
     51      -7.8806      1.00000
     52      -7.7834      1.00000
     53      -7.4189      1.00000
     54      -7.4167      1.00000
     55      -7.3220      1.00000
     56      -7.3124      1.00000
     57      -7.2686      1.00000
     58      -7.1526      1.00000
     59      -7.1452      1.00000
     60      -7.0397      1.00000
     61      -6.9025      1.00000
     62      -6.8250      1.00000
     63      -6.7481      1.00000
     64      -6.7435      1.00000
     65      -6.6141      1.00000
     66      -6.5899      1.00000
     67      -6.5190      1.00000
     68      -6.3321      1.00000
     69      -6.0058      1.00000
     70      -5.9373      1.00000
     71      -5.8917      1.00000
     72      -5.8536      1.00000
     73      -0.4892      1.00000
     74      -0.3740      1.00000
     75      -0.2161      1.00000
     76      -0.2016      1.00000
     77       1.5215      1.00000
     78       1.7916      1.00000
     79       1.8945      1.00000
     80       2.1922      1.00000
     81       2.3182      1.00000
     82       2.3343      1.00000
     83       2.3608      1.00000
     84       2.4129      1.00000
     85       2.6323      1.00000
     86       2.7207      1.00000
     87       2.8451      1.00000
     88       2.9121      1.00000
     89       3.2078      1.00000
     90       3.2678      1.00000
     91       3.3475      1.00000
     92       3.4458      1.00000
     93       3.5248      1.00000
     94       3.6413      1.00000
     95       3.7399      1.00000
     96       3.8924      1.00000
     97       3.9759      1.00000
     98       4.0335      1.00000
     99       4.1503      1.00000
    100       4.2731      1.00000
    101       4.2945      1.00000
    102       4.3540      1.00000
    103       4.4035      1.00000
    104       4.4369      1.00000
    105       4.4449      1.00000
    106       4.5657      1.00000
    107       4.5988      1.00000
    108       4.6544      1.00000
    109       4.7188      1.00000
    110       4.7851      1.00000
    111       4.8047      1.00000
    112       4.8663      1.00000
    113       4.8840      1.00000
    114       4.9517      1.00000
    115       5.0720      1.00000
    116       5.1530      1.00000
    117       5.1608      1.00000
    118       5.2224      1.00000
    119       5.2576      1.00000
    120       5.2764      1.00000
    121       5.8207      1.00000
    122       5.9004      1.00000
    123       5.9532      1.00000
    124       5.9817      1.00000
    125       5.9983      1.00000
    126       6.0445      1.00000
    127       8.0496      0.00000
    128       8.3008      0.00000
    129       9.8406      0.00000
    130      10.0590      0.00000
    131      10.1067      0.00000
    132      10.3476      0.00000
    133      10.4945      0.00000
    134      10.5066      0.00000
    135      10.6201      0.00000
    136      10.6771      0.00000
    137      10.7374      0.00000
    138      10.7510      0.00000
    139      10.7824      0.00000
    140      10.8361      0.00000
    141      10.8563      0.00000
    142      10.8868      0.00000
    143      10.9126      0.00000
    144      10.9271      0.00000
    145      10.9365      0.00000
    146      10.9724      0.00000
    147      11.0404      0.00000
    148      11.0417      0.00000
    149      11.0468      0.00000
    150      11.0573      0.00000
    151      11.0994      0.00000
    152      11.1212      0.00000
    153      11.1252      0.00000
    154      11.1380      0.00000
    155      11.1619      0.00000
    156      11.1627      0.00000
    157      11.1735      0.00000
    158      11.1846      0.00000
    159      11.2058      0.00000
    160      11.2137      0.00000
    161      11.2274      0.00000
    162      11.2382      0.00000
    163      11.2567      0.00000
    164      11.2868      0.00000
    165      11.3108      0.00000
    166      11.3317      0.00000
    167      11.3567      0.00000
    168      11.3568      0.00000
    169      11.3895      0.00000
    170      11.4688      0.00000
    171      11.4705      0.00000
    172      11.5090      0.00000
    173      11.5197      0.00000
    174      11.5365      0.00000
    175      11.6990      0.00000
    176      11.7308      0.00000
    177      11.8039      0.00000
    178      11.8442      0.00000
    179      11.8838      0.00000
    180      11.9929      0.00000
    181      12.0186      0.00000
    182      12.0745      0.00000
    183      12.1493      0.00000
    184      12.3006      0.00000
    185      12.3125      0.00000
    186      12.3242      0.00000
    187      12.4456      0.00000
    188      12.5196      0.00000
    189      12.5320      0.00000
    190      12.6187      0.00000
    191      12.7600      0.00000
    192      12.9063      0.00000

 k-point     3 :       0.5000    0.0000    0.0000
  band No.  band energies     occupation 
      1     -23.7382      1.00000
      2     -23.7378      1.00000
      3     -23.7119      1.00000
      4     -23.7116      1.00000
      5     -23.6352      1.00000
      6     -23.6351      1.00000
      7     -14.5630      1.00000
      8     -14.5620      1.00000
      9     -14.5605      1.00000
     10     -14.5596      1.00000
     11     -14.5590      1.00000
     12     -14.5400      1.00000
     13     -14.4208      1.00000
     14     -14.4147      1.00000
     15     -14.4009      1.00000
     16     -14.4008      1.00000
     17     -14.3994      1.00000
     18     -14.3990      1.00000
     19     -14.3281      1.00000
     20     -14.3280      1.00000
     21     -14.3266      1.00000
     22     -14.3265      1.00000
     23     -14.2938      1.00000
     24     -14.2937      1.00000
     25     -14.2508      1.00000
     26     -14.2492      1.00000
     27     -14.2390      1.00000
     28     -14.2380      1.00000
     29     -14.2379      1.00000
     30     -14.2369      1.00000
     31     -14.2278      1.00000
     32     -14.2275      1.00000
     33     -14.2270      1.00000
     34     -14.2270      1.00000
     35     -14.2253      1.00000
     36     -14.2252      1.00000
     37     -11.8250      1.00000
     38     -11.6143      1.00000
     39     -11.3192      1.00000
     40     -11.2690      1.00000
     41     -11.1917      1.00000
     42     -11.1038      1.00000
     43     -11.0179      1.00000
     44     -10.9543      1.00000
     45     -10.8836      1.00000
     46     -10.8674      1.00000
     47     -10.6147      1.00000
     48     -10.4002      1.00000
     49      -7.9592      1.00000
     50      -7.8944      1.00000
     51      -7.6910      1.00000
     52      -7.6899      1.00000
     53      -7.6619      1.00000
     54      -7.6448      1.00000
     55      -7.6061      1.00000
     56      -7.6009      1.00000
     57      -7.2275      1.00000
     58      -7.1119      1.00000
     59      -7.0271      1.00000
     60      -6.9821      1.00000
     61      -6.9599      1.00000
     62      -6.8637      1.00000
     63      -6.7983      1.00000
     64      -6.7967      1.00000
     65      -6.4008      1.00000
     66      -6.3521      1.00000
     67      -6.3397      1.00000
     68      -6.2830      1.00000
     69      -6.0920      1.00000
     70      -6.0592      1.00000
     71      -5.9470      1.00000
     72      -5.9064      1.00000
     73      -0.6087      1.00000
     74      -0.3942      1.00000
     75       0.5292      1.00000
     76       0.5453      1.00000
     77       0.6052      1.00000
     78       0.6235      1.00000
     79       2.2825      1.00000
     80       2.3260      1.00000
     81       2.4209      1.00000
     82       2.4234      1.00000
     83       2.4379      1.00000
     84       2.6558      1.00000
     85       2.8619      1.00000
     86       2.9736      1.00000
     87       3.2630      1.00000
     88       3.2752      1.00000
     89       3.3250      1.00000
     90       3.3520      1.00000
     91       3.3964      1.00000
     92       3.4709      1.00000
     93       3.5349      1.00000
     94       3.5588      1.00000
     95       3.7089      1.00000
     96       3.7352      1.00000
     97       3.7412      1.00000
     98       3.9743      1.00000
     99       3.9919      1.00000
    100       4.0564      1.00000
    101       4.0651      1.00000
    102       4.2116      1.00000
    103       4.3274      1.00000
    104       4.3484      1.00000
    105       4.4744      1.00000
    106       4.4802      1.00000
    107       4.4954      1.00000
    108       4.5529      1.00000
    109       4.6029      1.00000
    110       4.6621      1.00000
    111       4.7548      1.00000
    112       4.8159      1.00000
    113       4.8481      1.00000
    114       4.9322      1.00000
    115       5.0111      1.00000
    116       5.1638      1.00000
    117       5.1669      1.00000
    118       5.3073      1.00000
    119       5.3310      1.00000
    120       5.3532      1.00000
    121       5.8214      1.00000
    122       5.8870      1.00000
    123       5.9598      1.00000
    124       5.9715      1.00000
    125       6.0970      1.00000
    126       6.1331      1.00000
    127       9.2209      0.00000
    128       9.2621      0.00000
    129       9.3052      0.00000
    130       9.3234      0.00000
    131       9.3403      0.00000
    132       9.3645      0.00000
    133      10.2918      0.00000
    134      10.3345      0.00000
    135      10.4509      0.00000
    136      10.4725      0.00000
    137      10.7854      0.00000
    138      10.8123      0.00000
    139      10.8451      0.00000
    140      10.8650      0.00000
    141      10.8729      0.00000
    142      10.9295      0.00000
    143      10.9477      0.00000
    144      10.9494      0.00000
    145      10.9662      0.00000
    146      10.9747      0.00000
    147      11.0065      0.00000
    148      11.0234      0.00000
    149      11.0419      0.00000
    150      11.0671      0.00000
    151      11.0678      0.00000
    152      11.0852      0.00000
    153      11.1011      0.00000
    154      11.1324      0.00000
    155      11.1329      0.00000
    156      11.1444      0.00000
    157      11.1470      0.00000
    158      11.1510      0.00000
    159      11.1518      0.00000
    160      11.1619      0.00000
    161      11.1684      0.00000
    162      11.1698      0.00000
    163      11.1924      0.00000
    164      11.2492      0.00000
    165      11.2515      0.00000
    166      11.2883      0.00000
    167      11.2925      0.00000
    168      11.3214      0.00000
    169      11.4311      0.00000
    170      11.4448      0.00000
    171      11.4619      0.00000
    172      11.4858      0.00000
    173      11.4986      0.00000
    174      11.5008      0.00000
    175      11.5146      0.00000
    176      11.5831      0.00000
    177      11.6818      0.00000
    178      11.7221      0.00000
    179      11.8276      0.00000
    180      11.8318      0.00000
    181      11.9857      0.00000
    182      12.0271      0.00000
    183      12.2350      0.00000
    184      12.2397      0.00000
    185      12.2586      0.00000
    186      12.3153      0.00000
    187      12.5262      0.00000
    188      12.5746      0.00000
    189      12.6485      0.00000
    190      12.6555      0.00000
    191      12.9264      0.00000
    192      13.0035      0.00000

 k-point     4 :       0.2500    0.2500    0.0000
  band No.  band energies     occupation 
      1     -23.7528      1.00000
      2     -23.7524      1.00000
      3     -23.6840      1.00000
      4     -23.6838      1.00000
      5     -23.6485      1.00000
      6     -23.6483      1.00000
      7     -14.5696      1.00000
      8     -14.5636      1.00000
      9     -14.5607      1.00000
     10     -14.5605      1.00000
     11     -14.5448      1.00000
     12     -14.5448      1.00000
     13     -14.4142      1.00000
     14     -14.4120      1.00000
     15     -14.4109      1.00000
     16     -14.4087      1.00000
     17     -14.3970      1.00000
     18     -14.3949      1.00000
     19     -14.3340      1.00000
     20     -14.3339      1.00000
     21     -14.3070      1.00000
     22     -14.3069      1.00000
     23     -14.3062      1.00000
     24     -14.3061      1.00000
     25     -14.2470      1.00000
     26     -14.2456      1.00000
     27     -14.2455      1.00000
     28     -14.2440      1.00000
     29     -14.2361      1.00000
     30     -14.2352      1.00000
     31     -14.2281      1.00000
     32     -14.2279      1.00000
     33     -14.2272      1.00000
     34     -14.2270      1.00000
     35     -14.2240      1.00000
     36     -14.2238      1.00000
     37     -11.6523      1.00000
     38     -11.6316      1.00000
     39     -11.4609      1.00000
     40     -11.4349      1.00000
     41     -11.2722      1.00000
     42     -11.1383      1.00000
     43     -10.9188      1.00000
     44     -10.8092      1.00000
     45     -10.7924      1.00000
     46     -10.7714      1.00000
     47     -10.6082      1.00000
     48     -10.5964      1.00000
     49      -7.8380      1.00000
     50      -7.8161      1.00000
     51      -7.8094      1.00000
     52      -7.7370      1.00000
     53      -7.4734      1.00000
     54      -7.4721      1.00000
     55      -7.4584      1.00000
     56      -7.4361      1.00000
     57      -7.4106      1.00000
     58      -7.3951      1.00000
     59      -7.3666      1.00000
     60      -7.2974      1.00000
     61      -6.8743      1.00000
     62      -6.7797      1.00000
     63      -6.6423      1.00000
     64      -6.5396      1.00000
     65      -6.5125      1.00000
     66      -6.4645      1.00000
     67      -6.4327      1.00000
     68      -6.1049      1.00000
     69      -6.0886      1.00000
     70      -6.0542      1.00000
     71      -5.9745      1.00000
     72      -5.8973      1.00000
     73      -0.1965      1.00000
     74      -0.1803      1.00000
     75      -0.0365      1.00000
     76      -0.0346      1.00000
     77       0.7326      1.00000
     78       0.9105      1.00000
     79       2.1150      1.00000
     80       2.3386      1.00000
     81       2.3745      1.00000
     82       2.4356      1.00000
     83       2.6168      1.00000
     84       2.7526      1.00000
     85       2.9063      1.00000
     86       3.1174      1.00000
     87       3.1430      1.00000
     88       3.1571      1.00000
     89       3.2814      1.00000
     90       3.3656      1.00000
     91       3.4567      1.00000
     92       3.5507      1.00000
     93       3.5678      1.00000
     94       3.7153      1.00000
     95       3.7250      1.00000
     96       3.7426      1.00000
     97       3.7610      1.00000
     98       3.8260      1.00000
     99       3.8480      1.00000
    100       4.0106      1.00000
    101       4.0485      1.00000
    102       4.0687      1.00000
    103       4.2484      1.00000
    104       4.4149      1.00000
    105       4.4636      1.00000
    106       4.5010      1.00000
    107       4.5636      1.00000
    108       4.6711      1.00000
    109       4.6990      1.00000
    110       4.7958      1.00000
    111       4.8002      1.00000
    112       4.8408      1.00000
    113       4.8751      1.00000
    114       4.8831      1.00000
    115       5.1352      1.00000
    116       5.1551      1.00000
    117       5.1628      1.00000
    118       5.2655      1.00000
    119       5.2682      1.00000
    120       5.2737      1.00000
    121       5.8161      1.00000
    122       5.9078      1.00000
    123       5.9177      1.00000
    124       6.0180      1.00000
    125       6.0678      1.00000
    126       6.0872      1.00000
    127       8.8011      0.00000
    128       8.9094      0.00000
    129       9.3638      0.00000
    130       9.4111      0.00000
    131       9.4368      0.00000
    132       9.4975      0.00000
    133      10.3492      0.00000
    134      10.4443      0.00000
    135      10.5258      0.00000
    136      10.5416      0.00000
    137      10.5798      0.00000
    138      10.6224      0.00000
    139      10.8850      0.00000
    140      10.8882      0.00000
    141      10.8976      0.00000
    142      10.9074      0.00000
    143      10.9529      0.00000
    144      10.9653      0.00000
    145      10.9783      0.00000
    146      10.9972      0.00000
    147      11.0232      0.00000
    148      11.0424      0.00000
    149      11.0580      0.00000
    150      11.0861      0.00000
    151      11.0995      0.00000
    152      11.1102      0.00000
    153      11.1129      0.00000
    154      11.1196      0.00000
    155      11.1312      0.00000
    156      11.1413      0.00000
    157      11.1504      0.00000
    158      11.1603      0.00000
    159      11.1617      0.00000
    160      11.1694      0.00000
    161      11.1734      0.00000
    162      11.1823      0.00000
    163      11.2415      0.00000
    164      11.2598      0.00000
    165      11.2762      0.00000
    166      11.2822      0.00000
    167      11.3274      0.00000
    168      11.3770      0.00000
    169      11.4491      0.00000
    170      11.4561      0.00000
    171      11.4721      0.00000
    172      11.4971      0.00000
    173      11.5039      0.00000
    174      11.5131      0.00000
    175      11.5741      0.00000
    176      11.6421      0.00000
    177      11.6426      0.00000
    178      11.7019      0.00000
    179      11.7071      0.00000
    180      11.8190      0.00000
    181      11.9347      0.00000
    182      12.0111      0.00000
    183      12.0118      0.00000
    184      12.2329      0.00000
    185      12.3917      0.00000
    186      12.3984      0.00000
    187      12.5067      0.00000
    188      12.5274      0.00000
    189      12.7870      0.00000
    190      12.7904      0.00000
    191      12.8630      0.00000
    192      12.9787      0.00000

 k-point     5 :       0.0000    0.0000    0.5000
  band No.  band energies     occupation 
      1     -23.8508      1.00000
      2     -23.8508      1.00000
      3     -23.6291      1.00000
      4     -23.6291      1.00000
      5     -23.5997      1.00000
      6     -23.5996      1.00000
      7     -14.5860      1.00000
      8     -14.5859      1.00000
      9     -14.5301      1.00000
     10     -14.5301      1.00000
     11     -14.5301      1.00000
     12     -14.5300      1.00000
     13     -14.4194      1.00000
     14     -14.4193      1.00000
     15     -14.4144      1.00000
     16     -14.4143      1.00000
     17     -14.3738      1.00000
     18     -14.3738      1.00000
     19     -14.3738      1.00000
     20     -14.3738      1.00000
     21     -14.3072      1.00000
     22     -14.3072      1.00000
     23     -14.3071      1.00000
     24     -14.3071      1.00000
     25     -14.2458      1.00000
     26     -14.2458      1.00000
     27     -14.2425      1.00000
     28     -14.2424      1.00000
     29     -14.2377      1.00000
     30     -14.2377      1.00000
     31     -14.2377      1.00000
     32     -14.2377      1.00000
     33     -14.2212      1.00000
     34     -14.2212      1.00000
     35     -14.2211      1.00000
     36     -14.2211      1.00000
     37     -11.3888      1.00000
     38     -11.3888      1.00000
     39     -11.3887      1.00000
     40     -11.3886      1.00000
     41     -11.2392      1.00000
     42     -11.2392      1.00000
     43     -11.1225      1.00000
     44     -11.1225      1.00000
     45     -11.1224      1.00000
     46     -11.1223      1.00000
     47     -10.4626      1.00000
     48     -10.4626      1.00000
     49      -7.8860      1.00000
     50      -7.8860      1.00000
     51      -7.7653      1.00000
     52      -7.7653      1.00000
     53      -7.2874      1.00000
     54      -7.2874      1.00000
     55      -7.2874      1.00000
     56      -7.2874      1.00000
     57      -7.0344      1.00000
     58      -7.0343      1.00000
     59      -7.0343      1.00000
     60      -7.0343      1.00000
     61      -6.8944      1.00000
     62      -6.8944      1.00000
     63      -6.8944      1.00000
     64      -6.8944      1.00000
     65      -6.8168      1.00000
     66      -6.8168      1.00000
     67      -6.5423      1.00000
     68      -6.5423      1.00000
     69      -5.8717      1.00000
     70      -5.8716      1.00000
     71      -5.7953      1.00000
     72      -5.7953      1.00000
     73      -0.3848      1.00000
     74      -0.3848      1.00000
     75      -0.3847      1.00000
     76      -0.3847      1.00000
     77       1.8131      1.00000
     78       1.8131      1.00000
     79       1.8442      1.00000
     80       1.8442      1.00000
     81       1.8443      1.00000
     82       1.8443      1.00000
     83       2.4834      1.00000
     84       2.4834      1.00000
     85       2.8740      1.00000
     86       2.8740      1.00000
     87       2.8741      1.00000
     88       2.8741      1.00000
     89       3.1838      1.00000
     90       3.1838      1.00000
     91       3.1839      1.00000
     92       3.1839      1.00000
     93       3.7887      1.00000
     94       3.7887      1.00000
     95       3.7938      1.00000
     96       3.7938      1.00000
     97       4.0996      1.00000
     98       4.0996      1.00000
     99       4.0996      1.00000
    100       4.0997      1.00000
    101       4.5449      1.00000
    102       4.5449      1.00000
    103       4.5449      1.00000
    104       4.5449      1.00000
    105       4.6888      1.00000
    106       4.6888      1.00000
    107       4.6888      1.00000
    108       4.6888      1.00000
    109       4.8011      1.00000
    110       4.8011      1.00000
    111       4.8214      1.00000
    112       4.8214      1.00000
    113       5.0641      1.00000
    114       5.0641      1.00000
    115       5.0642      1.00000
    116       5.0642      1.00000
    117       5.2551      1.00000
    118       5.2551      1.00000
    119       5.2551      1.00000
    120       5.2551      1.00000
    121       5.7419      1.00000
    122       5.7419      1.00000
    123       5.8436      1.00000
    124       5.8436      1.00000
    125       6.0141      1.00000
    126       6.0141      1.00000
    127       7.5378      0.00000
    128       7.5378      0.00000
    129      10.4391      0.00000
    130      10.4391      0.00000
    131      10.4391      0.00000
    132      10.4391      0.00000
    133      10.6466      0.00000
    134      10.6467      0.00000
    135      10.6467      0.00000
    136      10.6467      0.00000
    137      10.7098      0.00000
    138      10.7098      0.00000
    139      10.8756      0.00000
    140      10.8756      0.00000
    141      10.8773      0.00000
    142      10.8774      0.00000
    143      10.8984      0.00000
    144      10.8984      0.00000
    145      10.8984      0.00000
    146      10.8984      0.00000
    147      11.0750      0.00000
    148      11.0750      0.00000
    149      11.0750      0.00000
    150      11.0750      0.00000
    151      11.1449      0.00000
    152      11.1449      0.00000
    153      11.1449      0.00000
    154      11.1449      0.00000
    155      11.1533      0.00000
    156      11.1533      0.00000
    157      11.1866      0.00000
    158      11.1866      0.00000
    159      11.3433      0.00000
    160      11.3433      0.00000
    161      11.3645      0.00000
    162      11.3645      0.00000
    163      11.4015      0.00000
    164      11.4015      0.00000
    165      11.4343      0.00000
    166      11.4343      0.00000
    167      11.4343      0.00000
    168      11.4343      0.00000
    169      11.5583      0.00000
    170      11.5583      0.00000
    171      11.5583      0.00000
    172      11.5583      0.00000
    173      11.5686      0.00000
    174      11.5686      0.00000
    175      11.6574      0.00000
    176      11.6574      0.00000
    177      11.6574      0.00000
    178      11.6574      0.00000
    179      12.2297      0.00000
    180      12.2297      0.00000
    181      12.2297      0.00000
    182      12.2297      0.00000
    183      12.2771      0.00000
    184      12.2771      0.00000
    185      12.7494      0.00000
    186      12.7494      0.00000
    187      12.7494      0.00000
    188      12.7494      0.00000
    189      12.8639      0.00000
    190      12.8639      0.00000
    191      12.9742      0.00000
    192      12.9742      0.00000

 k-point     6 :       0.2500    0.0000    0.5000
  band No.  band energies     occupation 
      1     -23.8130      1.00000
      2     -23.8130      1.00000
      3     -23.6514      1.00000
      4     -23.6514      1.00000
      5     -23.6178      1.00000
      6     -23.6178      1.00000
      7     -14.5781      1.00000
      8     -14.5780      1.00000
      9     -14.5408      1.00000
     10     -14.5407      1.00000
     11     -14.5406      1.00000
     12     -14.5406      1.00000
     13     -14.4180      1.00000
     14     -14.4179      1.00000
     15     -14.4128      1.00000
     16     -14.4128      1.00000
     17     -14.3823      1.00000
     18     -14.3822      1.00000
     19     -14.3587      1.00000
     20     -14.3587      1.00000
     21     -14.3081      1.00000
     22     -14.3081      1.00000
     23     -14.3010      1.00000
     24     -14.3009      1.00000
     25     -14.2481      1.00000
     26     -14.2480      1.00000
     27     -14.2446      1.00000
     28     -14.2445      1.00000
     29     -14.2328      1.00000
     30     -14.2328      1.00000
     31     -14.2315      1.00000
     32     -14.2314      1.00000
     33     -14.2266      1.00000
     34     -14.2266      1.00000
     35     -14.2224      1.00000
     36     -14.2223      1.00000
     37     -11.6052      1.00000
     38     -11.6052      1.00000
     39     -11.4679      1.00000
     40     -11.4679      1.00000
     41     -11.1476      1.00000
     42     -11.1476      1.00000
     43     -10.9277      1.00000
     44     -10.9276      1.00000
     45     -10.7728      1.00000
     46     -10.7728      1.00000
     47     -10.7123      1.00000
     48     -10.7123      1.00000
     49      -7.9163      1.00000
     50      -7.9163      1.00000
     51      -7.8350      1.00000
     52      -7.8350      1.00000
     53      -7.4184      1.00000
     54      -7.4184      1.00000
     55      -7.3194      1.00000
     56      -7.3194      1.00000
     57      -7.2098      1.00000
     58      -7.2098      1.00000
     59      -7.1408      1.00000
     60      -7.1408      1.00000
     61      -6.8138      1.00000
     62      -6.8138      1.00000
     63      -6.7338      1.00000
     64      -6.7338      1.00000
     65      -6.6061      1.00000
     66      -6.6061      1.00000
     67      -6.4404      1.00000
     68      -6.4404      1.00000
     69      -5.9590      1.00000
     70      -5.9590      1.00000
     71      -5.8805      1.00000
     72      -5.8805      1.00000
     73      -0.4325      1.00000
     74      -0.4325      1.00000
     75      -0.2087      1.00000
     76      -0.2087      1.00000
     77       1.6115      1.00000
     78       1.6115      1.00000
     79       2.1276      1.00000
     80       2.1276      1.00000
     81       2.2654      1.00000
     82       2.2654      1.00000
     83       2.3691      1.00000
     84       2.3692      1.00000
     85       2.7846      1.00000
     86       2.7846      1.00000
     87       2.8476      1.00000
     88       2.8476      1.00000
     89       3.2238      1.00000
     90       3.2238      1.00000
     91       3.3603      1.00000
     92       3.3603      1.00000
     93       3.6096      1.00000
     94       3.6096      1.00000
     95       3.7189      1.00000
     96       3.7189      1.00000
     97       4.0668      1.00000
     98       4.0668      1.00000
     99       4.1681      1.00000
    100       4.1681      1.00000
    101       4.3215      1.00000
    102       4.3215      1.00000
    103       4.4303      1.00000
    104       4.4304      1.00000
    105       4.5800      1.00000
    106       4.5800      1.00000
    107       4.6571      1.00000
    108       4.6571      1.00000
    109       4.7274      1.00000
    110       4.7274      1.00000
    111       4.8778      1.00000
    112       4.8778      1.00000
    113       4.9534      1.00000
    114       4.9534      1.00000
    115       5.0934      1.00000
    116       5.0934      1.00000
    117       5.1633      1.00000
    118       5.1633      1.00000
    119       5.1854      1.00000
    120       5.1854      1.00000
    121       5.8611      1.00000
    122       5.8611      1.00000
    123       5.9920      1.00000
    124       5.9920      1.00000
    125       6.0218      1.00000
    126       6.0219      1.00000
    127       8.1601      0.00000
    128       8.1602      0.00000
    129       9.9490      0.00000
    130       9.9490      0.00000
    131      10.1825      0.00000
    132      10.1825      0.00000
    133      10.5321      0.00000
    134      10.5321      0.00000
    135      10.6612      0.00000
    136      10.6612      0.00000
    137      10.7628      0.00000
    138      10.7628      0.00000
    139      10.8143      0.00000
    140      10.8143      0.00000
    141      10.8886      0.00000
    142      10.8886      0.00000
    143      10.9002      0.00000
    144      10.9002      0.00000
    145      10.9339      0.00000
    146      10.9339      0.00000
    147      11.0443      0.00000
    148      11.0443      0.00000
    149      11.0573      0.00000
    150      11.0573      0.00000
    151      11.1098      0.00000
    152      11.1098      0.00000
    153      11.1355      0.00000
    154      11.1355      0.00000
    155      11.1572      0.00000
    156      11.1572      0.00000
    157      11.1764      0.00000
    158      11.1764      0.00000
    159      11.2077      0.00000
    160      11.2077      0.00000
    161      11.2421      0.00000
    162      11.2421      0.00000
    163      11.2724      0.00000
    164      11.2724      0.00000
    165      11.3100      0.00000
    166      11.3100      0.00000
    167      11.3730      0.00000
    168      11.3730      0.00000
    169      11.4201      0.00000
    170      11.4201      0.00000
    171      11.4797      0.00000
    172      11.4797      0.00000
    173      11.5223      0.00000
    174      11.5223      0.00000
    175      11.7460      0.00000
    176      11.7460      0.00000
    177      11.8353      0.00000
    178      11.8353      0.00000
    179      11.9276      0.00000
    180      11.9276      0.00000
    181      12.1505      0.00000
    182      12.1505      0.00000
    183      12.2505      0.00000
    184      12.2505      0.00000
    185      12.3734      0.00000
    186      12.3734      0.00000
    187      12.5268      0.00000
    188      12.5268      0.00000
    189      12.5872      0.00000
    190      12.5872      0.00000
    191      12.8462      0.00000
    192      12.8462      0.00000

 k-point     7 :       0.5000    0.0000    0.5000
  band No.  band energies     occupation 
      1     -23.7380      1.00000
      2     -23.7380      1.00000
      3     -23.7118      1.00000
      4     -23.7118      1.00000
      5     -23.6351      1.00000
      6     -23.6351      1.00000
      7     -14.5608      1.00000
      8     -14.5608      1.00000
      9     -14.5598      1.00000
     10     -14.5597      1.00000
     11     -14.5514      1.00000
     12     -14.5513      1.00000
     13     -14.4180      1.00000
     14     -14.4180      1.00000
     15     -14.4002      1.00000
     16     -14.4001      1.00000
     17     -14.3999      1.00000
     18     -14.3999      1.00000
     19     -14.3281      1.00000
     20     -14.3280      1.00000
     21     -14.3266      1.00000
     22     -14.3265      1.00000
     23     -14.2938      1.00000
     24     -14.2937      1.00000
     25     -14.2500      1.00000
     26     -14.2500      1.00000
     27     -14.2385      1.00000
     28     -14.2384      1.00000
     29     -14.2375      1.00000
     30     -14.2374      1.00000
     31     -14.2277      1.00000
     32     -14.2276      1.00000
     33     -14.2270      1.00000
     34     -14.2269      1.00000
     35     -14.2253      1.00000
     36     -14.2252      1.00000
     37     -11.7282      1.00000
     38     -11.7282      1.00000
     39     -11.2654      1.00000
     40     -11.2654      1.00000
     41     -11.2093      1.00000
     42     -11.2092      1.00000
     43     -10.9279      1.00000
     44     -10.9279      1.00000
     45     -10.8999      1.00000
     46     -10.8999      1.00000
     47     -10.4999      1.00000
     48     -10.4999      1.00000
     49      -7.9237      1.00000
     50      -7.9237      1.00000
     51      -7.6790      1.00000
     52      -7.6789      1.00000
     53      -7.6741      1.00000
     54      -7.6741      1.00000
     55      -7.6032      1.00000
     56      -7.6032      1.00000
     57      -7.0982      1.00000
     58      -7.0982      1.00000
     59      -7.0549      1.00000
     60      -7.0549      1.00000
     61      -6.9209      1.00000
     62      -6.9209      1.00000
     63      -6.8466      1.00000
     64      -6.8466      1.00000
     65      -6.3869      1.00000
     66      -6.3869      1.00000
     67      -6.1984      1.00000
     68      -6.1984      1.00000
     69      -6.1052      1.00000
     70      -6.1052      1.00000
     71      -5.9570      1.00000
     72      -5.9569      1.00000
     73      -0.5071      1.00000
     74      -0.5071      1.00000
     75       0.5747      1.00000
     76       0.5747      1.00000
     77       0.5751      1.00000
     78       0.5751      1.00000
     79       2.3668      1.00000
     80       2.3669      1.00000
     81       2.3831      1.00000
     82       2.3831      1.00000
     83       2.6162      1.00000
     84       2.6162      1.00000
     85       2.9791      1.00000
     86       2.9791      1.00000
     87       3.0928      1.00000
     88       3.0928      1.00000
     89       3.1005      1.00000
     90       3.1005      1.00000
     91       3.4353      1.00000
     92       3.4353      1.00000
     93       3.8398      1.00000
     94       3.8398      1.00000
     95       3.9002      1.00000
     96       3.9002      1.00000
     97       3.9504      1.00000
     98       3.9504      1.00000
     99       3.9953      1.00000
    100       3.9953      1.00000
    101       4.1102      1.00000
    102       4.1103      1.00000
    103       4.1287      1.00000
    104       4.1287      1.00000
    105       4.2575      1.00000
    106       4.2575      1.00000
    107       4.4710      1.00000
    108       4.4710      1.00000
    109       4.8151      1.00000
    110       4.8151      1.00000
    111       4.8232      1.00000
    112       4.8232      1.00000
    113       4.9678      1.00000
    114       4.9679      1.00000
    115       5.1517      1.00000
    116       5.1517      1.00000
    117       5.1578      1.00000
    118       5.1578      1.00000
    119       5.2724      1.00000
    120       5.2724      1.00000
    121       5.9110      1.00000
    122       5.9110      1.00000
    123       6.0006      1.00000
    124       6.0006      1.00000
    125       6.0415      1.00000
    126       6.0415      1.00000
    127       9.2595      0.00000
    128       9.2595      0.00000
    129       9.3172      0.00000
    130       9.3172      0.00000
    131       9.3381      0.00000
    132       9.3381      0.00000
    133      10.3682      0.00000
    134      10.3682      0.00000
    135      10.3994      0.00000
    136      10.3994      0.00000
    137      10.7940      0.00000
    138      10.7940      0.00000
    139      10.8684      0.00000
    140      10.8684      0.00000
    141      10.8700      0.00000
    142      10.8700      0.00000
    143      10.9452      0.00000
    144      10.9452      0.00000
    145      10.9623      0.00000
    146      10.9623      0.00000
    147      11.0482      0.00000
    148      11.0482      0.00000
    149      11.0526      0.00000
    150      11.0526      0.00000
    151      11.0750      0.00000
    152      11.0750      0.00000
    153      11.0965      0.00000
    154      11.0965      0.00000
    155      11.1322      0.00000
    156      11.1322      0.00000
    157      11.1366      0.00000
    158      11.1366      0.00000
    159      11.1580      0.00000
    160      11.1580      0.00000
    161      11.1748      0.00000
    162      11.1748      0.00000
    163      11.2605      0.00000
    164      11.2605      0.00000
    165      11.2815      0.00000
    166      11.2815      0.00000
    167      11.3198      0.00000
    168      11.3198      0.00000
    169      11.4440      0.00000
    170      11.4440      0.00000
    171      11.4823      0.00000
    172      11.4823      0.00000
    173      11.4946      0.00000
    174      11.4946      0.00000
    175      11.5314      0.00000
    176      11.5314      0.00000
    177      11.7806      0.00000
    178      11.7806      0.00000
    179      11.7918      0.00000
    180      11.7918      0.00000
    181      12.1491      0.00000
    182      12.1491      0.00000
    183      12.2108      0.00000
    184      12.2108      0.00000
    185      12.3099      0.00000
    186      12.3099      0.00000
    187      12.4750      0.00000
    188      12.4750      0.00000
    189      12.5965      0.00000
    190      12.5965      0.00000
    191      12.9537      0.00000
    192      12.9537      0.00000

 k-point     8 :       0.2500    0.2500    0.5000
  band No.  band energies     occupation 
      1     -23.7526      1.00000
      2     -23.7526      1.00000
      3     -23.6839      1.00000
      4     -23.6839      1.00000
      5     -23.6484      1.00000
      6     -23.6484      1.00000
      7     -14.5666      1.00000
      8     -14.5666      1.00000
      9     -14.5528      1.00000
     10     -14.5528      1.00000
     11     -14.5527      1.00000
     12     -14.5526      1.00000
     13     -14.4127      1.00000
     14     -14.4126      1.00000
     15     -14.4105      1.00000
     16     -14.4104      1.00000
     17     -14.3960      1.00000
     18     -14.3959      1.00000
     19     -14.3340      1.00000
     20     -14.3339      1.00000
     21     -14.3070      1.00000
     22     -14.3069      1.00000
     23     -14.3062      1.00000
     24     -14.3061      1.00000
     25     -14.2463      1.00000
     26     -14.2462      1.00000
     27     -14.2448      1.00000
     28     -14.2447      1.00000
     29     -14.2357      1.00000
     30     -14.2356      1.00000
     31     -14.2280      1.00000
     32     -14.2280      1.00000
     33     -14.2271      1.00000
     34     -14.2271      1.00000
     35     -14.2240      1.00000
     36     -14.2239      1.00000
     37     -11.5663      1.00000
     38     -11.5663      1.00000
     39     -11.5440      1.00000
     40     -11.5440      1.00000
     41     -11.2163      1.00000
     42     -11.2163      1.00000
     43     -10.8521      1.00000
     44     -10.8520      1.00000
     45     -10.6893      1.00000
     46     -10.6893      1.00000
     47     -10.6750      1.00000
     48     -10.6749      1.00000
     49      -7.8221      1.00000
     50      -7.8221      1.00000
     51      -7.7757      1.00000
     52      -7.7757      1.00000
     53      -7.4648      1.00000
     54      -7.4648      1.00000
     55      -7.4440      1.00000
     56      -7.4440      1.00000
     57      -7.4040      1.00000
     58      -7.4040      1.00000
     59      -7.3547      1.00000
     60      -7.3547      1.00000
     61      -6.8053      1.00000
     62      -6.8053      1.00000
     63      -6.6033      1.00000
     64      -6.6033      1.00000
     65      -6.5593      1.00000
     66      -6.5593      1.00000
     67      -6.2057      1.00000
     68      -6.2057      1.00000
     69      -6.0317      1.00000
     70      -6.0316      1.00000
     71      -5.9676      1.00000
     72      -5.9676      1.00000
     73      -0.1209      1.00000
     74      -0.1209      1.00000
     75      -0.1113      1.00000
     76      -0.1113      1.00000
     77       0.8142      1.00000
     78       0.8143      1.00000
     79       2.2812      1.00000
     80       2.2812      1.00000
     81       2.5190      1.00000
     82       2.5191      1.00000
     83       2.5971      1.00000
     84       2.5972      1.00000
     85       2.9671      1.00000
     86       2.9671      1.00000
     87       3.0183      1.00000
     88       3.0183      1.00000
     89       3.2860      1.00000
     90       3.2860      1.00000
     91       3.5839      1.00000
     92       3.5839      1.00000
     93       3.5981      1.00000
     94       3.5981      1.00000
     95       3.7050      1.00000
     96       3.7050      1.00000
     97       3.9772      1.00000
     98       3.9772      1.00000
     99       4.0546      1.00000
    100       4.0546      1.00000
    101       4.0874      1.00000
    102       4.0875      1.00000
    103       4.2092      1.00000
    104       4.2092      1.00000
    105       4.3875      1.00000
    106       4.3875      1.00000
    107       4.4349      1.00000
    108       4.4349      1.00000
    109       4.9146      1.00000
    110       4.9146      1.00000
    111       4.9492      1.00000
    112       4.9492      1.00000
    113       4.9697      1.00000
    114       4.9697      1.00000
    115       5.0809      1.00000
    116       5.0809      1.00000
    117       5.1395      1.00000
    118       5.1395      1.00000
    119       5.1930      1.00000
    120       5.1930      1.00000
    121       5.9116      1.00000
    122       5.9116      1.00000
    123       6.0031      1.00000
    124       6.0031      1.00000
    125       6.0133      1.00000
    126       6.0134      1.00000
    127       8.8592      0.00000
    128       8.8592      0.00000
    129       9.4084      0.00000
    130       9.4084      0.00000
    131       9.4622      0.00000
    132       9.4622      0.00000
    133      10.4008      0.00000
    134      10.4008      0.00000
    135      10.5452      0.00000
    136      10.5452      0.00000
    137      10.5659      0.00000
    138      10.5659      0.00000
    139      10.8934      0.00000
    140      10.8934      0.00000
    141      10.8981      0.00000
    142      10.8981      0.00000
    143      10.9637      0.00000
    144      10.9637      0.00000
    145      10.9957      0.00000
    146      10.9957      0.00000
    147      11.0173      0.00000
    148      11.0173      0.00000
    149      11.0649      0.00000
    150      11.0649      0.00000
    151      11.1129      0.00000
    152      11.1129      0.00000
    153      11.1155      0.00000
    154      11.1155      0.00000
    155      11.1326      0.00000
    156      11.1326      0.00000
    157      11.1452      0.00000
    158      11.1452      0.00000
    159      11.1618      0.00000
    160      11.1618      0.00000
    161      11.1688      0.00000
    162      11.1688      0.00000
    163      11.2760      0.00000
    164      11.2760      0.00000
    165      11.2832      0.00000
    166      11.2832      0.00000
    167      11.3620      0.00000
    168      11.3620      0.00000
    169      11.4642      0.00000
    170      11.4642      0.00000
    171      11.4677      0.00000
    172      11.4677      0.00000
    173      11.4969      0.00000
    174      11.4969      0.00000
    175      11.6203      0.00000
    176      11.6204      0.00000
    177      11.6864      0.00000
    178      11.6864      0.00000
    179      11.7211      0.00000
    180      11.7211      0.00000
    181      12.1570      0.00000
    182      12.1570      0.00000
    183      12.2265      0.00000
    184      12.2265      0.00000
    185      12.2277      0.00000
    186      12.2277      0.00000
    187      12.6268      0.00000
    188      12.6268      0.00000
    189      12.7129      0.00000
    190      12.7129      0.00000
    191      13.0343      0.00000
    192      13.0344      0.00000


--------------------------------------------------------------------------------------------------------


 soft charge-density along one line, spin component           1
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 soft charge-density along one line, spin component           2
         0         1         2         3         4         5         6         7         8         9
 total charge-density along one line
 
 pseudopotential strength for first ion, spin component:           1
 52.368  -3.890  -0.000  -0.003  -0.000  -0.000  -0.001  -0.000
 -3.890  -4.509   0.000   0.001   0.000  -0.000  -0.001  -0.000
 -0.000   0.000  18.688  -0.000  -0.000   0.794  -0.000  -0.000
 -0.003   0.001  -0.000  18.683  -0.000  -0.000   0.794  -0.000
 -0.000   0.000  -0.000  -0.000  18.688  -0.000  -0.000   0.794
 -0.000  -0.000   0.794  -0.000  -0.000   6.924  -0.000  -0.000
 -0.001  -0.001  -0.000   0.794  -0.000  -0.000   6.923  -0.000
 -0.000  -0.000  -0.000  -0.000   0.794  -0.000  -0.000   6.924
 -0.000   0.000  -0.001   0.000  -0.000  -0.001   0.000  -0.000
 -0.000   0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000
 -0.002  -0.000   0.000  -0.003   0.000  -0.000  -0.001  -0.000
 -0.000   0.000   0.000  -0.000   0.000   0.000  -0.000  -0.000
  0.000  -0.000   0.000  -0.000   0.001   0.000  -0.000   0.001
  0.000   0.000   0.000  -0.000   0.000   0.000  -0.000   0.000
  0.000   0.000   0.000   0.000  -0.000  -0.000   0.000  -0.000
 -0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
  0.000   0.000  -0.000   0.000   0.000  -0.000   0.000  -0.000
 -0.000  -0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000
  0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000
 -0.000  -0.000   0.009   0.000   0.000   0.001  -0.000  -0.000
  0.000   0.000   0.003   0.000  -0.000  -0.000  -0.000   0.000
 -0.001   0.000  -0.000   0.002  -0.000   0.000  -0.000   0.000
  0.000   0.000  -0.000   0.000   0.003   0.000  -0.000  -0.000
  0.000   0.000  -0.000  -0.000  -0.009   0.000   0.000  -0.001
 -0.000  -0.001  -0.000  -0.008  -0.000   0.000  -0.001   0.000
 -0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000  -0.000
  0.000   0.000  -0.008  -0.000  -0.000  -0.002   0.000  -0.000
 -0.000  -0.000  -0.002  -0.000   0.000  -0.000   0.000   0.000
  0.000  -0.000   0.000  -0.002   0.000   0.000   0.000   0.000
 -0.000  -0.000   0.000  -0.000  -0.002   0.000   0.000  -0.000
 -0.000  -0.000   0.000   0.000   0.008   0.000  -0.000   0.002
  0.000   0.000   0.000   0.006   0.000   0.000   0.002   0.000
 pseudopotential strength for first ion, spin component:           2
 52.368  -3.890  -0.000  -0.003  -0.000  -0.000  -0.001  -0.000
 -3.890  -4.509   0.000   0.001   0.000  -0.000  -0.001  -0.000
 -0.000   0.000  18.688  -0.000  -0.000   0.794  -0.000  -0.000
 -0.003   0.001  -0.000  18.683  -0.000  -0.000   0.794  -0.000
 -0.000   0.000  -0.000  -0.000  18.688  -0.000  -0.000   0.794
 -0.000  -0.000   0.794  -0.000  -0.000   6.924  -0.000  -0.000
 -0.001  -0.001  -0.000   0.794  -0.000  -0.000   6.923  -0.000
 -0.000  -0.000  -0.000  -0.000   0.794  -0.000  -0.000   6.924
 -0.000   0.000  -0.001   0.000  -0.000  -0.001   0.000  -0.000
 -0.000   0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000
 -0.002  -0.000   0.000  -0.003   0.000  -0.000  -0.001  -0.000
 -0.000   0.000   0.000  -0.000   0.000   0.000  -0.000  -0.000
  0.000  -0.000   0.000  -0.000   0.001   0.000  -0.000   0.001
  0.000   0.000   0.000  -0.000   0.000   0.000  -0.000   0.000
  0.000   0.000   0.000   0.000  -0.000  -0.000   0.000  -0.000
 -0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000
  0.000   0.000  -0.000   0.000   0.000  -0.000   0.000  -0.000
 -0.000  -0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000
  0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000
 -0.000  -0.000   0.009   0.000   0.000   0.001  -0.000  -0.000
  0.000   0.000   0.003   0.000  -0.000  -0.000  -0.000   0.000
 -0.001   0.000  -0.000   0.002  -0.000   0.000  -0.000   0.000
  0.000   0.000  -0.000   0.000   0.003   0.000  -0.000  -0.000
  0.000   0.000  -0.000  -0.000  -0.009   0.000   0.000  -0.001
 -0.000  -0.001  -0.000  -0.008  -0.000   0.000  -0.001   0.000
 -0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000  -0.000
  0.000   0.000  -0.008  -0.000  -0.000  -0.002   0.000  -0.000
 -0.000  -0.000  -0.002  -0.000   0.000  -0.000   0.000   0.000
  0.000  -0.000   0.000  -0.002   0.000   0.000   0.000   0.000
 -0.000  -0.000   0.000  -0.000  -0.002   0.000   0.000  -0.000
 -0.000  -0.000   0.000   0.000   0.008   0.000  -0.000   0.002
  0.000   0.000   0.000   0.006   0.000   0.000   0.002   0.000
 total augmentation occupancy for first ion, spin component:           1
  2.288  -0.544   0.000  -0.004   0.000   0.000   0.001   0.000  -0.000   0.000   0.019   0.000  -0.000   0.000   0.000  -0.039
 -0.544   1.001   0.000   0.006   0.000   0.000  -0.001   0.000   0.000   0.000  -0.030   0.000   0.000  -0.000   0.000   0.072
  0.000   0.000   2.053   0.000   0.000  -0.023   0.000  -0.000  -0.020   0.021   0.000   0.000  -0.000   0.021  -0.025  -0.000
 -0.004   0.006   0.000   2.046   0.000   0.000  -0.018   0.000   0.000  -0.000  -0.028   0.000   0.000  -0.000  -0.000   0.021
  0.000  -0.000   0.000   0.000   2.053  -0.000  -0.000  -0.023   0.000   0.000  -0.000   0.021   0.020  -0.000  -0.000   0.000
  0.000   0.000  -0.023   0.000  -0.000   0.025   0.000   0.000   0.007  -0.008  -0.000  -0.000   0.000  -0.008   0.011   0.000
  0.001  -0.001   0.000  -0.018   0.000   0.000   0.021   0.000  -0.000  -0.000   0.011   0.000  -0.000   0.000   0.000  -0.011
  0.000  -0.000   0.000   0.000  -0.023   0.000   0.000   0.025  -0.000  -0.000   0.000  -0.008  -0.007   0.000   0.000  -0.000
 -0.000   0.000  -0.020   0.000   0.000   0.007  -0.000  -0.000   0.170   0.076  -0.000  -0.000   0.000  -0.232  -0.113   0.000
  0.000   0.000   0.021  -0.000   0.000  -0.008  -0.000  -0.000   0.076   0.164  -0.000   0.000   0.000  -0.113  -0.221   0.000
  0.019  -0.030   0.000  -0.028  -0.000  -0.000   0.011   0.000  -0.000  -0.000   0.152   0.000  -0.000   0.000   0.000  -0.182
  0.000   0.000  -0.000   0.000   0.021  -0.000   0.000  -0.008  -0.000   0.000   0.000   0.164  -0.076   0.000  -0.000  -0.000
 -0.000   0.000  -0.000   0.000   0.020   0.000  -0.000  -0.007   0.000   0.000  -0.000  -0.076   0.170  -0.000  -0.000   0.000
  0.000  -0.000   0.021  -0.000  -0.000  -0.008   0.000   0.000  -0.232  -0.113   0.000   0.000  -0.000   0.345   0.189  -0.000
  0.000   0.000  -0.025  -0.000  -0.000   0.011   0.000   0.000  -0.113  -0.221   0.000   0.000  -0.000   0.189   0.327  -0.000
 -0.039   0.072  -0.000   0.021   0.000   0.000  -0.011  -0.000   0.000   0.000  -0.182  -0.000   0.000  -0.000  -0.000   0.239
 -0.000  -0.000  -0.000  -0.000  -0.025  -0.000   0.000   0.011   0.000   0.000  -0.000  -0.221   0.113  -0.000   0.000   0.000
  0.000  -0.000   0.000  -0.000  -0.021  -0.000   0.000   0.008  -0.000  -0.000   0.000   0.113  -0.232   0.000   0.000  -0.000
  0.000  -0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000
  0.000  -0.000  -0.089  -0.000  -0.000   0.037   0.000   0.000  -0.004  -0.020  -0.000  -0.000   0.000   0.011   0.030   0.000
 -0.000   0.000  -0.001  -0.000  -0.000  -0.002   0.000   0.000  -0.021  -0.033  -0.000  -0.000  -0.000   0.032   0.042   0.000
 -0.041   0.063   0.000   0.027  -0.000  -0.000  -0.014   0.000  -0.000   0.000   0.069  -0.000  -0.000   0.000  -0.000  -0.076
 -0.000   0.000   0.000  -0.000  -0.001  -0.000   0.000  -0.002   0.000  -0.000  -0.000  -0.033   0.021  -0.000   0.000   0.000
  0.000  -0.000  -0.000  -0.000   0.089   0.000   0.000  -0.037  -0.000  -0.000  -0.000   0.020  -0.004   0.000   0.000   0.000
  0.020  -0.028  -0.000   0.081   0.000   0.000  -0.034   0.000  -0.000   0.000   0.032   0.000  -0.000  -0.000   0.000  -0.044
  0.000  -0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000
  0.000  -0.000  -0.063  -0.000  -0.000   0.029   0.000   0.000  -0.005  -0.017  -0.000  -0.000  -0.000   0.010   0.026   0.000
 -0.000   0.000  -0.002  -0.000  -0.000  -0.001   0.000   0.000  -0.017  -0.025  -0.000   0.000  -0.000   0.027   0.034   0.000
 -0.032   0.048   0.000   0.019  -0.000  -0.000  -0.012   0.000  -0.000   0.000   0.050  -0.000  -0.000   0.000  -0.000  -0.057
 -0.000   0.000   0.000  -0.000  -0.002  -0.000   0.000  -0.001   0.000  -0.000  -0.000  -0.025   0.017  -0.000   0.000   0.000
  0.000  -0.000  -0.000  -0.000   0.063   0.000   0.000  -0.029  -0.000  -0.000  -0.000   0.017  -0.005   0.000   0.000   0.000
  0.017  -0.022  -0.000   0.056   0.000   0.000  -0.026   0.000  -0.000  -0.000   0.026   0.000  -0.000   0.000  -0.000  -0.036
 total augmentation occupancy for first ion, spin component:           2
 -0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000   0.000   0.000
  0.000  -0.000  -0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000   0.000   0.000
  0.000   0.000  -0.000   0.000  -0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000   0.000
  0.000  -0.000  -0.000  -0.000   0.000  -0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000  -0.000
  0.000   0.000  -0.000  -0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000
  0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000
 -0.000   0.000  -0.000   0.000   0.000  -0.000  -0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000  -0.000   0.000
  0.000   0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000   0.000
  0.000  -0.000  -0.000  -0.000   0.000   0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000   0.000   0.000
 -0.000  -0.000  -0.000  -0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000   0.000
 -0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000
  0.000   0.000  -0.000   0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000  -0.000  -0.000  -0.000
  0.000  -0.000  -0.000   0.000   0.000   0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000   0.000   0.000
 -0.000  -0.000  -0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000   0.000   0.000  -0.000
  0.000  -0.000  -0.000  -0.000  -0.000   0.000  -0.000   0.000   0.000  -0.000   0.000  -0.000   0.000   0.000   0.000   0.000
  0.000   0.000   0.000  -0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000  -0.000   0.000  -0.000
 -0.000   0.000  -0.000   0.000  -0.000   0.000  -0.000   0.000   0.000  -0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000
 -0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000   0.000   0.000  -0.000
  0.000  -0.000   0.000   0.000  -0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000  -0.000   0.000  -0.000   0.000
 -0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000
 -0.000   0.000  -0.000  -0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000   0.000  -0.000
  0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000   0.000   0.000
 -0.000   0.000   0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000
 -0.000   0.000  -0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000   0.000   0.000  -0.000
  0.000  -0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000
 -0.000  -0.000   0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000
 -0.000   0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000
  0.000   0.000  -0.000   0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000  -0.000   0.000  -0.000   0.000   0.000
  0.000  -0.000   0.000  -0.000  -0.000   0.000  -0.000  -0.000  -0.000  -0.000   0.000   0.000  -0.000  -0.000   0.000   0.000
  0.000   0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000  -0.000   0.000  -0.000  -0.000   0.000   0.000  -0.000
 -0.000   0.000  -0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000   0.000   0.000  -0.000
 -0.000  -0.000  -0.000   0.000   0.000   0.000   0.000   0.000  -0.000   0.000   0.000   0.000   0.000  -0.000  -0.000   0.000


------------------------ aborting loop because EDIFF is reached ----------------------------------------


     CALCP:  cpu time   65.8052: real time   65.7551
 
            Ionic dipole moment: p[ion]=(     0.00000     0.00000    29.75988 ) electrons Angst
 
    Spin resolved dipole moment: p[sp1]=(     0.00004     0.00000    -4.29028 ) electrons Angst
                                 p[sp2]=(     0.00004     0.00000    -4.29027 ) electrons Angst
 
 Total electronic dipole moment: p[elc]=(     0.00008     0.00000    -8.58055 ) electrons Angst
 
 
 


 total charge     
 
# of ion       s       p       d       f       tot
--------------------------------------------------
    1        2.021   5.716   0.624   0.249   8.611
    2        2.021   5.716   0.624   0.249   8.611
    3        2.020   5.715   0.617   0.253   8.605
    4        2.020   5.715   0.617   0.253   8.605
    5        2.020   5.715   0.617   0.253   8.605
    6        2.020   5.715   0.617   0.253   8.605
    7        1.243   2.388   0.000   0.000   3.631
    8        1.243   2.388   0.000   0.000   3.631
    9        1.244   2.383   0.000   0.000   3.627
   10        1.244   2.383   0.000   0.000   3.627
   11        1.244   2.383   0.000   0.000   3.627
   12        1.244   2.383   0.000   0.000   3.627
   13        1.570   3.463   0.000   0.000   5.033
   14        1.570   3.457   0.000   0.000   5.028
   15        1.571   3.465   0.000   0.000   5.036
   16        1.570   3.463   0.000   0.000   5.033
   17        1.570   3.457   0.000   0.000   5.028
   18        1.571   3.465   0.000   0.000   5.036
   19        1.569   3.452   0.000   0.000   5.021
   20        1.569   3.459   0.000   0.000   5.028
   21        1.569   3.449   0.000   0.000   5.018
   22        1.569   3.452   0.000   0.000   5.021
   23        1.569   3.459   0.000   0.000   5.028
   24        1.569   3.449   0.000   0.000   5.018
   25        0.964   1.021  10.219   0.000  12.205
   26        0.965   1.021  10.219   0.000  12.204
   27        0.964   1.021  10.219   0.000  12.205
   28        0.964   1.021  10.219   0.000  12.205
   29        0.965   1.021  10.219   0.000  12.204
   30        0.964   1.021  10.219   0.000  12.205
--------------------------------------------------
tot         44.206  96.217  65.032   1.511 206.966

 


 magnetization (x)
 
# of ion       s       p       d       f       tot
--------------------------------------------------
    1        0.000  -0.000   0.000   0.000   0.000
    2        0.000  -0.000   0.000   0.000   0.000
    3        0.000  -0.000   0.000   0.000   0.000
    4        0.000  -0.000   0.000   0.000   0.000
    5        0.000  -0.000   0.000   0.000   0.000
    6        0.000  -0.000   0.000   0.000   0.000
    7        0.000  -0.000   0.000   0.000  -0.000
    8        0.000  -0.000   0.000   0.000  -0.000
    9        0.000  -0.000   0.000   0.000  -0.000
   10        0.000  -0.000   0.000   0.000  -0.000
   11        0.000  -0.000   0.000   0.000  -0.000
   12        0.000  -0.000   0.000   0.000  -0.000
   13       -0.000  -0.000   0.000   0.000  -0.000
   14       -0.000  -0.000   0.000   0.000  -0.000
   15       -0.000  -0.000   0.000   0.000  -0.000
   16       -0.000  -0.000   0.000   0.000  -0.000
   17       -0.000  -0.000   0.000   0.000  -0.000
   18       -0.000  -0.000   0.000   0.000  -0.000
   19       -0.000  -0.000   0.000   0.000  -0.000
   20       -0.000  -0.000   0.000   0.000  -0.000
   21       -0.000  -0.000   0.000   0.000  -0.000
   22       -0.000  -0.000   0.000   0.000  -0.000
   23       -0.000  -0.000   0.000   0.000  -0.000
   24       -0.000  -0.000   0.000   0.000  -0.000
   25        0.000   0.000   0.000   0.000   0.000
   26        0.000   0.000   0.000   0.000   0.000
   27        0.000   0.000   0.000   0.000   0.000
   28        0.000   0.000   0.000   0.000   0.000
   29        0.000   0.000   0.000   0.000   0.000
   30        0.000   0.000   0.000   0.000   0.000
--------------------------------------------------
tot          0.000  -0.000   0.000   0.000  -0.000

 
    CHARGE:  cpu time    0.5783: real time    0.5778
    FORLOC:  cpu time    0.0250: real time    0.0260
    FORNL :  cpu time    2.9667: real time    2.9654
    STRESS:  cpu time    8.9153: real time    8.9094
    FORCOR:  cpu time    0.1897: real time    0.1900
    FORHAR:  cpu time    0.0387: real time    0.0389
    MIXING:  cpu time    0.0061: real time    0.0061
    OFIELD:  cpu time    0.0004: real time    0.0009

  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  --------------------------------------------------------------------------------------
  Alpha Z   865.42119   865.42119   865.42119
  Ewald   -6396.02499 -6396.02468 -5370.45366    -0.00002    -0.00000    -0.00000
  Hartree  1634.86736  1634.86738  2379.31193     0.00002     0.00003     0.00002
  E(xc)   -1153.01968 -1153.01972 -1152.59537    -0.00003     0.00000     0.00000
  Local    -581.22210  -581.22248 -2366.55382    -0.00039    -0.00006    -0.00004
  n-local   925.55459   925.11772   935.47573     0.98041    -0.00000    -0.00000
  augment   923.28285   923.28313   926.81202     0.00030     0.00001     0.00001
  Kinetic  3781.95199  3781.59103  3781.70594    -0.99625     0.00000    -0.00000
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  -------------------------------------------------------------------------------------
  Total       0.41238     0.41238    -0.87605     0.00000     0.00000     0.00000
  in kB       1.40512     1.40512    -2.98502     0.00000     0.00000     0.00000
  external pressure =       -0.06 kB  Pullay stress =        0.00 kB


 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
  energy-cutoff  :      520.00
  volume of cell :      470.21
      direct lattice vectors                 reciprocal lattice vectors
     6.506500000  0.000000000  0.000000000     0.153692461  0.088734389  0.000000000
    -3.253250000  5.634794000  0.000000000     0.000000000  0.177468777  0.000000000
     0.000000000  0.000000000 12.825300000     0.000000000  0.000000000  0.077970886

  length of vectors
     6.506500000  6.506499749 12.825300000     0.177468770  0.177468777  0.077970886


 FORCES acting on ions
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction
 -----------------------------------------------------------------------------------------------
   0.349E-04 0.554E-04 0.113E+02   -.336E-05 -.213E-13 -.109E+02   -.322E-18 0.148E-17 -.217E+00   -.538E-08 -.937E-08 0.115E-03
   -.420E-03 -.722E-03 0.113E+02   0.336E-05 -.178E-13 -.109E+02   0.716E-19 -.442E-18 -.217E+00   -.157E-07 -.271E-07 0.115E-03
   0.368E-02 -.208E-02 -.351E+02   0.151E-02 -.869E-03 0.344E+02   -.145E-18 -.378E-18 0.558E+00   -.186E-04 0.107E-04 -.101E-02
   -.433E-02 0.151E-02 -.351E+02   -.151E-02 0.869E-03 0.344E+02   -.165E-17 0.776E-18 0.558E+00   0.186E-04 -.107E-04 -.101E-02
   0.348E-02 -.299E-02 -.351E+02   0.150E-02 -.869E-03 0.344E+02   -.493E-18 0.239E-18 0.558E+00   -.186E-04 0.107E-04 -.101E-02
   -.364E-02 0.215E-02 -.351E+02   -.150E-02 0.869E-03 0.344E+02   -.156E-17 0.512E-18 0.558E+00   0.186E-04 -.108E-04 -.101E-02
   -.880E-04 -.181E-03 0.565E+01   -.229E-04 0.000E+00 -.103E+02   0.612E-17 0.131E-17 0.452E+01   -.439E-07 -.749E-07 0.851E-03
   -.147E-03 -.225E-03 0.565E+01   0.229E-04 0.444E-15 -.103E+02   -.968E-18 0.156E-17 0.452E+01   0.243E-07 0.410E-07 0.851E-03
   0.171E-02 -.125E-02 0.201E+01   0.103E-02 -.583E-03 0.308E+01   0.436E-17 -.101E-17 -.502E+01   -.960E-05 0.548E-05 -.283E-03
   -.195E-02 0.838E-03 0.201E+01   -.103E-02 0.583E-03 0.308E+01   -.156E-17 -.506E-19 -.502E+01   0.958E-05 -.552E-05 -.283E-03
   0.172E-02 -.128E-02 0.201E+01   0.987E-03 -.583E-03 0.308E+01   0.311E-17 -.102E-17 -.502E+01   -.957E-05 0.553E-05 -.283E-03
   -.196E-02 0.868E-03 0.201E+01   -.987E-03 0.583E-03 0.308E+01   -.482E-17 0.394E-18 -.502E+01   0.955E-05 -.557E-05 -.283E-03
   0.242E+02 -.221E-03 -.754E+02   -.261E+02 0.638E-03 0.804E+02   0.185E+01 0.216E-18 -.499E+01   -.169E-03 -.847E-05 -.212E-02
   -.121E+02 0.210E+02 -.754E+02   0.130E+02 -.226E+02 0.804E+02   -.924E+00 0.160E+01 -.499E+01   0.773E-04 -.151E-03 -.212E-02
   -.121E+02 -.210E+02 -.754E+02   0.130E+02 0.226E+02 0.804E+02   -.924E+00 -.160E+01 -.499E+01   0.919E-04 0.159E-03 -.213E-02
   -.242E+02 -.451E-03 -.754E+02   0.261E+02 -.638E-03 0.804E+02   -.185E+01 0.371E-18 -.499E+01   0.169E-03 0.837E-05 -.212E-02
   0.121E+02 -.210E+02 -.754E+02   -.130E+02 0.226E+02 0.804E+02   0.924E+00 -.160E+01 -.499E+01   -.773E-04 0.151E-03 -.212E-02
   0.121E+02 0.210E+02 -.754E+02   -.130E+02 -.226E+02 0.804E+02   0.924E+00 0.160E+01 -.499E+01   -.920E-04 -.159E-03 -.213E-02
   0.234E+02 -.423E-03 0.692E+02   -.251E+02 -.461E-03 -.750E+02   0.168E+01 -.137E-18 0.576E+01   0.408E-04 0.692E-05 0.148E-02
   -.117E+02 0.203E+02 0.692E+02   0.125E+02 -.217E+02 -.750E+02   -.838E+00 0.145E+01 0.576E+01   -.144E-04 0.388E-04 0.148E-02
   -.117E+02 -.203E+02 0.692E+02   0.125E+02 0.217E+02 -.750E+02   -.838E+00 -.145E+01 0.576E+01   -.264E-04 -.458E-04 0.148E-02
   -.234E+02 -.338E-03 0.692E+02   0.251E+02 0.461E-03 -.750E+02   -.168E+01 -.221E-18 0.576E+01   -.408E-04 -.703E-05 0.148E-02
   0.117E+02 -.203E+02 0.692E+02   -.125E+02 0.217E+02 -.750E+02   0.838E+00 -.145E+01 0.576E+01   0.143E-04 -.389E-04 0.148E-02
   0.117E+02 0.203E+02 0.692E+02   -.125E+02 -.217E+02 -.750E+02   0.838E+00 0.145E+01 0.576E+01   0.264E-04 0.457E-04 0.148E-02
   -.901E+00 0.351E-03 0.233E+02   0.109E+01 0.106E-02 -.234E+02   -.169E+00 0.301E-17 0.144E+00   -.251E-03 -.149E-04 0.110E-02
   0.451E+00 -.780E+00 0.233E+02   -.546E+00 0.948E+00 -.234E+02   0.843E-01 -.146E+00 0.144E+00   0.113E-03 -.225E-03 0.110E-02
   0.450E+00 0.779E+00 0.233E+02   -.548E+00 -.949E+00 -.234E+02   0.843E-01 0.146E+00 0.144E+00   0.139E-03 0.240E-03 0.110E-02
   0.901E+00 -.341E-03 0.233E+02   -.109E+01 -.106E-02 -.234E+02   0.169E+00 0.266E-18 0.144E+00   0.252E-03 0.150E-04 0.110E-02
   -.451E+00 0.780E+00 0.233E+02   0.546E+00 -.948E+00 -.234E+02   -.843E-01 0.146E+00 0.144E+00   -.113E-03 0.225E-03 0.110E-02
   -.450E+00 -.779E+00 0.233E+02   0.548E+00 0.949E+00 -.234E+02   -.843E-01 -.146E+00 0.144E+00   -.139E-03 -.240E-03 0.110E-02
 -----------------------------------------------------------------------------------------------
   -.432E-02 -.747E-02 0.382E+01   -.496E-13 0.344E-14 -.405E-12   0.805E-15 0.638E-15 -.375E+01   -.339E-06 -.586E-06 -.504E-03
 
 
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      3.44562        -0.000000     -0.000000      0.151025
      0.00000      0.00000      9.85827        -0.000000     -0.000000      0.151025
     -0.00000      3.75653      3.05801        -0.000000     -0.000000     -0.072386
      3.25325      1.87826      9.47066        -0.000000     -0.000000     -0.072386
     -0.00000      3.75653      9.47066        -0.000000     -0.000000     -0.072386
      3.25325      1.87826      3.05801        -0.000000     -0.000000     -0.072386
      0.00000      0.00000      6.16755        -0.000000     -0.000000     -0.147207
      0.00000      0.00000     12.58020        -0.000000     -0.000000     -0.147207
     -0.00000      3.75653      0.23325        -0.000000     -0.000000      0.076472
      3.25325      1.87826      6.64590        -0.000000     -0.000000      0.076472
     -0.00000      3.75653      6.64590        -0.000000     -0.000000      0.076472
      3.25325      1.87826      0.23325        -0.000000     -0.000000      0.076472
      2.03394      0.00000      2.13110        -0.017464     -0.000000      0.037310
     -1.01697      1.76144      2.13110         0.008732     -0.015124      0.037310
      2.23628      3.87335      2.13110         0.008732      0.015124      0.037310
      4.47256      0.00000      8.54375         0.017464     -0.000000      0.037310
     -2.23628      3.87335      8.54375        -0.008732      0.015124      0.037310
      1.01697      1.76144      8.54375        -0.008732     -0.015124      0.037310
      4.20548      0.00000      4.30935        -0.007880     -0.000000     -0.041490
     -2.10274      3.64205      4.30935         0.003940     -0.006825     -0.041490
      1.15051      1.99274      4.30935         0.003940      0.006825     -0.041490
      2.30102      0.00000     10.72200         0.007880     -0.000000     -0.041490
     -1.15051      1.99274     10.72200        -0.003940      0.006825     -0.041490
      2.10274      3.64205     10.72200        -0.003940     -0.006825     -0.041490
      2.16883      0.00000      0.00832         0.025612     -0.000000      0.000183
     -1.08442      1.87826      0.00832        -0.012806      0.022180      0.000183
      2.16883      3.75653      0.00832        -0.012806     -0.022180      0.000183
      4.33767      0.00000      6.42097        -0.025612     -0.000000      0.000183
     -2.16883      3.75653      6.42097         0.012806     -0.022180      0.000183
      1.08442      1.87826      6.42097         0.012806      0.022180      0.000183
 -----------------------------------------------------------------------------------
    total drift:                               -0.004316     -0.007475      0.061055


--------------------------------------------------------------------------------------------------------



  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
  ---------------------------------------------------
  free  energy   TOTEN  =      -219.53574997 eV

  energy  without entropy=     -219.53574997  energy(sigma->0) =     -219.53574997
 


--------------------------------------------------------------------------------------------------------


    POTLOK:  cpu time    0.2701: real time    0.2705


--------------------------------------------------------------------------------------------------------


     LOOP+:  cpu time  502.7103: real time  503.7685
    4ORBIT:  cpu time    0.0000: real time    0.0000
 


 total charge     
 
# of ion       s       p       d       f       tot
--------------------------------------------------
    1        2.021   5.716   0.624   0.249   8.611
    2        2.021   5.716   0.624   0.249   8.611
    3        2.020   5.715   0.617   0.253   8.605
    4        2.020   5.715   0.617   0.253   8.605
    5        2.020   5.715   0.617   0.253   8.605
    6        2.020   5.715   0.617   0.253   8.605
    7        1.243   2.388   0.000   0.000   3.631
    8        1.243   2.388   0.000   0.000   3.631
    9        1.244   2.383   0.000   0.000   3.627
   10        1.244   2.383   0.000   0.000   3.627
   11        1.244   2.383   0.000   0.000   3.627
   12        1.244   2.383   0.000   0.000   3.627
   13        1.570   3.463   0.000   0.000   5.033
   14        1.570   3.457   0.000   0.000   5.028
   15        1.571   3.465   0.000   0.000   5.036
   16        1.570   3.463   0.000   0.000   5.033
   17        1.570   3.457   0.000   0.000   5.028
   18        1.571   3.465   0.000   0.000   5.036
   19        1.569   3.452   0.000   0.000   5.021
   20        1.569   3.459   0.000   0.000   5.028
   21        1.569   3.449   0.000   0.000   5.018
   22        1.569   3.452   0.000   0.000   5.021
   23        1.569   3.459   0.000   0.000   5.028
   24        1.569   3.449   0.000   0.000   5.018
   25        0.964   1.021  10.219   0.000  12.205
   26        0.965   1.021  10.219   0.000  12.204
   27        0.964   1.021  10.219   0.000  12.205
   28        0.964   1.021  10.219   0.000  12.205
   29        0.965   1.021  10.219   0.000  12.204
   30        0.964   1.021  10.219   0.000  12.205
--------------------------------------------------
tot         44.206  96.217  65.032   1.511 206.966

 


 magnetization (x)
 
# of ion       s       p       d       f       tot
--------------------------------------------------
    1        0.000  -0.000   0.000   0.000   0.000
    2        0.000  -0.000   0.000   0.000   0.000
    3        0.000  -0.000   0.000   0.000   0.000
    4        0.000  -0.000   0.000   0.000   0.000
    5        0.000  -0.000   0.000   0.000   0.000
    6        0.000  -0.000   0.000   0.000   0.000
    7        0.000  -0.000   0.000   0.000  -0.000
    8        0.000  -0.000   0.000   0.000  -0.000
    9        0.000  -0.000   0.000   0.000  -0.000
   10        0.000  -0.000   0.000   0.000  -0.000
   11        0.000  -0.000   0.000   0.000  -0.000
   12        0.000  -0.000   0.000   0.000  -0.000
   13       -0.000  -0.000   0.000   0.000  -0.000
   14       -0.000  -0.000   0.000   0.000  -0.000
   15       -0.000  -0.000   0.000   0.000  -0.000
   16       -0.000  -0.000   0.000   0.000  -0.000
   17       -0.000  -0.000   0.000   0.000  -0.000
   18       -0.000  -0.000   0.000   0.000  -0.000
   19       -0.000  -0.000   0.000   0.000  -0.000
   20       -0.000  -0.000   0.000   0.000  -0.000
   21       -0.000  -0.000   0.000   0.000  -0.000
   22       -0.000  -0.000   0.000   0.000  -0.000
   23       -0.000  -0.000   0.000   0.000  -0.000
   24       -0.000  -0.000   0.000   0.000  -0.000
   25        0.000   0.000   0.000   0.000   0.000
   26        0.000   0.000   0.000   0.000   0.000
   27        0.000   0.000   0.000   0.000   0.000
   28        0.000   0.000   0.000   0.000   0.000
   29        0.000   0.000   0.000   0.000   0.000
   30        0.000   0.000   0.000   0.000   0.000
--------------------------------------------------
tot          0.000  -0.000   0.000   0.000  -0.000

 
 BZINTS: Fermi energy:  6.140520;252.000000 electrons
         Band energy:***********;  BLOECHL correction:  0.000000

 total amount of memory used by VASP MPI-rank0    93604. kBytes
=======================================================================

   base      :      30000. kBytes
   nonlr-proj:      27111. kBytes
   fftplans  :       2081. kBytes
   grid      :      21338. kBytes
   one-center:       2949. kBytes
   wavefun   :      10125. kBytes
 
  
  
 General timing and accounting informations for this job:
 ========================================================
  
                  Total CPU time used (sec):      618.098
                            User time (sec):      604.577
                          System time (sec):       13.521
                         Elapsed time (sec):      620.743
  
                   Maximum memory used (kb):      769716.
                   Average memory used (kb):           0.
  
                          Minor page faults:        44184
                          Major page faults:            0
                 Voluntary context switches:         7138
