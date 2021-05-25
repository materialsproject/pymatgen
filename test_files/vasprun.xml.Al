<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <generator>
  <i name="program" type="string">vasp </i>
  <i name="version" type="string">5.4.4.18Apr17-6-g9f103f2a35  </i>
  <i name="subversion" type="string">(build Sep 10 2019 12:14:21) complex            parallel </i>
  <i name="platform" type="string">LinuxIFC </i>
  <i name="date" type="string">2020 10 23 </i>
  <i name="time" type="string">11:27:57 </i>
 </generator>
 <incar>
  <i type="string" name="PREC">accurate</i>
  <i type="string" name="ALGO"> All</i>
  <i type="int" name="ISPIN">     2</i>
  <i type="int" name="ICHARG">     1</i>
  <i type="int" name="NELM">   200</i>
  <i type="int" name="IBRION">     2</i>
  <i name="EDIFF">      0.00001000</i>
  <i name="EDIFFG">     -0.02000000</i>
  <i type="int" name="NSW">    99</i>
  <i type="int" name="ISIF">     3</i>
  <i name="ENCUT">    680.00000000</i>
  <i name="ENAUG">   1360.00000000</i>
  <v name="MAGMOM">      0.60000000</v>
  <i type="string" name="LREAL"> Auto</i>
  <i type="int" name="ISMEAR">     2</i>
  <i name="SIGMA">      0.50000000</i>
  <i name="KSPACING">      0.22000000</i>
  <i type="logical" name="LWAVE"> F  </i>
  <i type="logical" name="LCHARG"> T  </i>
  <i type="logical" name="LVTOT"> T  </i>
  <i type="int" name="LORBIT">    11</i>
  <i type="logical" name="LELF"> T  </i>
  <i type="logical" name="LASPH"> T  </i>
  <i type="logical" name="LAECHG"> T  </i>
  <v type="int" name="KPOINT_BSE">    -1     0     0     0</v>
  <i type="int" name="NELM">   200</i>
  <i type="string" name="METAGGA"> R2scan</i>
  <i type="logical" name="LASPH"> T  </i>
  <i type="int" name="LMIXTAU">    -1</i>
 </incar>
 <structure name="primitive_cell" >
  <crystal>
   <varray name="basis" >
    <v>       0.00000000       2.00554200       2.00554200 </v>
    <v>       2.00554200       0.00000000       2.00554200 </v>
    <v>       2.00554200       2.00554200       0.00000000 </v>
   </varray>
   <i name="volume">     16.13337691 </i>
   <varray name="rec_basis" >
    <v>      -0.24930916       0.24930916       0.24930916 </v>
    <v>       0.24930916      -0.24930916       0.24930916 </v>
    <v>       0.24930916       0.24930916      -0.24930916 </v>
   </varray>
  </crystal>
  <varray name="positions" >
   <v>       0.00000000       0.00000000       0.00000000 </v>
  </varray>
 </structure>
 <varray name="primitive_index" >
  <v type="int" >        1 </v>
 </varray>
 <kpoints>
  <generation param="Gamma">
   <v type="int" name="divisions">      13       13       13 </v>
   <v name="usershift">      0.00000000       0.00000000       0.00000000 </v>
   <v name="genvec1">      0.07692308       0.00000000       0.00000000 </v>
   <v name="genvec2">      0.00000000       0.07692308       0.00000000 </v>
   <v name="genvec3">      0.00000000       0.00000000       0.07692308 </v>
   <v name="shift">      0.00000000       0.00000000       0.00000000 </v>
  </generation>
  <varray name="kpointlist" >
   <v>       0.00000000       0.00000000       0.00000000 </v>
   <v>       0.07692308       0.00000000       0.00000000 </v>
   <v>       0.15384615       0.00000000       0.00000000 </v>
   <v>       0.23076923       0.00000000       0.00000000 </v>
   <v>       0.30769231       0.00000000       0.00000000 </v>
   <v>       0.38461538       0.00000000       0.00000000 </v>
   <v>       0.46153846       0.00000000       0.00000000 </v>
   <v>       0.07692308       0.07692308       0.00000000 </v>
   <v>       0.15384615       0.07692308       0.00000000 </v>
   <v>       0.23076923       0.07692308       0.00000000 </v>
   <v>       0.30769231       0.07692308       0.00000000 </v>
   <v>       0.38461538       0.07692308       0.00000000 </v>
   <v>       0.46153846       0.07692308       0.00000000 </v>
   <v>      -0.46153846       0.07692308       0.00000000 </v>
   <v>      -0.38461538       0.07692308       0.00000000 </v>
   <v>      -0.30769231       0.07692308       0.00000000 </v>
   <v>      -0.23076923       0.07692308       0.00000000 </v>
   <v>      -0.15384615       0.07692308       0.00000000 </v>
   <v>      -0.07692308       0.07692308       0.00000000 </v>
   <v>       0.15384615       0.15384615       0.00000000 </v>
   <v>       0.23076923       0.15384615       0.00000000 </v>
   <v>       0.30769231       0.15384615       0.00000000 </v>
   <v>       0.38461538       0.15384615       0.00000000 </v>
   <v>       0.46153846       0.15384615       0.00000000 </v>
   <v>      -0.46153846       0.15384615       0.00000000 </v>
   <v>      -0.38461538       0.15384615       0.00000000 </v>
   <v>      -0.30769231       0.15384615       0.00000000 </v>
   <v>      -0.23076923       0.15384615       0.00000000 </v>
   <v>      -0.15384615       0.15384615       0.00000000 </v>
   <v>       0.23076923       0.23076923       0.00000000 </v>
   <v>       0.30769231       0.23076923       0.00000000 </v>
   <v>       0.38461538       0.23076923       0.00000000 </v>
   <v>       0.46153846       0.23076923       0.00000000 </v>
   <v>      -0.46153846       0.23076923       0.00000000 </v>
   <v>      -0.38461538       0.23076923       0.00000000 </v>
   <v>      -0.30769231       0.23076923       0.00000000 </v>
   <v>      -0.23076923       0.23076923       0.00000000 </v>
   <v>       0.30769231       0.30769231       0.00000000 </v>
   <v>       0.38461538       0.30769231       0.00000000 </v>
   <v>       0.46153846       0.30769231       0.00000000 </v>
   <v>      -0.46153846       0.30769231       0.00000000 </v>
   <v>      -0.38461538       0.30769231       0.00000000 </v>
   <v>      -0.30769231       0.30769231       0.00000000 </v>
   <v>       0.38461538       0.38461538       0.00000000 </v>
   <v>       0.46153846       0.38461538       0.00000000 </v>
   <v>      -0.46153846       0.38461538       0.00000000 </v>
   <v>      -0.38461538       0.38461538       0.00000000 </v>
   <v>       0.46153846       0.46153846       0.00000000 </v>
   <v>      -0.46153846       0.46153846       0.00000000 </v>
   <v>       0.23076923       0.15384615       0.07692308 </v>
   <v>       0.30769231       0.15384615       0.07692308 </v>
   <v>       0.38461538       0.15384615       0.07692308 </v>
   <v>       0.46153846       0.15384615       0.07692308 </v>
   <v>      -0.46153846       0.15384615       0.07692308 </v>
   <v>       0.30769231       0.23076923       0.07692308 </v>
   <v>       0.38461538       0.23076923       0.07692308 </v>
   <v>       0.46153846       0.23076923       0.07692308 </v>
   <v>      -0.46153846       0.23076923       0.07692308 </v>
   <v>      -0.38461538       0.23076923       0.07692308 </v>
   <v>      -0.30769231       0.23076923       0.07692308 </v>
   <v>      -0.23076923       0.23076923       0.07692308 </v>
   <v>      -0.15384615       0.23076923       0.07692308 </v>
   <v>       0.38461538       0.30769231       0.07692308 </v>
   <v>       0.46153846       0.30769231       0.07692308 </v>
   <v>      -0.46153846       0.30769231       0.07692308 </v>
   <v>      -0.38461538       0.30769231       0.07692308 </v>
   <v>      -0.30769231       0.30769231       0.07692308 </v>
   <v>      -0.23076923       0.30769231       0.07692308 </v>
   <v>       0.46153846       0.38461538       0.07692308 </v>
   <v>      -0.46153846       0.38461538       0.07692308 </v>
   <v>      -0.38461538       0.38461538       0.07692308 </v>
   <v>      -0.30769231       0.38461538       0.07692308 </v>
   <v>      -0.46153846       0.46153846       0.07692308 </v>
   <v>      -0.38461538       0.46153846       0.07692308 </v>
   <v>       0.46153846       0.30769231       0.15384615 </v>
   <v>      -0.46153846       0.30769231       0.15384615 </v>
   <v>      -0.38461538       0.30769231       0.15384615 </v>
   <v>      -0.46153846       0.38461538       0.15384615 </v>
   <v>      -0.38461538       0.38461538       0.15384615 </v>
   <v>      -0.30769231       0.38461538       0.15384615 </v>
   <v>      -0.23076923       0.38461538       0.15384615 </v>
   <v>      -0.38461538       0.46153846       0.15384615 </v>
   <v>      -0.30769231       0.46153846       0.15384615 </v>
   <v>      -0.30769231       0.46153846       0.23076923 </v>
  </varray>
  <varray name="weights" >
   <v>       0.00045517 </v>
   <v>       0.00364133 </v>
   <v>       0.00364133 </v>
   <v>       0.00364133 </v>
   <v>       0.00364133 </v>
   <v>       0.00364133 </v>
   <v>       0.00364133 </v>
   <v>       0.00273100 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.00546199 </v>
   <v>       0.00273100 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.00546199 </v>
   <v>       0.00273100 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.00546199 </v>
   <v>       0.00273100 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.00546199 </v>
   <v>       0.00273100 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.00546199 </v>
   <v>       0.00273100 </v>
   <v>       0.00546199 </v>
   <v>       0.01092399 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.01092399 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.01092399 </v>
   <v>       0.02184797 </v>
   <v>       0.02184797 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
   <v>       0.01092399 </v>
  </varray>
 </kpoints>
 <parameters>
  <separator name="general" >
   <i type="string" name="SYSTEM">unknown system</i>
   <i type="logical" name="LCOMPAT"> F  </i>
  </separator>
  <separator name="electronic" >
   <i type="string" name="PREC">accura</i>
   <i name="ENMAX">    680.00000000</i>
   <i name="ENAUG">   1360.00000000</i>
   <i name="EDIFF">      0.00001000</i>
   <i type="int" name="IALGO">    58</i>
   <i type="int" name="IWAVPR">    11</i>
   <i type="int" name="NBANDS">     5</i>
   <i name="NELECT">      3.00000000</i>
   <i type="int" name="TURBO">     0</i>
   <i type="int" name="IRESTART">     0</i>
   <i type="int" name="NREBOOT">     0</i>
   <i type="int" name="NMIN">     0</i>
   <i name="EREF">      0.00000000</i>
   <separator name="electronic smearing" >
    <i type="int" name="ISMEAR">     2</i>
    <i name="SIGMA">      0.50000000</i>
    <i name="KSPACING">      0.22000000</i>
    <i type="logical" name="KGAMMA"> T  </i>
   </separator>
   <separator name="electronic projectors" >
    <i type="logical" name="LREAL"> T  </i>
    <v name="ROPT">     -0.00025000</v>
    <i type="int" name="LMAXPAW">  -100</i>
    <i type="int" name="LMAXMIX">     2</i>
    <i type="logical" name="NLSPLINE"> F  </i>
   </separator>
   <separator name="electronic startup" >
    <i type="int" name="ISTART">     0</i>
    <i type="int" name="ICHARG">     1</i>
    <i type="int" name="INIWAV">     1</i>
   </separator>
   <separator name="electronic spin" >
    <i type="int" name="ISPIN">     2</i>
    <i type="logical" name="LNONCOLLINEAR"> F  </i>
    <v name="MAGMOM">      0.60000000</v>
    <i name="NUPDOWN">     -1.00000000</i>
    <i type="logical" name="LSORBIT"> F  </i>
    <v name="SAXIS">      0.00000000      0.00000000      1.00000000</v>
    <i type="logical" name="LSPIRAL"> F  </i>
    <v name="QSPIRAL">      0.00000000      0.00000000      0.00000000</v>
    <i type="logical" name="LZEROZ"> F  </i>
   </separator>
   <separator name="electronic exchange-correlation" >
    <i type="logical" name="LASPH"> T  </i>
    <i type="logical" name="LMETAGGA"> F  </i>
   </separator>
   <separator name="electronic convergence" >
    <i type="int" name="NELM">   200</i>
    <i type="int" name="NELMDL">    -5</i>
    <i type="int" name="NELMIN">     2</i>
    <i name="ENINI">    680.00000000</i>
    <separator name="electronic convergence detail" >
     <i type="logical" name="LDIAG"> T  </i>
     <i type="logical" name="LSUBROT"> F  </i>
     <i name="WEIMIN">      0.00100000</i>
     <i name="EBREAK">      0.00000050</i>
     <i name="DEPER">      0.30000000</i>
     <i type="int" name="NRMM">     4</i>
     <i name="TIME">      0.40000000</i>
    </separator>
   </separator>
   <separator name="electronic mixer" >
    <i name="AMIX">      0.40000000</i>
    <i name="BMIX">      1.00000000</i>
    <i name="AMIN">      0.10000000</i>
    <i name="AMIX_MAG">      1.60000000</i>
    <i name="BMIX_MAG">      1.00000000</i>
    <separator name="electronic mixer details" >
     <i type="int" name="IMIX">     4</i>
     <i type="logical" name="MIXFIRST"> F  </i>
     <i type="int" name="MAXMIX">   -45</i>
     <i name="WC">    100.00000000</i>
     <i type="int" name="INIMIX">     1</i>
     <i type="int" name="MIXPRE">     1</i>
     <i type="int" name="MREMOVE">     5</i>
    </separator>
   </separator>
   <separator name="electronic dipolcorrection" >
    <i type="logical" name="LDIPOL"> F  </i>
    <i type="logical" name="LMONO"> F  </i>
    <i type="int" name="IDIPOL">     0</i>
    <i name="EPSILON">      1.00000000</i>
    <v name="DIPOL">   -100.00000000   -100.00000000   -100.00000000</v>
    <i name="EFIELD">      0.00000000</i>
   </separator>
  </separator>
  <separator name="grids" >
   <i type="int" name="NGX">    24</i>
   <i type="int" name="NGY">    24</i>
   <i type="int" name="NGZ">    24</i>
   <i type="int" name="NGXF">    48</i>
   <i type="int" name="NGYF">    48</i>
   <i type="int" name="NGZF">    48</i>
   <i type="logical" name="ADDGRID"> F  </i>
  </separator>
  <separator name="ionic" >
   <i type="int" name="NSW">    99</i>
   <i type="int" name="IBRION">     2</i>
   <i type="int" name="MDALGO">     0</i>
   <i type="int" name="ISIF">     3</i>
   <i name="PSTRESS">      0.00000000</i>
   <i name="EDIFFG">     -0.02000000</i>
   <i type="int" name="NFREE">     1</i>
   <i name="POTIM">      0.50000000</i>
   <i name="SMASS">     -3.00000000</i>
   <i name="SCALEE">      1.00000000</i>
  </separator>
  <separator name="ionic md" >
   <i name="TEBEG">      0.00010000</i>
   <i name="TEEND">      0.00010000</i>
   <i type="int" name="NBLOCK">     1</i>
   <i type="int" name="KBLOCK">    99</i>
   <i type="int" name="NPACO">   256</i>
   <i name="APACO">     16.00000000</i>
  </separator>
  <separator name="symmetry" >
   <i type="int" name="ISYM">     2</i>
   <i name="SYMPREC">      0.00001000</i>
  </separator>
  <separator name="dos" >
   <i type="int" name="LORBIT">    11</i>
   <v name="RWIGS">     -1.00000000</v>
   <i type="int" name="NEDOS">   301</i>
   <i name="EMIN">     10.00000000</i>
   <i name="EMAX">    -10.00000000</i>
   <i name="EFERMI">      0.00000000</i>
  </separator>
  <separator name="writing" >
   <i type="int" name="NWRITE">     2</i>
   <i type="logical" name="LWAVE"> F  </i>
   <i type="logical" name="LDOWNSAMPLE"> F  </i>
   <i type="logical" name="LCHARG"> T  </i>
   <i type="logical" name="LPARD"> F  </i>
   <i type="logical" name="LVTOT"> T  </i>
   <i type="logical" name="LVHAR"> F  </i>
   <i type="logical" name="LELF"> T  </i>
   <i type="logical" name="LOPTICS"> F  </i>
   <v name="STM">      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000</v>
  </separator>
  <separator name="performance" >
   <i type="int" name="NPAR">     1</i>
   <i type="int" name="NSIM">     4</i>
   <i type="int" name="NBLK">    -1</i>
   <i type="logical" name="LPLANE"> T  </i>
   <i type="logical" name="LSCALAPACK"> F  </i>
   <i type="logical" name="LSCAAWARE"> F  </i>
   <i type="logical" name="LSCALU"> F  </i>
   <i type="logical" name="LASYNC"> F  </i>
   <i type="logical" name="LORBITALREAL"> F  </i>
  </separator>
  <separator name="miscellaneous" >
   <i type="int" name="IDIOT">     3</i>
   <i type="int" name="ISPECIAL">     0</i>
   <i type="logical" name="LMUSIC"> F  </i>
   <v name="POMASS">     26.98100000</v>
   <v name="DARWINR">      0.00000000</v>
   <v name="DARWINV">      1.00000000</v>
   <i type="logical" name="LCORR"> T  </i>
  </separator>
  <i type="logical" name="GGA_COMPAT"> T  </i>
  <i type="logical" name="LBERRY"> F  </i>
  <i type="int" name="ICORELEVEL">     0</i>
  <i type="logical" name="LDAU"> F  </i>
  <i type="int" name="I_CONSTRAINED_M">     0</i>
  <separator name="electronic exchange-correlation" >
   <i type="string" name="GGA">--</i>
   <i type="int" name="VOSKOWN">     0</i>
   <i type="logical" name="LHFCALC"> F  </i>
   <i type="string" name="PRECFOCK"></i>
   <i type="logical" name="LSYMGRAD"> F  </i>
   <i type="logical" name="LHFONE"> F  </i>
   <i type="logical" name="LRHFCALC"> F  </i>
   <i type="logical" name="LTHOMAS"> F  </i>
   <i type="logical" name="LMODELHF"> F  </i>
   <i name="ENCUT4O">     -1.00000000</i>
   <i type="int" name="EXXOEP">     0</i>
   <i type="int" name="FOURORBIT">     0</i>
   <i name="AEXX">      0.00000000</i>
   <i name="HFALPHA">      0.00000000</i>
   <i name="MCALPHA">      0.00000000</i>
   <i name="ALDAX">      1.00000000</i>
   <i name="AGGAX">      1.00000000</i>
   <i name="ALDAC">      1.00000000</i>
   <i name="AGGAC">      1.00000000</i>
   <i type="int" name="NKREDX">     1</i>
   <i type="int" name="NKREDY">     1</i>
   <i type="int" name="NKREDZ">     1</i>
   <i type="logical" name="SHIFTRED"> F  </i>
   <i type="logical" name="ODDONLY"> F  </i>
   <i type="logical" name="EVENONLY"> F  </i>
   <i type="int" name="LMAXFOCK">     0</i>
   <i type="int" name="NMAXFOCKAE">     0</i>
   <i type="logical" name="LFOCKAEDFT"> F  </i>
   <i name="HFSCREEN">      0.00000000</i>
   <i name="HFSCREENC">      0.00000000</i>
   <i type="int" name="NBANDSGWLOW">     0</i>
  </separator>
  <separator name="vdW DFT" >
   <i type="logical" name="LUSE_VDW"> F  </i>
   <i name="Zab_VDW">     -0.84910000</i>
   <i name="PARAM1">      0.12340000</i>
   <i name="PARAM2">      1.00000000</i>
   <i name="PARAM3">      0.00000000</i>
  </separator>
  <separator name="model GW" >
   <i type="int" name="MODEL_GW">     0</i>
   <i name="MODEL_EPS0">      9.68002614</i>
   <i name="MODEL_ALPHA">      1.00000000</i>
  </separator>
  <separator name="linear response parameters" >
   <i type="logical" name="LEPSILON"> F  </i>
   <i type="logical" name="LRPA"> F  </i>
   <i type="logical" name="LNABLA"> F  </i>
   <i type="logical" name="LVEL"> F  </i>
   <i type="int" name="KINTER">     0</i>
   <i name="CSHIFT">      0.10000000</i>
   <i name="OMEGAMAX">     -1.00000000</i>
   <i name="DEG_THRESHOLD">      0.00200000</i>
   <i name="RTIME">     -0.10000000</i>
   <i name="WPLASMAI">      0.00000000</i>
   <v name="DFIELD">      0.00000000      0.00000000      0.00000000</v>
   <v name="WPLASMA">      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000      0.00000000</v>
  </separator>
  <separator name="orbital magnetization" >
   <i type="logical" name="NUCIND"> F  </i>
   <v name="MAGPOS">      0.00000000      0.00000000      0.00000000</v>
   <i type="logical" name="LNICSALL"> T  </i>
   <i type="logical" name="ORBITALMAG"> F  </i>
   <i type="logical" name="LMAGBLOCH"> F  </i>
   <i type="logical" name="LCHIMAG"> F  </i>
   <i type="logical" name="LGAUGE"> T  </i>
   <i type="int" name="MAGATOM">     0</i>
   <v name="MAGDIPOL">      0.00000000      0.00000000      0.00000000</v>
   <v name="AVECCONST">      0.00000000      0.00000000      0.00000000</v>
  </separator>
  <separator name="response functions" >
   <i type="logical" name="LADDER"> F  </i>
   <i type="logical" name="LRPAFORCE"> F  </i>
   <i type="logical" name="LFXC"> F  </i>
   <i type="logical" name="LHARTREE"> T  </i>
   <i type="int" name="IBSE">     0</i>
   <v type="int" name="KPOINT">    -1     0     0     0</v>
   <i type="logical" name="LTCTE"> F  </i>
   <i type="logical" name="LTETE"> F  </i>
   <i type="logical" name="LTRIPLET"> F  </i>
   <i type="logical" name="LFXCEPS"> F  </i>
   <i type="logical" name="LFXHEG"> F  </i>
   <i type="int" name="NATURALO">     2</i>
   <i type="logical" name="L2ORDER"> F  </i>
   <i type="logical" name="LMP2LT"> F  </i>
   <i type="logical" name="LUSEW"> F  </i>
   <i name="ENCUTGW">     -2.00000000</i>
   <i name="ENCUTGWSOFT">     -2.00000000</i>
   <i name="ENCUTLF">     -1.00000000</i>
   <i type="int" name="LMAXMP2">    -1</i>
   <i name="SCISSOR">      0.00000000</i>
   <i type="int" name="NOMEGA">     0</i>
   <i type="int" name="NOMEGAR">     0</i>
   <i type="int" name="NBANDSGW">    -1</i>
   <i type="int" name="NBANDSO">    -1</i>
   <i type="int" name="NBANDSV">    -1</i>
   <i type="int" name="NELM">   200</i>
   <i type="int" name="NELMHF">     1</i>
   <i type="int" name="DIM">     3</i>
   <i type="int" name="ANTIRES">     0</i>
   <i name="OMEGAMAX">    -30.00000000</i>
   <i name="OMEGAMIN">    -30.00000000</i>
   <i name="OMEGATL">   -200.00000000</i>
   <i type="int" name="OMEGAGRID">   140</i>
   <i name="CSHIFT">     -0.10000000</i>
   <i type="logical" name="SELFENERGY"> F  </i>
   <i type="logical" name="LSPECTRAL"> F  </i>
   <i type="logical" name="LSPECTRALGW"> F  </i>
   <i type="logical" name="LSINGLES"> F  </i>
   <i type="logical" name="LFERMIGW"> F  </i>
   <i type="logical" name="ODDONLYGW"> F  </i>
   <i type="logical" name="EVENONLYGW"> F  </i>
   <i type="int" name="NKREDLFX">     1</i>
   <i type="int" name="NKREDLFY">     1</i>
   <i type="int" name="NKREDLFZ">     1</i>
   <i type="int" name="MAXMEM">  2800</i>
   <i type="int" name="TELESCOPE">     0</i>
   <i type="int" name="TAUPAR">     1</i>
   <i type="int" name="OMEGAPAR">    -1</i>
   <i name="LAMBDA">      1.00000000</i>
  </separator>
  <separator name="External order field" >
   <i name="OFIELD_KAPPA">      0.00000000</i>
   <v name="OFIELD_K">      0.00000000      0.00000000      0.00000000</v>
   <i name="OFIELD_Q6_NEAR">      0.00000000</i>
   <i name="OFIELD_Q6_FAR">      0.00000000</i>
   <i name="OFIELD_A">      0.00000000</i>
  </separator>
 </parameters>
 <atominfo>
  <atoms>       1 </atoms>
  <types>       1 </types>
  <array name="atoms" >
   <dimension dim="1">ion</dimension>
   <field type="string">element</field>
   <field type="int">atomtype</field>
   <set>
    <rc><c>Al</c><c>   1</c></rc>
   </set>
  </array>
  <array name="atomtypes" >
   <dimension dim="1">type</dimension>
   <field type="int">atomspertype</field>
   <field type="string">element</field>
   <field>mass</field>
   <field>valence</field>
   <field type="string">pseudopotential</field>
   <set>
    <rc><c>   1</c><c>Al</c><c>     26.98100000</c><c>      3.00000000</c><c>  PAW_PBE Al 04Jan2001                  </c></rc>
   </set>
  </array>
 </atominfo>
 <structure name="initialpos" >
  <crystal>
   <varray name="basis" >
    <v>       0.00000000       2.00554200       2.00554200 </v>
    <v>       2.00554200       0.00000000       2.00554200 </v>
    <v>       2.00554200       2.00554200       0.00000000 </v>
   </varray>
   <i name="volume">     16.13337691 </i>
   <varray name="rec_basis" >
    <v>      -0.24930916       0.24930916       0.24930916 </v>
    <v>       0.24930916      -0.24930916       0.24930916 </v>
    <v>       0.24930916       0.24930916      -0.24930916 </v>
   </varray>
  </crystal>
  <varray name="positions" >
   <v>       0.00000000       0.00000000       0.00000000 </v>
  </varray>
 </structure>
 <calculation>
  <scstep>
   <time name="total">    1.02    1.02</time>
   <energy>
    <i name="alphaZ">     -0.20819013 </i>
    <i name="ewald">    -74.06842420 </i>
    <i name="hartreedc">     -0.10831478 </i>
    <i name="XCdc">    -26.08230095 </i>
    <i name="pawpsdc">    152.40535911 </i>
    <i name="pawaedc">   -120.46669884 </i>
    <i name="eentropy">      0.00117381 </i>
    <i name="bandstr">     13.88318686 </i>
    <i name="atom">     53.53856935 </i>
    <i name="e_fr_energy">     -1.10563977 </i>
    <i name="e_wo_entrp">     -1.10681358 </i>
    <i name="e_0_energy">     -1.10593322 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    0.98    0.98</time>
   <energy>
    <i name="e_fr_energy">     -3.71296439 </i>
    <i name="e_wo_entrp">     -3.71277505 </i>
    <i name="e_0_energy">     -3.71291705 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.05    1.05</time>
   <energy>
    <i name="e_fr_energy">     -3.74033193 </i>
    <i name="e_wo_entrp">     -3.74035137 </i>
    <i name="e_0_energy">     -3.74033679 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.02    1.03</time>
   <energy>
    <i name="e_fr_energy">     -3.74043421 </i>
    <i name="e_wo_entrp">     -3.74045489 </i>
    <i name="e_0_energy">     -3.74043938 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.16    1.16</time>
   <energy>
    <i name="e_fr_energy">     -3.74043475 </i>
    <i name="e_wo_entrp">     -3.74045544 </i>
    <i name="e_0_energy">     -3.74043992 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.44    1.45</time>
   <energy>
    <i name="alphaZ">     -0.20819013 </i>
    <i name="ewald">    -74.06842420 </i>
    <i name="hartreedc">     -0.10932828 </i>
    <i name="XCdc">    -26.08362149 </i>
    <i name="pawpsdc">    146.46739395 </i>
    <i name="pawaedc">   -114.52883190 </i>
    <i name="eentropy">      0.00002059 </i>
    <i name="bandstr">     11.25197125 </i>
    <i name="atom">     53.53856935 </i>
    <i name="e_fr_energy">     -3.74044075 </i>
    <i name="e_wo_entrp">     -3.74046135 </i>
    <i name="e_0_energy">     -3.74044590 </i>
   </energy>
  </scstep>
  <structure>
   <crystal>
    <varray name="basis" >
     <v>       0.00000000       2.00554200       2.00554200 </v>
     <v>       2.00554200       0.00000000       2.00554200 </v>
     <v>       2.00554200       2.00554200       0.00000000 </v>
    </varray>
    <i name="volume">     16.13337691 </i>
    <varray name="rec_basis" >
     <v>      -0.24930916       0.24930916       0.24930916 </v>
     <v>       0.24930916      -0.24930916       0.24930916 </v>
     <v>       0.24930916       0.24930916      -0.24930916 </v>
    </varray>
   </crystal>
   <varray name="positions" >
    <v>       0.00000000       0.00000000       0.00000000 </v>
   </varray>
  </structure>
  <varray name="forces" >
   <v>       0.00000000       0.00000000       0.00000000 </v>
  </varray>
  <varray name="stress" >
   <v>      15.40962538       0.00000000      -0.00000000 </v>
   <v>       0.00000000      15.40962538       0.00000000 </v>
   <v>       0.00000000       0.00000000      15.40962538 </v>
  </varray>
  <energy>
   <i name="e_fr_energy">     -3.74044075 </i>
   <i name="e_wo_entrp">     -3.74044590 </i>
   <i name="e_0_energy">      0.00002059 </i>
  </energy>
  <time name="totalsc">    8.90    8.98</time>
 </calculation>
 <calculation>
  <scstep>
   <time name="total">    1.22    1.23</time>
   <energy>
    <i name="alphaZ">     -0.20358386 </i>
    <i name="ewald">    -73.51808285 </i>
    <i name="hartreedc">     -0.10775027 </i>
    <i name="XCdc">    -26.12203371 </i>
    <i name="pawpsdc">    142.93126255 </i>
    <i name="pawaedc">   -110.98245168 </i>
    <i name="eentropy">     -0.00006690 </i>
    <i name="bandstr">     10.72444116 </i>
    <i name="atom">     53.53856935 </i>
    <i name="e_fr_energy">     -3.73960873 </i>
    <i name="e_wo_entrp">     -3.73954183 </i>
    <i name="e_0_energy">     -3.73959200 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.24    1.25</time>
   <energy>
    <i name="e_fr_energy">     -3.74168038 </i>
    <i name="e_wo_entrp">     -3.74148322 </i>
    <i name="e_0_energy">     -3.74163109 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.17    1.18</time>
   <energy>
    <i name="e_fr_energy">     -3.74198837 </i>
    <i name="e_wo_entrp">     -3.74174439 </i>
    <i name="e_0_energy">     -3.74192737 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.18    1.19</time>
   <energy>
    <i name="e_fr_energy">     -3.74204064 </i>
    <i name="e_wo_entrp">     -3.74179360 </i>
    <i name="e_0_energy">     -3.74197888 </i>
   </energy>
  </scstep>
  <scstep>
   <time name="total">    1.36    1.36</time>
   <energy>
    <i name="alphaZ">     -0.20358386 </i>
    <i name="ewald">    -73.51808285 </i>
    <i name="hartreedc">     -0.10566225 </i>
    <i name="XCdc">    -26.12907663 </i>
    <i name="pawpsdc">    140.51259603 </i>
    <i name="pawaedc">   -108.57025386 </i>
    <i name="eentropy">     -0.00024694 </i>
    <i name="bandstr">     10.73369816 </i>
    <i name="atom">     53.53856935 </i>
    <i name="e_fr_energy">     -3.74204295 </i>
    <i name="e_wo_entrp">     -3.74179601 </i>
    <i name="e_0_energy">     -3.74198121 </i>
   </energy>
  </scstep>
  <structure>
   <crystal>
    <varray name="basis" >
     <v>       0.00000000       2.02055508       2.02055508 </v>
     <v>       2.02055508       0.00000000       2.02055508 </v>
     <v>       2.02055508       2.02055508      -0.00000000 </v>
    </varray>
    <i name="volume">     16.49840937 </i>
    <varray name="rec_basis" >
     <v>      -0.24745675       0.24745675       0.24745675 </v>
     <v>       0.24745675      -0.24745675       0.24745675 </v>
     <v>       0.24745675       0.24745675      -0.24745675 </v>
    </varray>
   </crystal>
   <varray name="positions" >
    <v>       0.00000000       0.00000000       0.00000000 </v>
   </varray>
  </structure>
  <varray name="forces" >
   <v>      -0.00000000       0.00000000       0.00000000 </v>
  </varray>
  <varray name="stress" >
   <v>      -0.83964176       0.00000000      -0.00000000 </v>
   <v>      -0.00000000      -0.83964176       0.00000000 </v>
   <v>       0.00000000       0.00000000      -0.83964176 </v>
  </varray>
  <energy>
   <i name="e_fr_energy">     -3.74204295 </i>
   <i name="e_wo_entrp">     -3.74198121 </i>
   <i name="e_0_energy">     -0.00024694 </i>
  </energy>
  <time name="totalsc">    7.85    7.88</time>
  <eigenvalues>
   <array>
    <dimension dim="1">band</dimension>
    <dimension dim="2">kpoint</dimension>
    <dimension dim="3">spin</dimension>
    <field>eigene</field>
    <field>occ</field>
    <set>
     <set comment="spin 1">
      <set comment="kpoint 1">
       <r>   -3.1341    1.0000 </r>
       <r>   20.6208    0.0000 </r>
       <r>   20.6208    0.0000 </r>
       <r>   20.6208    0.0000 </r>
       <r>   22.1131    0.0000 </r>
      </set>
      <set comment="kpoint 2">
       <r>   -2.9732    1.0000 </r>
       <r>   18.2625    0.0000 </r>
       <r>   20.1120    0.0000 </r>
       <r>   20.1120    0.0000 </r>
       <r>   22.7761    0.0000 </r>
      </set>
      <set comment="kpoint 3">
       <r>   -2.4918    1.0000 </r>
       <r>   15.1545    0.0000 </r>
       <r>   19.4447    0.0000 </r>
       <r>   19.4447    0.0000 </r>
       <r>   22.6316    0.0000 </r>
      </set>
      <set comment="kpoint 4">
       <r>   -1.6930    1.0000 </r>
       <r>   12.1586    0.0000 </r>
       <r>   18.9854    0.0000 </r>
       <r>   18.9854    0.0000 </r>
       <r>   22.3991    0.0000 </r>
      </set>
      <set comment="kpoint 5">
       <r>   -0.5822    1.0000 </r>
       <r>    9.3722    0.0005 </r>
       <r>   18.7803    0.0000 </r>
       <r>   18.7805    0.0000 </r>
       <r>   21.7676    0.0000 </r>
      </set>
      <set comment="kpoint 6">
       <r>    0.8335    1.0000 </r>
       <r>    6.8308    0.9959 </r>
       <r>   18.8406    0.0000 </r>
       <r>   18.8412    0.0000 </r>
       <r>   20.6402    0.0000 </r>
      </set>
      <set comment="kpoint 7">
       <r>    2.5434    1.0000 </r>
       <r>    4.5561    1.0000 </r>
       <r>   19.1696    0.0000 </r>
       <r>   19.1756    0.0000 </r>
       <r>   19.7719    0.0000 </r>
      </set>
      <set comment="kpoint 8">
       <r>   -2.9196    1.0000 </r>
       <r>   19.0625    0.0000 </r>
       <r>   19.0625    0.0000 </r>
       <r>   19.3530    0.0000 </r>
       <r>   21.6182    0.0000 </r>
      </set>
      <set comment="kpoint 9">
       <r>   -2.5452    1.0000 </r>
       <r>   16.2579    0.0000 </r>
       <r>   18.2305    0.0000 </r>
       <r>   19.4493    0.0000 </r>
       <r>   21.9393    0.0000 </r>
      </set>
      <set comment="kpoint 10">
       <r>   -1.8524    1.0000 </r>
       <r>   13.2453    0.0000 </r>
       <r>   17.6567    0.0000 </r>
       <r>   19.1235    0.0000 </r>
       <r>   21.3429    0.0000 </r>
      </set>
      <set comment="kpoint 11">
       <r>   -0.8459    1.0000 </r>
       <r>   10.3980    0.0000 </r>
       <r>   17.3504    0.0000 </r>
       <r>   18.8807    0.0000 </r>
       <r>   19.8568    0.0000 </r>
      </set>
      <set comment="kpoint 12">
       <r>    0.4677    1.0000 </r>
       <r>    7.7792    0.8110 </r>
       <r>   17.3148    0.0000 </r>
       <r>   18.5785    0.0000 </r>
       <r>   18.9007    0.0000 </r>
      </set>
      <set comment="kpoint 13">
       <r>    2.0796    1.0000 </r>
       <r>    5.4172    1.0000 </r>
       <r>   17.5514    0.0000 </r>
       <r>   17.6063    0.0000 </r>
       <r>   19.1519    0.0000 </r>
      </set>
      <set comment="kpoint 14">
       <r>    3.3035    1.0000 </r>
       <r>    4.0052    1.0000 </r>
       <r>   16.8983    0.0000 </r>
       <r>   18.0612    0.0000 </r>
       <r>   19.6932    0.0000 </r>
      </set>
      <set comment="kpoint 15">
       <r>    1.5100    1.0000 </r>
       <r>    6.1739    1.0000 </r>
       <r>   16.4703    0.0000 </r>
       <r>   18.8428    0.0000 </r>
       <r>   20.5110    0.0000 </r>
      </set>
      <set comment="kpoint 16">
       <r>   -0.0038    1.0000 </r>
       <r>    8.6224   -0.0218 </r>
       <r>   16.3249    0.0000 </r>
       <r>   19.8771    0.0000 </r>
       <r>   20.5167    0.0000 </r>
      </set>
      <set comment="kpoint 17">
       <r>   -1.2159    1.0000 </r>
       <r>   11.3163    0.0000 </r>
       <r>   16.4641    0.0000 </r>
       <r>   20.5551    0.0000 </r>
       <r>   21.2450    0.0000 </r>
      </set>
      <set comment="kpoint 18">
       <r>   -2.1185    1.0000 </r>
       <r>   14.2062    0.0000 </r>
       <r>   16.8998    0.0000 </r>
       <r>   20.8873    0.0000 </r>
       <r>   22.6024    0.0000 </r>
      </set>
      <set comment="kpoint 19">
       <r>   -2.7056    1.0000 </r>
       <r>   16.9661    0.0000 </r>
       <r>   17.8588    0.0000 </r>
       <r>   21.0378    0.0000 </r>
       <r>   22.3633    0.0000 </r>
      </set>
      <set comment="kpoint 20">
       <r>   -2.2783    1.0000 </r>
       <r>   17.2814    0.0000 </r>
       <r>   17.2814    0.0000 </r>
       <r>   17.4619    0.0000 </r>
       <r>   20.5141    0.0000 </r>
      </set>
      <set comment="kpoint 21">
       <r>   -1.6930    1.0000 </r>
       <r>   14.5332    0.0000 </r>
       <r>   16.6041    0.0000 </r>
       <r>   17.8754    0.0000 </r>
       <r>   19.2376    0.0000 </r>
      </set>
      <set comment="kpoint 22">
       <r>   -0.7932    1.0000 </r>
       <r>   11.6813    0.0000 </r>
       <r>   16.2000    0.0000 </r>
       <r>   17.5061    0.0000 </r>
       <r>   17.8750    0.0000 </r>
      </set>
      <set comment="kpoint 23">
       <r>    0.4150    1.0000 </r>
       <r>    9.0083    0.0046 </r>
       <r>   16.0694    0.0000 </r>
       <r>   16.1799    0.0000 </r>
       <r>   17.7810    0.0000 </r>
      </set>
      <set comment="kpoint 24">
       <r>    1.9234    1.0000 </r>
       <r>    6.5725    0.9991 </r>
       <r>   15.0818    0.0000 </r>
       <r>   16.2127    0.0000 </r>
       <r>   17.9735    0.0000 </r>
      </set>
      <set comment="kpoint 25">
       <r>    3.7068    1.0000 </r>
       <r>    4.4152    1.0000 </r>
       <r>   14.2658    0.0000 </r>
       <r>   16.6303    0.0000 </r>
       <r>   18.4467    0.0000 </r>
      </set>
      <set comment="kpoint 26">
       <r>    2.4919    1.0000 </r>
       <r>    5.8125    1.0000 </r>
       <r>   13.7414    0.0000 </r>
       <r>   17.3224    0.0000 </r>
       <r>   19.2038    0.0000 </r>
      </set>
      <set comment="kpoint 27">
       <r>    0.8854    1.0000 </r>
       <r>    8.1494    0.1093 </r>
       <r>   13.5200    0.0000 </r>
       <r>   18.2870    0.0000 </r>
       <r>   20.2350    0.0000 </r>
      </set>
      <set comment="kpoint 28">
       <r>   -0.4243    1.0000 </r>
       <r>   10.6937    0.0000 </r>
       <r>   13.6386    0.0000 </r>
       <r>   19.5157    0.0000 </r>
       <r>   21.5178    0.0000 </r>
      </set>
      <set comment="kpoint 29">
       <r>   -1.4278    1.0000 </r>
       <r>   13.0656    0.0000 </r>
       <r>   14.4433    0.0000 </r>
       <r>   20.9757    0.0000 </r>
       <r>   22.2850    0.0000 </r>
      </set>
      <set comment="kpoint 30">
       <r>   -1.2160    1.0000 </r>
       <r>   15.8284    0.0000 </r>
       <r>   15.8284    0.0000 </r>
       <r>   15.8293    0.0000 </r>
       <r>   17.2063    0.0000 </r>
      </set>
      <set comment="kpoint 31">
       <r>   -0.4246    1.0000 </r>
       <r>   13.1660    0.0000 </r>
       <r>   15.3287    0.0000 </r>
       <r>   15.4253    0.0000 </r>
       <r>   16.5538    0.0000 </r>
      </set>
      <set comment="kpoint 32">
       <r>    0.6755    1.0000 </r>
       <r>   10.4935    0.0000 </r>
       <r>   13.9343    0.0000 </r>
       <r>   15.1040    0.0000 </r>
       <r>   16.6700    0.0000 </r>
      </set>
      <set comment="kpoint 33">
       <r>    2.0762    1.0000 </r>
       <r>    8.0069    0.3565 </r>
       <r>   12.7214    0.0000 </r>
       <r>   15.1539    0.0000 </r>
       <r>   16.8304    0.0000 </r>
      </set>
      <set comment="kpoint 34">
       <r>    3.7624    1.0000 </r>
       <r>    5.7667    1.0000 </r>
       <r>   11.8000    0.0000 </r>
       <r>   15.4786    0.0000 </r>
       <r>   17.2330    0.0000 </r>
      </set>
      <set comment="kpoint 35">
       <r>    3.7703    1.0000 </r>
       <r>    5.7430    1.0000 </r>
       <r>   11.1866    0.0000 </r>
       <r>   16.0784    0.0000 </r>
       <r>   17.9154    0.0000 </r>
      </set>
      <set comment="kpoint 36">
       <r>    2.0791    1.0000 </r>
       <r>    7.9166    0.5452 </r>
       <r>   10.9288    0.0000 </r>
       <r>   16.9535    0.0000 </r>
       <r>   18.8771    0.0000 </r>
      </set>
      <set comment="kpoint 37">
       <r>    0.6766    1.0000 </r>
       <r>    9.9061    0.0000 </r>
       <r>   11.4181    0.0000 </r>
       <r>   18.1019    0.0000 </r>
       <r>   20.0968    0.0000 </r>
      </set>
      <set comment="kpoint 38">
       <r>    0.2572    1.0000 </r>
       <r>   13.5342    0.0000 </r>
       <r>   14.5772    0.0000 </r>
       <r>   14.7348    0.0000 </r>
       <r>   14.7348    0.0000 </r>
      </set>
      <set comment="kpoint 39">
       <r>    1.2467    1.0000 </r>
       <r>   11.8647    0.0000 </r>
       <r>   12.2444    0.0000 </r>
       <r>   14.4167    0.0000 </r>
       <r>   15.5526    0.0000 </r>
      </set>
      <set comment="kpoint 40">
       <r>    2.5349    1.0000 </r>
       <r>    9.6698    0.0000 </r>
       <r>   10.6384    0.0000 </r>
       <r>   14.3736    0.0000 </r>
       <r>   15.8675    0.0000 </r>
      </set>
      <set comment="kpoint 41">
       <r>    4.1038    1.0000 </r>
       <r>    7.3956    1.0471 </r>
       <r>    9.6030    0.0001 </r>
       <r>   14.6052    0.0000 </r>
       <r>   16.2338    0.0000 </r>
      </set>
      <set comment="kpoint 42">
       <r>    5.3203    1.0000 </r>
       <r>    5.9384    1.0000 </r>
       <r>    8.9330    0.0049 </r>
       <r>   15.1119    0.0000 </r>
       <r>   16.8439    0.0000 </r>
      </set>
      <set comment="kpoint 43">
       <r>    3.5681    1.0000 </r>
       <r>    7.5369    1.0415 </r>
       <r>    9.0263    0.0044 </r>
       <r>   15.8946    0.0000 </r>
       <r>   17.7290    0.0000 </r>
      </set>
      <set comment="kpoint 44">
       <r>    2.1225    1.0000 </r>
       <r>   10.1899    0.0000 </r>
       <r>   13.7349    0.0000 </r>
       <r>   14.0057    0.0000 </r>
       <r>   14.0057    0.0000 </r>
      </set>
      <set comment="kpoint 45">
       <r>    3.2895    1.0000 </r>
       <r>    8.7585   -0.0025 </r>
       <r>   11.6066    0.0000 </r>
       <r>   13.8699    0.0000 </r>
       <r>   14.9490    0.0000 </r>
      </set>
      <set comment="kpoint 46">
       <r>    4.7029    1.0000 </r>
       <r>    7.6607    0.9671 </r>
       <r>    9.3425    0.0007 </r>
       <r>   14.0086    0.0000 </r>
       <r>   15.4595    0.0000 </r>
      </set>
      <set comment="kpoint 47">
       <r>    5.9615    1.0000 </r>
       <r>    7.0099    0.9961 </r>
       <r>    7.5610    1.0335 </r>
       <r>   14.4224    0.0000 </r>
       <r>   16.0315    0.0000 </r>
      </set>
      <set comment="kpoint 48">
       <r>    4.3007    1.0000 </r>
       <r>    7.2695    1.0247 </r>
       <r>   13.3121    0.0000 </r>
       <r>   13.6415    0.0000 </r>
       <r>   13.6415    0.0000 </r>
      </set>
      <set comment="kpoint 49">
       <r>    5.1754    1.0000 </r>
       <r>    6.5036    0.9995 </r>
       <r>   11.4071    0.0000 </r>
       <r>   13.6877    0.0000 </r>
       <r>   14.7469    0.0000 </r>
      </set>
      <set comment="kpoint 50">
       <r>   -2.0652    1.0000 </r>
       <r>   15.0315    0.0000 </r>
       <r>   16.1512    0.0000 </r>
       <r>   20.0199    0.0000 </r>
       <r>   21.1683    0.0000 </r>
      </set>
      <set comment="kpoint 51">
       <r>   -1.2689    1.0000 </r>
       <r>   12.4574    0.0000 </r>
       <r>   15.2797    0.0000 </r>
       <r>   19.5419    0.0000 </r>
       <r>   20.1691    0.0000 </r>
      </set>
      <set comment="kpoint 52">
       <r>   -0.1614    1.0000 </r>
       <r>    9.7382    0.0000 </r>
       <r>   15.0008    0.0000 </r>
       <r>   18.0988    0.0000 </r>
       <r>   20.0874    0.0000 </r>
      </set>
      <set comment="kpoint 53">
       <r>    1.2498    1.0000 </r>
       <r>    7.2214    1.0161 </r>
       <r>   15.0354    0.0000 </r>
       <r>   16.9219    0.0000 </r>
       <r>   20.1896    0.0000 </r>
      </set>
      <set comment="kpoint 54">
       <r>    2.9531    1.0000 </r>
       <r>    4.9602    1.0000 </r>
       <r>   15.3444    0.0000 </r>
       <r>   16.0412    0.0000 </r>
       <r>   20.5287    0.0000 </r>
      </set>
      <set comment="kpoint 55">
       <r>   -1.0045    1.0000 </r>
       <r>   13.4633    0.0000 </r>
       <r>   14.6801    0.0000 </r>
       <r>   17.4148    0.0000 </r>
       <r>   18.5677    0.0000 </r>
      </set>
      <set comment="kpoint 56">
       <r>   -0.0042    1.0000 </r>
       <r>   11.0792    0.0000 </r>
       <r>   13.9496    0.0000 </r>
       <r>   15.8519    0.0000 </r>
       <r>   18.9514    0.0000 </r>
      </set>
      <set comment="kpoint 57">
       <r>    1.3008    1.0000 </r>
       <r>    8.5425   -0.0369 </r>
       <r>   13.8129    0.0000 </r>
       <r>   14.5909    0.0000 </r>
       <r>   19.0846    0.0000 </r>
      </set>
      <set comment="kpoint 58">
       <r>    2.9004    1.0000 </r>
       <r>    6.2143    1.0000 </r>
       <r>   13.4410    0.0000 </r>
       <r>   14.1743    0.0000 </r>
       <r>   19.3910    0.0000 </r>
      </set>
      <set comment="kpoint 59">
       <r>    4.1120    1.0000 </r>
       <r>    4.8161    1.0000 </r>
       <r>   12.7617    0.0000 </r>
       <r>   14.6560    0.0000 </r>
       <r>   19.9480    0.0000 </r>
      </set>
      <set comment="kpoint 60">
       <r>    2.3367    1.0000 </r>
       <r>    6.9472    0.9951 </r>
       <r>   12.3544    0.0000 </r>
       <r>   15.4615    0.0000 </r>
       <r>   20.7617    0.0000 </r>
      </set>
      <set comment="kpoint 61">
       <r>    0.8329    1.0000 </r>
       <r>    9.3076    0.0009 </r>
       <r>   12.2910    0.0000 </r>
       <r>   16.5547    0.0000 </r>
       <r>   21.7689    0.0000 </r>
      </set>
      <set comment="kpoint 62">
       <r>   -0.3717    1.0000 </r>
       <r>   11.4872    0.0000 </r>
       <r>   12.9443    0.0000 </r>
       <r>   17.9221    0.0000 </r>
       <r>   21.6206    0.0000 </r>
      </set>
      <set comment="kpoint 63">
       <r>    0.4665    1.0000 </r>
       <r>   12.2802    0.0000 </r>
       <r>   13.4000    0.0000 </r>
       <r>   13.9191    0.0000 </r>
       <r>   17.4206    0.0000 </r>
      </set>
      <set comment="kpoint 64">
       <r>    1.6620    1.0000 </r>
       <r>   10.0835    0.0000 </r>
       <r>   12.2632    0.0000 </r>
       <r>   13.1030    0.0000 </r>
       <r>   18.0118    0.0000 </r>
      </set>
      <set comment="kpoint 65">
       <r>    3.1512    1.0000 </r>
       <r>    7.7385    0.8741 </r>
       <r>   11.1916    0.0000 </r>
       <r>   13.1430    0.0000 </r>
       <r>   18.3471    0.0000 </r>
      </set>
      <set comment="kpoint 66">
       <r>    4.8972    1.0000 </r>
       <r>    5.6246    1.0000 </r>
       <r>   10.3862    0.0000 </r>
       <r>   13.5521    0.0000 </r>
       <r>   18.8516    0.0000 </r>
      </set>
      <set comment="kpoint 67">
       <r>    3.7194    1.0000 </r>
       <r>    6.9194    0.9951 </r>
       <r>    9.9299    0.0000 </r>
       <r>   14.2668    0.0000 </r>
       <r>   19.6050    0.0000 </r>
      </set>
      <set comment="kpoint 68">
       <r>    2.1301    1.0000 </r>
       <r>    8.7209   -0.0065 </r>
       <r>   10.2275    0.0000 </r>
       <r>   15.2714    0.0000 </r>
       <r>   20.6107    0.0000 </r>
      </set>
      <set comment="kpoint 69">
       <r>    2.3288    1.0000 </r>
       <r>   10.3743    0.0000 </r>
       <r>   11.4891    0.0000 </r>
       <r>   12.8494    0.0000 </r>
       <r>   16.6458    0.0000 </r>
      </set>
      <set comment="kpoint 70">
       <r>    3.6976    1.0000 </r>
       <r>    9.0878    0.0035 </r>
       <r>    9.5783    0.0001 </r>
       <r>   12.4683    0.0000 </r>
       <r>   17.4349    0.0000 </r>
      </set>
      <set comment="kpoint 71">
       <r>    5.3026    1.0000 </r>
       <r>    7.3024    1.0310 </r>
       <r>    8.3514   -0.0446 </r>
       <r>   12.7441    0.0000 </r>
       <r>   17.9677    0.0000 </r>
      </set>
      <set comment="kpoint 72">
       <r>    5.3813    1.0000 </r>
       <r>    6.7487    0.9971 </r>
       <r>    8.2121    0.0360 </r>
       <r>   13.3561    0.0000 </r>
       <r>   18.6664    0.0000 </r>
      </set>
      <set comment="kpoint 73">
       <r>    4.5020    1.0000 </r>
       <r>    7.4705    1.0514 </r>
       <r>   11.0930    0.0000 </r>
       <r>   12.4574    0.0000 </r>
       <r>   16.2568    0.0000 </r>
      </set>
      <set comment="kpoint 74">
       <r>    5.5683    1.0000 </r>
       <r>    6.9102    0.9951 </r>
       <r>    9.3209    0.0008 </r>
       <r>   12.2701    0.0000 </r>
       <r>   17.2411    0.0000 </r>
      </set>
      <set comment="kpoint 75">
       <r>    1.0927    1.0000 </r>
       <r>   10.3020    0.0000 </r>
       <r>   11.7654    0.0000 </r>
       <r>   14.3713    0.0000 </r>
       <r>   20.5530    0.0000 </r>
      </set>
      <set comment="kpoint 76">
       <r>    2.4891    1.0000 </r>
       <r>    8.3117   -0.0321 </r>
       <r>   11.2845    0.0000 </r>
       <r>   13.1729    0.0000 </r>
       <r>   21.2332    0.0000 </r>
      </set>
      <set comment="kpoint 77">
       <r>    4.1696    1.0000 </r>
       <r>    6.1448    1.0000 </r>
       <r>   11.4850    0.0000 </r>
       <r>   12.3146    0.0000 </r>
       <r>   21.6515    0.0000 </r>
      </set>
      <set comment="kpoint 78">
       <r>    2.9462    1.0000 </r>
       <r>    9.5115    0.0001 </r>
       <r>   10.6699    0.0000 </r>
       <r>   11.3690    0.0000 </r>
       <r>   19.7986    0.0000 </r>
      </set>
      <set comment="kpoint 79">
       <r>    4.5096    1.0000 </r>
       <r>    7.7074    0.9158 </r>
       <r>    9.8927    0.0000 </r>
       <r>   10.8643    0.0000 </r>
       <r>   20.6961    0.0000 </r>
      </set>
      <set comment="kpoint 80">
       <r>    5.7016    1.0000 </r>
       <r>    6.3354    0.9999 </r>
       <r>    9.2983    0.0010 </r>
       <r>   11.2554    0.0000 </r>
       <r>   21.3526    0.0000 </r>
      </set>
      <set comment="kpoint 81">
       <r>    3.9686    1.0000 </r>
       <r>    7.9315    0.5137 </r>
       <r>    9.4006    0.0004 </r>
       <r>   12.0547    0.0000 </r>
       <r>   22.1590    0.0000 </r>
      </set>
      <set comment="kpoint 82">
       <r>    5.1044    1.0000 </r>
       <r>    8.0578    0.2575 </r>
       <r>    9.1163    0.0031 </r>
       <r>   10.6891    0.0000 </r>
       <r>   19.4155    0.0000 </r>
      </set>
      <set comment="kpoint 83">
       <r>    6.3550    0.9999 </r>
       <r>    7.3825    1.0453 </r>
       <r>    7.8767    0.6285 </r>
       <r>   10.5937    0.0000 </r>
       <r>   20.5118    0.0000 </r>
      </set>
      <set comment="kpoint 84">
       <r>    6.1012    1.0000 </r>
       <r>    7.5371    1.0415 </r>
       <r>    8.6857   -0.0112 </r>
       <r>    9.5044    0.0002 </r>
       <r>   22.9474    0.0000 </r>
      </set>
     </set>
     <set comment="spin 2">
      <set comment="kpoint 1">
       <r>   -3.1341    1.0000 </r>
       <r>   20.6208    0.0000 </r>
       <r>   20.6208    0.0000 </r>
       <r>   20.6208    0.0000 </r>
       <r>   22.1131    0.0000 </r>
      </set>
      <set comment="kpoint 2">
       <r>   -2.9732    1.0000 </r>
       <r>   18.2625    0.0000 </r>
       <r>   20.1120    0.0000 </r>
       <r>   20.1120    0.0000 </r>
       <r>   22.7714    0.0000 </r>
      </set>
      <set comment="kpoint 3">
       <r>   -2.4918    1.0000 </r>
       <r>   15.1545    0.0000 </r>
       <r>   19.4447    0.0000 </r>
       <r>   19.4447    0.0000 </r>
       <r>   22.6316    0.0000 </r>
      </set>
      <set comment="kpoint 4">
       <r>   -1.6930    1.0000 </r>
       <r>   12.1586    0.0000 </r>
       <r>   18.9854    0.0000 </r>
       <r>   18.9854    0.0000 </r>
       <r>   22.3980    0.0000 </r>
      </set>
      <set comment="kpoint 5">
       <r>   -0.5822    1.0000 </r>
       <r>    9.3722    0.0005 </r>
       <r>   18.7803    0.0000 </r>
       <r>   18.7805    0.0000 </r>
       <r>   21.7676    0.0000 </r>
      </set>
      <set comment="kpoint 6">
       <r>    0.8335    1.0000 </r>
       <r>    6.8308    0.9959 </r>
       <r>   18.8406    0.0000 </r>
       <r>   18.8412    0.0000 </r>
       <r>   20.6402    0.0000 </r>
      </set>
      <set comment="kpoint 7">
       <r>    2.5434    1.0000 </r>
       <r>    4.5561    1.0000 </r>
       <r>   19.1696    0.0000 </r>
       <r>   19.1741    0.0000 </r>
       <r>   19.7719    0.0000 </r>
      </set>
      <set comment="kpoint 8">
       <r>   -2.9196    1.0000 </r>
       <r>   19.0625    0.0000 </r>
       <r>   19.0625    0.0000 </r>
       <r>   19.3530    0.0000 </r>
       <r>   21.6182    0.0000 </r>
      </set>
      <set comment="kpoint 9">
       <r>   -2.5452    1.0000 </r>
       <r>   16.2579    0.0000 </r>
       <r>   18.2305    0.0000 </r>
       <r>   19.4493    0.0000 </r>
       <r>   21.9398    0.0000 </r>
      </set>
      <set comment="kpoint 10">
       <r>   -1.8524    1.0000 </r>
       <r>   13.2453    0.0000 </r>
       <r>   17.6567    0.0000 </r>
       <r>   19.1235    0.0000 </r>
       <r>   21.3428    0.0000 </r>
      </set>
      <set comment="kpoint 11">
       <r>   -0.8459    1.0000 </r>
       <r>   10.3980    0.0000 </r>
       <r>   17.3504    0.0000 </r>
       <r>   18.8807    0.0000 </r>
       <r>   19.8568    0.0000 </r>
      </set>
      <set comment="kpoint 12">
       <r>    0.4677    1.0000 </r>
       <r>    7.7792    0.8110 </r>
       <r>   17.3148    0.0000 </r>
       <r>   18.5785    0.0000 </r>
       <r>   18.9007    0.0000 </r>
      </set>
      <set comment="kpoint 13">
       <r>    2.0796    1.0000 </r>
       <r>    5.4172    1.0000 </r>
       <r>   17.5514    0.0000 </r>
       <r>   17.6063    0.0000 </r>
       <r>   19.1519    0.0000 </r>
      </set>
      <set comment="kpoint 14">
       <r>    3.3035    1.0000 </r>
       <r>    4.0052    1.0000 </r>
       <r>   16.8983    0.0000 </r>
       <r>   18.0612    0.0000 </r>
       <r>   19.6931    0.0000 </r>
      </set>
      <set comment="kpoint 15">
       <r>    1.5100    1.0000 </r>
       <r>    6.1739    1.0000 </r>
       <r>   16.4703    0.0000 </r>
       <r>   18.8428    0.0000 </r>
       <r>   20.5099    0.0000 </r>
      </set>
      <set comment="kpoint 16">
       <r>   -0.0038    1.0000 </r>
       <r>    8.6224   -0.0218 </r>
       <r>   16.3249    0.0000 </r>
       <r>   19.8771    0.0000 </r>
       <r>   20.5167    0.0000 </r>
      </set>
      <set comment="kpoint 17">
       <r>   -1.2159    1.0000 </r>
       <r>   11.3163    0.0000 </r>
       <r>   16.4641    0.0000 </r>
       <r>   20.5551    0.0000 </r>
       <r>   21.2451    0.0000 </r>
      </set>
      <set comment="kpoint 18">
       <r>   -2.1185    1.0000 </r>
       <r>   14.2062    0.0000 </r>
       <r>   16.8998    0.0000 </r>
       <r>   20.8874    0.0000 </r>
       <r>   22.6175    0.0000 </r>
      </set>
      <set comment="kpoint 19">
       <r>   -2.7056    1.0000 </r>
       <r>   16.9661    0.0000 </r>
       <r>   17.8588    0.0000 </r>
       <r>   21.0378    0.0000 </r>
       <r>   22.3631    0.0000 </r>
      </set>
      <set comment="kpoint 20">
       <r>   -2.2783    1.0000 </r>
       <r>   17.2814    0.0000 </r>
       <r>   17.2814    0.0000 </r>
       <r>   17.4619    0.0000 </r>
       <r>   20.5135    0.0000 </r>
      </set>
      <set comment="kpoint 21">
       <r>   -1.6930    1.0000 </r>
       <r>   14.5332    0.0000 </r>
       <r>   16.6041    0.0000 </r>
       <r>   17.8754    0.0000 </r>
       <r>   19.2377    0.0000 </r>
      </set>
      <set comment="kpoint 22">
       <r>   -0.7932    1.0000 </r>
       <r>   11.6813    0.0000 </r>
       <r>   16.2000    0.0000 </r>
       <r>   17.5061    0.0000 </r>
       <r>   17.8750    0.0000 </r>
      </set>
      <set comment="kpoint 23">
       <r>    0.4150    1.0000 </r>
       <r>    9.0083    0.0046 </r>
       <r>   16.0694    0.0000 </r>
       <r>   16.1799    0.0000 </r>
       <r>   17.7810    0.0000 </r>
      </set>
      <set comment="kpoint 24">
       <r>    1.9234    1.0000 </r>
       <r>    6.5725    0.9991 </r>
       <r>   15.0818    0.0000 </r>
       <r>   16.2127    0.0000 </r>
       <r>   17.9735    0.0000 </r>
      </set>
      <set comment="kpoint 25">
       <r>    3.7068    1.0000 </r>
       <r>    4.4152    1.0000 </r>
       <r>   14.2658    0.0000 </r>
       <r>   16.6303    0.0000 </r>
       <r>   18.4467    0.0000 </r>
      </set>
      <set comment="kpoint 26">
       <r>    2.4919    1.0000 </r>
       <r>    5.8125    1.0000 </r>
       <r>   13.7414    0.0000 </r>
       <r>   17.3224    0.0000 </r>
       <r>   19.2038    0.0000 </r>
      </set>
      <set comment="kpoint 27">
       <r>    0.8854    1.0000 </r>
       <r>    8.1494    0.1093 </r>
       <r>   13.5200    0.0000 </r>
       <r>   18.2870    0.0000 </r>
       <r>   20.2350    0.0000 </r>
      </set>
      <set comment="kpoint 28">
       <r>   -0.4243    1.0000 </r>
       <r>   10.6937    0.0000 </r>
       <r>   13.6386    0.0000 </r>
       <r>   19.5157    0.0000 </r>
       <r>   21.5173    0.0000 </r>
      </set>
      <set comment="kpoint 29">
       <r>   -1.4278    1.0000 </r>
       <r>   13.0656    0.0000 </r>
       <r>   14.4433    0.0000 </r>
       <r>   20.9757    0.0000 </r>
       <r>   22.2850    0.0000 </r>
      </set>
      <set comment="kpoint 30">
       <r>   -1.2160    1.0000 </r>
       <r>   15.8284    0.0000 </r>
       <r>   15.8284    0.0000 </r>
       <r>   15.8293    0.0000 </r>
       <r>   17.2063    0.0000 </r>
      </set>
      <set comment="kpoint 31">
       <r>   -0.4246    1.0000 </r>
       <r>   13.1660    0.0000 </r>
       <r>   15.3287    0.0000 </r>
       <r>   15.4253    0.0000 </r>
       <r>   16.5538    0.0000 </r>
      </set>
      <set comment="kpoint 32">
       <r>    0.6755    1.0000 </r>
       <r>   10.4935    0.0000 </r>
       <r>   13.9343    0.0000 </r>
       <r>   15.1040    0.0000 </r>
       <r>   16.6700    0.0000 </r>
      </set>
      <set comment="kpoint 33">
       <r>    2.0762    1.0000 </r>
       <r>    8.0069    0.3565 </r>
       <r>   12.7214    0.0000 </r>
       <r>   15.1539    0.0000 </r>
       <r>   16.8304    0.0000 </r>
      </set>
      <set comment="kpoint 34">
       <r>    3.7624    1.0000 </r>
       <r>    5.7667    1.0000 </r>
       <r>   11.8000    0.0000 </r>
       <r>   15.4786    0.0000 </r>
       <r>   17.2330    0.0000 </r>
      </set>
      <set comment="kpoint 35">
       <r>    3.7703    1.0000 </r>
       <r>    5.7430    1.0000 </r>
       <r>   11.1866    0.0000 </r>
       <r>   16.0784    0.0000 </r>
       <r>   17.9154    0.0000 </r>
      </set>
      <set comment="kpoint 36">
       <r>    2.0791    1.0000 </r>
       <r>    7.9166    0.5452 </r>
       <r>   10.9288    0.0000 </r>
       <r>   16.9535    0.0000 </r>
       <r>   18.8771    0.0000 </r>
      </set>
      <set comment="kpoint 37">
       <r>    0.6766    1.0000 </r>
       <r>    9.9061    0.0000 </r>
       <r>   11.4181    0.0000 </r>
       <r>   18.1019    0.0000 </r>
       <r>   20.0968    0.0000 </r>
      </set>
      <set comment="kpoint 38">
       <r>    0.2572    1.0000 </r>
       <r>   13.5342    0.0000 </r>
       <r>   14.5772    0.0000 </r>
       <r>   14.7348    0.0000 </r>
       <r>   14.7348    0.0000 </r>
      </set>
      <set comment="kpoint 39">
       <r>    1.2467    1.0000 </r>
       <r>   11.8647    0.0000 </r>
       <r>   12.2444    0.0000 </r>
       <r>   14.4167    0.0000 </r>
       <r>   15.5526    0.0000 </r>
      </set>
      <set comment="kpoint 40">
       <r>    2.5349    1.0000 </r>
       <r>    9.6698    0.0000 </r>
       <r>   10.6384    0.0000 </r>
       <r>   14.3736    0.0000 </r>
       <r>   15.8675    0.0000 </r>
      </set>
      <set comment="kpoint 41">
       <r>    4.1038    1.0000 </r>
       <r>    7.3956    1.0471 </r>
       <r>    9.6030    0.0001 </r>
       <r>   14.6052    0.0000 </r>
       <r>   16.2338    0.0000 </r>
      </set>
      <set comment="kpoint 42">
       <r>    5.3203    1.0000 </r>
       <r>    5.9384    1.0000 </r>
       <r>    8.9330    0.0049 </r>
       <r>   15.1119    0.0000 </r>
       <r>   16.8439    0.0000 </r>
      </set>
      <set comment="kpoint 43">
       <r>    3.5681    1.0000 </r>
       <r>    7.5369    1.0415 </r>
       <r>    9.0263    0.0044 </r>
       <r>   15.8946    0.0000 </r>
       <r>   17.7290    0.0000 </r>
      </set>
      <set comment="kpoint 44">
       <r>    2.1225    1.0000 </r>
       <r>   10.1899    0.0000 </r>
       <r>   13.7349    0.0000 </r>
       <r>   14.0057    0.0000 </r>
       <r>   14.0057    0.0000 </r>
      </set>
      <set comment="kpoint 45">
       <r>    3.2895    1.0000 </r>
       <r>    8.7585   -0.0025 </r>
       <r>   11.6066    0.0000 </r>
       <r>   13.8699    0.0000 </r>
       <r>   14.9490    0.0000 </r>
      </set>
      <set comment="kpoint 46">
       <r>    4.7029    1.0000 </r>
       <r>    7.6607    0.9671 </r>
       <r>    9.3425    0.0007 </r>
       <r>   14.0086    0.0000 </r>
       <r>   15.4595    0.0000 </r>
      </set>
      <set comment="kpoint 47">
       <r>    5.9615    1.0000 </r>
       <r>    7.0099    0.9961 </r>
       <r>    7.5610    1.0335 </r>
       <r>   14.4224    0.0000 </r>
       <r>   16.0315    0.0000 </r>
      </set>
      <set comment="kpoint 48">
       <r>    4.3007    1.0000 </r>
       <r>    7.2695    1.0247 </r>
       <r>   13.3121    0.0000 </r>
       <r>   13.6415    0.0000 </r>
       <r>   13.6415    0.0000 </r>
      </set>
      <set comment="kpoint 49">
       <r>    5.1754    1.0000 </r>
       <r>    6.5036    0.9995 </r>
       <r>   11.4071    0.0000 </r>
       <r>   13.6877    0.0000 </r>
       <r>   14.7469    0.0000 </r>
      </set>
      <set comment="kpoint 50">
       <r>   -2.0652    1.0000 </r>
       <r>   15.0315    0.0000 </r>
       <r>   16.1512    0.0000 </r>
       <r>   20.0199    0.0000 </r>
       <r>   21.1664    0.0000 </r>
      </set>
      <set comment="kpoint 51">
       <r>   -1.2689    1.0000 </r>
       <r>   12.4574    0.0000 </r>
       <r>   15.2797    0.0000 </r>
       <r>   19.5420    0.0000 </r>
       <r>   20.1691    0.0000 </r>
      </set>
      <set comment="kpoint 52">
       <r>   -0.1614    1.0000 </r>
       <r>    9.7382    0.0000 </r>
       <r>   15.0008    0.0000 </r>
       <r>   18.0988    0.0000 </r>
       <r>   20.0874    0.0000 </r>
      </set>
      <set comment="kpoint 53">
       <r>    1.2498    1.0000 </r>
       <r>    7.2214    1.0161 </r>
       <r>   15.0354    0.0000 </r>
       <r>   16.9219    0.0000 </r>
       <r>   20.1896    0.0000 </r>
      </set>
      <set comment="kpoint 54">
       <r>    2.9531    1.0000 </r>
       <r>    4.9602    1.0000 </r>
       <r>   15.3444    0.0000 </r>
       <r>   16.0412    0.0000 </r>
       <r>   20.5273    0.0000 </r>
      </set>
      <set comment="kpoint 55">
       <r>   -1.0045    1.0000 </r>
       <r>   13.4633    0.0000 </r>
       <r>   14.6801    0.0000 </r>
       <r>   17.4148    0.0000 </r>
       <r>   18.5677    0.0000 </r>
      </set>
      <set comment="kpoint 56">
       <r>   -0.0042    1.0000 </r>
       <r>   11.0792    0.0000 </r>
       <r>   13.9496    0.0000 </r>
       <r>   15.8519    0.0000 </r>
       <r>   18.9514    0.0000 </r>
      </set>
      <set comment="kpoint 57">
       <r>    1.3008    1.0000 </r>
       <r>    8.5425   -0.0369 </r>
       <r>   13.8129    0.0000 </r>
       <r>   14.5909    0.0000 </r>
       <r>   19.0846    0.0000 </r>
      </set>
      <set comment="kpoint 58">
       <r>    2.9004    1.0000 </r>
       <r>    6.2143    1.0000 </r>
       <r>   13.4410    0.0000 </r>
       <r>   14.1743    0.0000 </r>
       <r>   19.3910    0.0000 </r>
      </set>
      <set comment="kpoint 59">
       <r>    4.1120    1.0000 </r>
       <r>    4.8161    1.0000 </r>
       <r>   12.7617    0.0000 </r>
       <r>   14.6560    0.0000 </r>
       <r>   19.9480    0.0000 </r>
      </set>
      <set comment="kpoint 60">
       <r>    2.3367    1.0000 </r>
       <r>    6.9472    0.9951 </r>
       <r>   12.3544    0.0000 </r>
       <r>   15.4615    0.0000 </r>
       <r>   20.7617    0.0000 </r>
      </set>
      <set comment="kpoint 61">
       <r>    0.8329    1.0000 </r>
       <r>    9.3076    0.0009 </r>
       <r>   12.2910    0.0000 </r>
       <r>   16.5547    0.0000 </r>
       <r>   21.7678    0.0000 </r>
      </set>
      <set comment="kpoint 62">
       <r>   -0.3717    1.0000 </r>
       <r>   11.4872    0.0000 </r>
       <r>   12.9443    0.0000 </r>
       <r>   17.9221    0.0000 </r>
       <r>   21.6206    0.0000 </r>
      </set>
      <set comment="kpoint 63">
       <r>    0.4665    1.0000 </r>
       <r>   12.2802    0.0000 </r>
       <r>   13.4000    0.0000 </r>
       <r>   13.9191    0.0000 </r>
       <r>   17.4206    0.0000 </r>
      </set>
      <set comment="kpoint 64">
       <r>    1.6620    1.0000 </r>
       <r>   10.0835    0.0000 </r>
       <r>   12.2632    0.0000 </r>
       <r>   13.1030    0.0000 </r>
       <r>   18.0118    0.0000 </r>
      </set>
      <set comment="kpoint 65">
       <r>    3.1512    1.0000 </r>
       <r>    7.7385    0.8741 </r>
       <r>   11.1916    0.0000 </r>
       <r>   13.1430    0.0000 </r>
       <r>   18.3471    0.0000 </r>
      </set>
      <set comment="kpoint 66">
       <r>    4.8972    1.0000 </r>
       <r>    5.6246    1.0000 </r>
       <r>   10.3862    0.0000 </r>
       <r>   13.5521    0.0000 </r>
       <r>   18.8516    0.0000 </r>
      </set>
      <set comment="kpoint 67">
       <r>    3.7194    1.0000 </r>
       <r>    6.9194    0.9951 </r>
       <r>    9.9299    0.0000 </r>
       <r>   14.2668    0.0000 </r>
       <r>   19.6050    0.0000 </r>
      </set>
      <set comment="kpoint 68">
       <r>    2.1301    1.0000 </r>
       <r>    8.7209   -0.0065 </r>
       <r>   10.2275    0.0000 </r>
       <r>   15.2714    0.0000 </r>
       <r>   20.6107    0.0000 </r>
      </set>
      <set comment="kpoint 69">
       <r>    2.3288    1.0000 </r>
       <r>   10.3743    0.0000 </r>
       <r>   11.4891    0.0000 </r>
       <r>   12.8494    0.0000 </r>
       <r>   16.6458    0.0000 </r>
      </set>
      <set comment="kpoint 70">
       <r>    3.6976    1.0000 </r>
       <r>    9.0878    0.0035 </r>
       <r>    9.5783    0.0001 </r>
       <r>   12.4683    0.0000 </r>
       <r>   17.4349    0.0000 </r>
      </set>
      <set comment="kpoint 71">
       <r>    5.3026    1.0000 </r>
       <r>    7.3024    1.0310 </r>
       <r>    8.3514   -0.0446 </r>
       <r>   12.7441    0.0000 </r>
       <r>   17.9677    0.0000 </r>
      </set>
      <set comment="kpoint 72">
       <r>    5.3813    1.0000 </r>
       <r>    6.7487    0.9971 </r>
       <r>    8.2121    0.0360 </r>
       <r>   13.3561    0.0000 </r>
       <r>   18.6664    0.0000 </r>
      </set>
      <set comment="kpoint 73">
       <r>    4.5020    1.0000 </r>
       <r>    7.4705    1.0514 </r>
       <r>   11.0930    0.0000 </r>
       <r>   12.4574    0.0000 </r>
       <r>   16.2568    0.0000 </r>
      </set>
      <set comment="kpoint 74">
       <r>    5.5683    1.0000 </r>
       <r>    6.9102    0.9951 </r>
       <r>    9.3209    0.0008 </r>
       <r>   12.2701    0.0000 </r>
       <r>   17.2411    0.0000 </r>
      </set>
      <set comment="kpoint 75">
       <r>    1.0927    1.0000 </r>
       <r>   10.3020    0.0000 </r>
       <r>   11.7654    0.0000 </r>
       <r>   14.3713    0.0000 </r>
       <r>   20.5530    0.0000 </r>
      </set>
      <set comment="kpoint 76">
       <r>    2.4891    1.0000 </r>
       <r>    8.3117   -0.0321 </r>
       <r>   11.2845    0.0000 </r>
       <r>   13.1729    0.0000 </r>
       <r>   21.2332    0.0000 </r>
      </set>
      <set comment="kpoint 77">
       <r>    4.1696    1.0000 </r>
       <r>    6.1448    1.0000 </r>
       <r>   11.4850    0.0000 </r>
       <r>   12.3146    0.0000 </r>
       <r>   21.6506    0.0000 </r>
      </set>
      <set comment="kpoint 78">
       <r>    2.9462    1.0000 </r>
       <r>    9.5115    0.0001 </r>
       <r>   10.6699    0.0000 </r>
       <r>   11.3690    0.0000 </r>
       <r>   19.7985    0.0000 </r>
      </set>
      <set comment="kpoint 79">
       <r>    4.5096    1.0000 </r>
       <r>    7.7074    0.9158 </r>
       <r>    9.8927    0.0000 </r>
       <r>   10.8643    0.0000 </r>
       <r>   20.6961    0.0000 </r>
      </set>
      <set comment="kpoint 80">
       <r>    5.7016    1.0000 </r>
       <r>    6.3354    0.9999 </r>
       <r>    9.2983    0.0010 </r>
       <r>   11.2554    0.0000 </r>
       <r>   21.3526    0.0000 </r>
      </set>
      <set comment="kpoint 81">
       <r>    3.9686    1.0000 </r>
       <r>    7.9315    0.5137 </r>
       <r>    9.4006    0.0004 </r>
       <r>   12.0547    0.0000 </r>
       <r>   22.1590    0.0000 </r>
      </set>
      <set comment="kpoint 82">
       <r>    5.1044    1.0000 </r>
       <r>    8.0578    0.2575 </r>
       <r>    9.1163    0.0031 </r>
       <r>   10.6891    0.0000 </r>
       <r>   19.4155    0.0000 </r>
      </set>
      <set comment="kpoint 83">
       <r>    6.3550    0.9999 </r>
       <r>    7.3825    1.0453 </r>
       <r>    7.8767    0.6285 </r>
       <r>   10.5937    0.0000 </r>
       <r>   20.5118    0.0000 </r>
      </set>
      <set comment="kpoint 84">
       <r>    6.1012    1.0000 </r>
       <r>    7.5371    1.0415 </r>
       <r>    8.6857   -0.0112 </r>
       <r>    9.5044    0.0002 </r>
       <r>   22.9401    0.0000 </r>
      </set>
     </set>
    </set>
   </array>
  </eigenvalues>
  <separator name="orbital magnetization" >
   <v name="MAGDIPOLOUT">      0.00000000      0.00000000      0.00000000</v>
  </separator>
  <dos>
   <i name="efermi">      7.93802192 </i>
   <total>
    <array>
     <dimension dim="1">gridpoints</dimension>
     <dimension dim="2">spin</dimension>
     <field>energy</field>
     <field>total</field>
     <field>integrated</field>
     <set>
      <set comment="spin 1">
       <r>    -8.1341     0.0000     0.0000 </r>
       <r>    -8.0138     0.0000     0.0000 </r>
       <r>    -7.8935     0.0000     0.0000 </r>
       <r>    -7.7732     0.0000     0.0000 </r>
       <r>    -7.6530     0.0000     0.0000 </r>
       <r>    -7.5327     0.0000     0.0000 </r>
       <r>    -7.4124     0.0000     0.0000 </r>
       <r>    -7.2922     0.0000     0.0000 </r>
       <r>    -7.1719     0.0000     0.0000 </r>
       <r>    -7.0516     0.0000     0.0000 </r>
       <r>    -6.9313     0.0000     0.0000 </r>
       <r>    -6.8111     0.0000     0.0000 </r>
       <r>    -6.6908     0.0000     0.0000 </r>
       <r>    -6.5705     0.0000     0.0000 </r>
       <r>    -6.4503     0.0000     0.0000 </r>
       <r>    -6.3300     0.0000     0.0000 </r>
       <r>    -6.2097     0.0000     0.0000 </r>
       <r>    -6.0894     0.0000     0.0000 </r>
       <r>    -5.9692     0.0000     0.0000 </r>
       <r>    -5.8489     0.0000     0.0000 </r>
       <r>    -5.7286     0.0000     0.0000 </r>
       <r>    -5.6084     0.0000     0.0000 </r>
       <r>    -5.4881     0.0000     0.0000 </r>
       <r>    -5.3678     0.0000     0.0000 </r>
       <r>    -5.2475     0.0000     0.0000 </r>
       <r>    -5.1273     0.0000     0.0000 </r>
       <r>    -5.0070     0.0000     0.0000 </r>
       <r>    -4.8867     0.0000     0.0000 </r>
       <r>    -4.7665     0.0000     0.0000 </r>
       <r>    -4.6462     0.0000     0.0000 </r>
       <r>    -4.5259     0.0000     0.0000 </r>
       <r>    -4.4056     0.0000     0.0000 </r>
       <r>    -4.2854     0.0001     0.0000 </r>
       <r>    -4.1651     0.0001     0.0000 </r>
       <r>    -4.0448     0.0001     0.0000 </r>
       <r>    -3.9246     0.0001     0.0001 </r>
       <r>    -3.8043    -0.0001     0.0000 </r>
       <r>    -3.6840    -0.0006    -0.0000 </r>
       <r>    -3.5637    -0.0011    -0.0002 </r>
       <r>    -3.4435    -0.0012    -0.0003 </r>
       <r>    -3.3232    -0.0002    -0.0003 </r>
       <r>    -3.2029     0.0030     0.0000 </r>
       <r>    -3.0827     0.0086     0.0011 </r>
       <r>    -2.9624     0.0165     0.0031 </r>
       <r>    -2.8421     0.0254     0.0061 </r>
       <r>    -2.7218     0.0336     0.0102 </r>
       <r>    -2.6016     0.0396     0.0149 </r>
       <r>    -2.4813     0.0435     0.0202 </r>
       <r>    -2.3610     0.0472     0.0258 </r>
       <r>    -2.2408     0.0526     0.0321 </r>
       <r>    -2.1205     0.0596     0.0393 </r>
       <r>    -2.0002     0.0653     0.0472 </r>
       <r>    -1.8799     0.0666     0.0552 </r>
       <r>    -1.7597     0.0637     0.0628 </r>
       <r>    -1.6394     0.0611     0.0702 </r>
       <r>    -1.5191     0.0640     0.0779 </r>
       <r>    -1.3989     0.0743     0.0868 </r>
       <r>    -1.2786     0.0874     0.0973 </r>
       <r>    -1.1583     0.0960     0.1089 </r>
       <r>    -1.0380     0.0951     0.1203 </r>
       <r>    -0.9178     0.0860     0.1307 </r>
       <r>    -0.7975     0.0752     0.1397 </r>
       <r>    -0.6772     0.0704     0.1482 </r>
       <r>    -0.5570     0.0764     0.1574 </r>
       <r>    -0.4367     0.0924     0.1685 </r>
       <r>    -0.3164     0.1120     0.1819 </r>
       <r>    -0.1961     0.1259     0.1971 </r>
       <r>    -0.0759     0.1269     0.2123 </r>
       <r>     0.0444     0.1148     0.2261 </r>
       <r>     0.1647     0.0973     0.2379 </r>
       <r>     0.2849     0.0853     0.2481 </r>
       <r>     0.4052     0.0854     0.2584 </r>
       <r>     0.5255     0.0962     0.2700 </r>
       <r>     0.6458     0.1102     0.2832 </r>
       <r>     0.7660     0.1214     0.2978 </r>
       <r>     0.8863     0.1291     0.3133 </r>
       <r>     1.0066     0.1365     0.3298 </r>
       <r>     1.1268     0.1449     0.3472 </r>
       <r>     1.2471     0.1499     0.3652 </r>
       <r>     1.3674     0.1447     0.3826 </r>
       <r>     1.4877     0.1271     0.3979 </r>
       <r>     1.6079     0.1040     0.4104 </r>
       <r>     1.7282     0.0878     0.4210 </r>
       <r>     1.8485     0.0890     0.4317 </r>
       <r>     1.9687     0.1095     0.4448 </r>
       <r>     2.0890     0.1412     0.4618 </r>
       <r>     2.2093     0.1707     0.4824 </r>
       <r>     2.3296     0.1856     0.5047 </r>
       <r>     2.4498     0.1814     0.5265 </r>
       <r>     2.5701     0.1640     0.5462 </r>
       <r>     2.6904     0.1466     0.5639 </r>
       <r>     2.8106     0.1394     0.5806 </r>
       <r>     2.9309     0.1418     0.5977 </r>
       <r>     3.0512     0.1434     0.6149 </r>
       <r>     3.1715     0.1353     0.6312 </r>
       <r>     3.2917     0.1211     0.6458 </r>
       <r>     3.4120     0.1145     0.6595 </r>
       <r>     3.5323     0.1275     0.6749 </r>
       <r>     3.6525     0.1583     0.6939 </r>
       <r>     3.7728     0.1922     0.7170 </r>
       <r>     3.8931     0.2133     0.7427 </r>
       <r>     4.0134     0.2145     0.7685 </r>
       <r>     4.1336     0.1993     0.7924 </r>
       <r>     4.2539     0.1772     0.8138 </r>
       <r>     4.3742     0.1583     0.8328 </r>
       <r>     4.4944     0.1493     0.8507 </r>
       <r>     4.6147     0.1518     0.8690 </r>
       <r>     4.7350     0.1615     0.8884 </r>
       <r>     4.8553     0.1706     0.9090 </r>
       <r>     4.9755     0.1736     0.9298 </r>
       <r>     5.0958     0.1709     0.9504 </r>
       <r>     5.2161     0.1690     0.9707 </r>
       <r>     5.3363     0.1740     0.9916 </r>
       <r>     5.4566     0.1864     1.0141 </r>
       <r>     5.5769     0.2002     1.0381 </r>
       <r>     5.6972     0.2089     1.0633 </r>
       <r>     5.8174     0.2114     1.0887 </r>
       <r>     5.9377     0.2119     1.1142 </r>
       <r>     6.0580     0.2135     1.1398 </r>
       <r>     6.1782     0.2119     1.1653 </r>
       <r>     6.2985     0.1986     1.1892 </r>
       <r>     6.4188     0.1716     1.2098 </r>
       <r>     6.5391     0.1412     1.2268 </r>
       <r>     6.6593     0.1233     1.2417 </r>
       <r>     6.7796     0.1262     1.2568 </r>
       <r>     6.8999     0.1446     1.2742 </r>
       <r>     7.0201     0.1663     1.2942 </r>
       <r>     7.1404     0.1837     1.3163 </r>
       <r>     7.2607     0.1977     1.3401 </r>
       <r>     7.3810     0.2128     1.3657 </r>
       <r>     7.5012     0.2302     1.3934 </r>
       <r>     7.6215     0.2453     1.4229 </r>
       <r>     7.7418     0.2512     1.4531 </r>
       <r>     7.8620     0.2441     1.4825 </r>
       <r>     7.9823     0.2270     1.5098 </r>
       <r>     8.1026     0.2074     1.5347 </r>
       <r>     8.2229     0.1923     1.5578 </r>
       <r>     8.3431     0.1832     1.5799 </r>
       <r>     8.4634     0.1763     1.6011 </r>
       <r>     8.5837     0.1676     1.6212 </r>
       <r>     8.7039     0.1577     1.6402 </r>
       <r>     8.8242     0.1531     1.6586 </r>
       <r>     8.9445     0.1625     1.6782 </r>
       <r>     9.0648     0.1894     1.7009 </r>
       <r>     9.1850     0.2282     1.7284 </r>
       <r>     9.3053     0.2652     1.7603 </r>
       <r>     9.4256     0.2871     1.7948 </r>
       <r>     9.5458     0.2875     1.8294 </r>
       <r>     9.6661     0.2696     1.8618 </r>
       <r>     9.7864     0.2423     1.8910 </r>
       <r>     9.9067     0.2153     1.9169 </r>
       <r>    10.0269     0.1963     1.9405 </r>
       <r>    10.1472     0.1894     1.9632 </r>
       <r>    10.2675     0.1941     1.9866 </r>
       <r>    10.3877     0.2037     2.0111 </r>
       <r>    10.5080     0.2090     2.0362 </r>
       <r>    10.6283     0.2051     2.0609 </r>
       <r>    10.7486     0.1972     2.0846 </r>
       <r>    10.8688     0.1993     2.1086 </r>
       <r>    10.9891     0.2229     2.1354 </r>
       <r>    11.1094     0.2663     2.1674 </r>
       <r>    11.2296     0.3109     2.2048 </r>
       <r>    11.3499     0.3317     2.2447 </r>
       <r>    11.4702     0.3126     2.2823 </r>
       <r>    11.5905     0.2568     2.3132 </r>
       <r>    11.7107     0.1868     2.3356 </r>
       <r>    11.8310     0.1345     2.3518 </r>
       <r>    11.9513     0.1277     2.3672 </r>
       <r>    12.0715     0.1745     2.3882 </r>
       <r>    12.1918     0.2548     2.4188 </r>
       <r>    12.3121     0.3280     2.4583 </r>
       <r>    12.4324     0.3568     2.5012 </r>
       <r>    12.5526     0.3324     2.5411 </r>
       <r>    12.6729     0.2778     2.5746 </r>
       <r>    12.7932     0.2293     2.6021 </r>
       <r>    12.9134     0.2114     2.6276 </r>
       <r>    13.0337     0.2254     2.6547 </r>
       <r>    13.1540     0.2557     2.6854 </r>
       <r>    13.2743     0.2834     2.7195 </r>
       <r>    13.3945     0.2974     2.7553 </r>
       <r>    13.5148     0.2975     2.7910 </r>
       <r>    13.6351     0.2896     2.8259 </r>
       <r>    13.7553     0.2806     2.8596 </r>
       <r>    13.8756     0.2742     2.8926 </r>
       <r>    13.9959     0.2717     2.9253 </r>
       <r>    14.1162     0.2735     2.9582 </r>
       <r>    14.2364     0.2783     2.9917 </r>
       <r>    14.3567     0.2819     3.0256 </r>
       <r>    14.4770     0.2781     3.0590 </r>
       <r>    14.5972     0.2649     3.0909 </r>
       <r>    14.7175     0.2501     3.1209 </r>
       <r>    14.8378     0.2487     3.1509 </r>
       <r>    14.9581     0.2706     3.1834 </r>
       <r>    15.0783     0.3092     3.2206 </r>
       <r>    15.1986     0.3421     3.2617 </r>
       <r>    15.3189     0.3448     3.3032 </r>
       <r>    15.4391     0.3085     3.3403 </r>
       <r>    15.5594     0.2482     3.3702 </r>
       <r>    15.6797     0.1953     3.3936 </r>
       <r>    15.8000     0.1774     3.4150 </r>
       <r>    15.9202     0.2011     3.4392 </r>
       <r>    16.0405     0.2484     3.4690 </r>
       <r>    16.1608     0.2916     3.5041 </r>
       <r>    16.2810     0.3124     3.5417 </r>
       <r>    16.4013     0.3105     3.5790 </r>
       <r>    16.5216     0.2968     3.6147 </r>
       <r>    16.6419     0.2801     3.6484 </r>
       <r>    16.7621     0.2620     3.6799 </r>
       <r>    16.8824     0.2420     3.7090 </r>
       <r>    17.0027     0.2250     3.7361 </r>
       <r>    17.1229     0.2194     3.7625 </r>
       <r>    17.2432     0.2286     3.7899 </r>
       <r>    17.3635     0.2457     3.8195 </r>
       <r>    17.4838     0.2607     3.8509 </r>
       <r>    17.6040     0.2719     3.8836 </r>
       <r>    17.7243     0.2865     3.9180 </r>
       <r>    17.8446     0.3080     3.9550 </r>
       <r>    17.9648     0.3259     3.9942 </r>
       <r>    18.0851     0.3216     4.0329 </r>
       <r>    18.2054     0.2861     4.0673 </r>
       <r>    18.3257     0.2318     4.0952 </r>
       <r>    18.4459     0.1861     4.1176 </r>
       <r>    18.5662     0.1740     4.1385 </r>
       <r>    18.6865     0.2020     4.1628 </r>
       <r>    18.8067     0.2547     4.1935 </r>
       <r>    18.9270     0.3051     4.2302 </r>
       <r>    19.0473     0.3315     4.2700 </r>
       <r>    19.1676     0.3299     4.3097 </r>
       <r>    19.2878     0.3105     4.3471 </r>
       <r>    19.4081     0.2870     4.3816 </r>
       <r>    19.5284     0.2668     4.4137 </r>
       <r>    19.6486     0.2529     4.4441 </r>
       <r>    19.7689     0.2469     4.4738 </r>
       <r>    19.8892     0.2493     4.5038 </r>
       <r>    20.0095     0.2571     4.5347 </r>
       <r>    20.1297     0.2650     4.5666 </r>
       <r>    20.2500     0.2714     4.5992 </r>
       <r>    20.3703     0.2783     4.6327 </r>
       <r>    20.4905     0.2842     4.6668 </r>
       <r>    20.6108     0.2795     4.7004 </r>
       <r>    20.7311     0.2544     4.7310 </r>
       <r>    20.8514     0.2123     4.7566 </r>
       <r>    20.9716     0.1713     4.7772 </r>
       <r>    21.0919     0.1510     4.7953 </r>
       <r>    21.2122     0.1565     4.8142 </r>
       <r>    21.3324     0.1764     4.8354 </r>
       <r>    21.4527     0.1935     4.8587 </r>
       <r>    21.5730     0.1961     4.8822 </r>
       <r>    21.6933     0.1814     4.9041 </r>
       <r>    21.8135     0.1528     4.9224 </r>
       <r>    21.9338     0.1180     4.9366 </r>
       <r>    22.0541     0.0861     4.9470 </r>
       <r>    22.1743     0.0644     4.9547 </r>
       <r>    22.2946     0.0545     4.9613 </r>
       <r>    22.4149     0.0533     4.9677 </r>
       <r>    22.5352     0.0550     4.9743 </r>
       <r>    22.6554     0.0554     4.9810 </r>
       <r>    22.7757     0.0521     4.9872 </r>
       <r>    22.8960     0.0446     4.9926 </r>
       <r>    23.0162     0.0337     4.9967 </r>
       <r>    23.1365     0.0213     4.9992 </r>
       <r>    23.2568     0.0100     5.0004 </r>
       <r>    23.3771     0.0020     5.0007 </r>
       <r>    23.4973    -0.0018     5.0005 </r>
       <r>    23.6176    -0.0025     5.0002 </r>
       <r>    23.7379    -0.0016     5.0000 </r>
       <r>    23.8581    -0.0005     4.9999 </r>
       <r>    23.9784     0.0001     4.9999 </r>
       <r>    24.0987     0.0003     5.0000 </r>
       <r>    24.2190     0.0002     5.0000 </r>
       <r>    24.3392     0.0001     5.0000 </r>
       <r>    24.4595     0.0001     5.0000 </r>
       <r>    24.5798     0.0000     5.0000 </r>
       <r>    24.7000     0.0000     5.0000 </r>
       <r>    24.8203     0.0000     5.0000 </r>
       <r>    24.9406     0.0000     5.0000 </r>
       <r>    25.0609     0.0000     5.0000 </r>
       <r>    25.1811     0.0000     5.0000 </r>
       <r>    25.3014     0.0000     5.0000 </r>
       <r>    25.4217     0.0000     5.0000 </r>
       <r>    25.5419     0.0000     5.0000 </r>
       <r>    25.6622     0.0000     5.0000 </r>
       <r>    25.7825     0.0000     5.0000 </r>
       <r>    25.9028     0.0000     5.0000 </r>
       <r>    26.0230     0.0000     5.0000 </r>
       <r>    26.1433     0.0000     5.0000 </r>
       <r>    26.2636     0.0000     5.0000 </r>
       <r>    26.3838     0.0000     5.0000 </r>
       <r>    26.5041     0.0000     5.0000 </r>
       <r>    26.6244     0.0000     5.0000 </r>
       <r>    26.7447     0.0000     5.0000 </r>
       <r>    26.8649     0.0000     5.0000 </r>
       <r>    26.9852     0.0000     5.0000 </r>
       <r>    27.1055     0.0000     5.0000 </r>
       <r>    27.2257     0.0000     5.0000 </r>
       <r>    27.3460     0.0000     5.0000 </r>
       <r>    27.4663     0.0000     5.0000 </r>
       <r>    27.5866     0.0000     5.0000 </r>
       <r>    27.7068     0.0000     5.0000 </r>
       <r>    27.8271     0.0000     5.0000 </r>
       <r>    27.9474     0.0000     5.0000 </r>
      </set>
      <set comment="spin 2">
       <r>    -8.1341     0.0000     0.0000 </r>
       <r>    -8.0138     0.0000     0.0000 </r>
       <r>    -7.8935     0.0000     0.0000 </r>
       <r>    -7.7732     0.0000     0.0000 </r>
       <r>    -7.6530     0.0000     0.0000 </r>
       <r>    -7.5327     0.0000     0.0000 </r>
       <r>    -7.4124     0.0000     0.0000 </r>
       <r>    -7.2922     0.0000     0.0000 </r>
       <r>    -7.1719     0.0000     0.0000 </r>
       <r>    -7.0516     0.0000     0.0000 </r>
       <r>    -6.9313     0.0000     0.0000 </r>
       <r>    -6.8111     0.0000     0.0000 </r>
       <r>    -6.6908     0.0000     0.0000 </r>
       <r>    -6.5705     0.0000     0.0000 </r>
       <r>    -6.4503     0.0000     0.0000 </r>
       <r>    -6.3300     0.0000     0.0000 </r>
       <r>    -6.2097     0.0000     0.0000 </r>
       <r>    -6.0894     0.0000     0.0000 </r>
       <r>    -5.9692     0.0000     0.0000 </r>
       <r>    -5.8489     0.0000     0.0000 </r>
       <r>    -5.7286     0.0000     0.0000 </r>
       <r>    -5.6084     0.0000     0.0000 </r>
       <r>    -5.4881     0.0000     0.0000 </r>
       <r>    -5.3678     0.0000     0.0000 </r>
       <r>    -5.2475     0.0000     0.0000 </r>
       <r>    -5.1273     0.0000     0.0000 </r>
       <r>    -5.0070     0.0000     0.0000 </r>
       <r>    -4.8867     0.0000     0.0000 </r>
       <r>    -4.7665     0.0000     0.0000 </r>
       <r>    -4.6462     0.0000     0.0000 </r>
       <r>    -4.5259     0.0000     0.0000 </r>
       <r>    -4.4056     0.0000     0.0000 </r>
       <r>    -4.2854     0.0001     0.0000 </r>
       <r>    -4.1651     0.0001     0.0000 </r>
       <r>    -4.0448     0.0001     0.0000 </r>
       <r>    -3.9246     0.0001     0.0001 </r>
       <r>    -3.8043    -0.0001     0.0000 </r>
       <r>    -3.6840    -0.0006    -0.0000 </r>
       <r>    -3.5637    -0.0011    -0.0002 </r>
       <r>    -3.4435    -0.0012    -0.0003 </r>
       <r>    -3.3232    -0.0002    -0.0003 </r>
       <r>    -3.2029     0.0030     0.0000 </r>
       <r>    -3.0827     0.0086     0.0011 </r>
       <r>    -2.9624     0.0165     0.0031 </r>
       <r>    -2.8421     0.0254     0.0061 </r>
       <r>    -2.7218     0.0336     0.0102 </r>
       <r>    -2.6016     0.0396     0.0149 </r>
       <r>    -2.4813     0.0435     0.0202 </r>
       <r>    -2.3610     0.0472     0.0258 </r>
       <r>    -2.2408     0.0526     0.0321 </r>
       <r>    -2.1205     0.0596     0.0393 </r>
       <r>    -2.0002     0.0653     0.0472 </r>
       <r>    -1.8799     0.0666     0.0552 </r>
       <r>    -1.7597     0.0637     0.0628 </r>
       <r>    -1.6394     0.0611     0.0702 </r>
       <r>    -1.5191     0.0640     0.0779 </r>
       <r>    -1.3989     0.0743     0.0868 </r>
       <r>    -1.2786     0.0874     0.0973 </r>
       <r>    -1.1583     0.0960     0.1089 </r>
       <r>    -1.0380     0.0951     0.1203 </r>
       <r>    -0.9178     0.0860     0.1307 </r>
       <r>    -0.7975     0.0752     0.1397 </r>
       <r>    -0.6772     0.0704     0.1482 </r>
       <r>    -0.5570     0.0764     0.1574 </r>
       <r>    -0.4367     0.0924     0.1685 </r>
       <r>    -0.3164     0.1120     0.1819 </r>
       <r>    -0.1961     0.1259     0.1971 </r>
       <r>    -0.0759     0.1269     0.2123 </r>
       <r>     0.0444     0.1148     0.2261 </r>
       <r>     0.1647     0.0973     0.2379 </r>
       <r>     0.2849     0.0853     0.2481 </r>
       <r>     0.4052     0.0854     0.2584 </r>
       <r>     0.5255     0.0962     0.2700 </r>
       <r>     0.6458     0.1102     0.2832 </r>
       <r>     0.7660     0.1214     0.2978 </r>
       <r>     0.8863     0.1291     0.3133 </r>
       <r>     1.0066     0.1365     0.3298 </r>
       <r>     1.1268     0.1449     0.3472 </r>
       <r>     1.2471     0.1499     0.3652 </r>
       <r>     1.3674     0.1447     0.3826 </r>
       <r>     1.4877     0.1271     0.3979 </r>
       <r>     1.6079     0.1040     0.4104 </r>
       <r>     1.7282     0.0878     0.4210 </r>
       <r>     1.8485     0.0890     0.4317 </r>
       <r>     1.9687     0.1095     0.4448 </r>
       <r>     2.0890     0.1412     0.4618 </r>
       <r>     2.2093     0.1707     0.4824 </r>
       <r>     2.3296     0.1856     0.5047 </r>
       <r>     2.4498     0.1814     0.5265 </r>
       <r>     2.5701     0.1640     0.5462 </r>
       <r>     2.6904     0.1466     0.5639 </r>
       <r>     2.8106     0.1394     0.5806 </r>
       <r>     2.9309     0.1418     0.5977 </r>
       <r>     3.0512     0.1434     0.6149 </r>
       <r>     3.1715     0.1353     0.6312 </r>
       <r>     3.2917     0.1211     0.6458 </r>
       <r>     3.4120     0.1145     0.6595 </r>
       <r>     3.5323     0.1275     0.6749 </r>
       <r>     3.6525     0.1583     0.6939 </r>
       <r>     3.7728     0.1922     0.7170 </r>
       <r>     3.8931     0.2133     0.7427 </r>
       <r>     4.0134     0.2145     0.7685 </r>
       <r>     4.1336     0.1993     0.7924 </r>
       <r>     4.2539     0.1772     0.8138 </r>
       <r>     4.3742     0.1583     0.8328 </r>
       <r>     4.4944     0.1493     0.8507 </r>
       <r>     4.6147     0.1518     0.8690 </r>
       <r>     4.7350     0.1615     0.8884 </r>
       <r>     4.8553     0.1706     0.9090 </r>
       <r>     4.9755     0.1736     0.9298 </r>
       <r>     5.0958     0.1709     0.9504 </r>
       <r>     5.2161     0.1690     0.9707 </r>
       <r>     5.3363     0.1740     0.9916 </r>
       <r>     5.4566     0.1864     1.0141 </r>
       <r>     5.5769     0.2002     1.0381 </r>
       <r>     5.6972     0.2089     1.0633 </r>
       <r>     5.8174     0.2114     1.0887 </r>
       <r>     5.9377     0.2119     1.1142 </r>
       <r>     6.0580     0.2135     1.1398 </r>
       <r>     6.1782     0.2119     1.1653 </r>
       <r>     6.2985     0.1986     1.1892 </r>
       <r>     6.4188     0.1716     1.2098 </r>
       <r>     6.5391     0.1412     1.2268 </r>
       <r>     6.6593     0.1233     1.2417 </r>
       <r>     6.7796     0.1262     1.2568 </r>
       <r>     6.8999     0.1446     1.2742 </r>
       <r>     7.0201     0.1663     1.2942 </r>
       <r>     7.1404     0.1837     1.3163 </r>
       <r>     7.2607     0.1977     1.3401 </r>
       <r>     7.3810     0.2128     1.3657 </r>
       <r>     7.5012     0.2302     1.3934 </r>
       <r>     7.6215     0.2453     1.4229 </r>
       <r>     7.7418     0.2512     1.4531 </r>
       <r>     7.8620     0.2441     1.4825 </r>
       <r>     7.9823     0.2270     1.5098 </r>
       <r>     8.1026     0.2074     1.5347 </r>
       <r>     8.2229     0.1923     1.5578 </r>
       <r>     8.3431     0.1832     1.5799 </r>
       <r>     8.4634     0.1763     1.6011 </r>
       <r>     8.5837     0.1676     1.6212 </r>
       <r>     8.7039     0.1577     1.6402 </r>
       <r>     8.8242     0.1531     1.6586 </r>
       <r>     8.9445     0.1625     1.6782 </r>
       <r>     9.0648     0.1894     1.7009 </r>
       <r>     9.1850     0.2282     1.7284 </r>
       <r>     9.3053     0.2652     1.7603 </r>
       <r>     9.4256     0.2871     1.7948 </r>
       <r>     9.5458     0.2875     1.8294 </r>
       <r>     9.6661     0.2696     1.8618 </r>
       <r>     9.7864     0.2423     1.8910 </r>
       <r>     9.9067     0.2153     1.9169 </r>
       <r>    10.0269     0.1963     1.9405 </r>
       <r>    10.1472     0.1894     1.9632 </r>
       <r>    10.2675     0.1941     1.9866 </r>
       <r>    10.3877     0.2037     2.0111 </r>
       <r>    10.5080     0.2090     2.0362 </r>
       <r>    10.6283     0.2051     2.0609 </r>
       <r>    10.7486     0.1972     2.0846 </r>
       <r>    10.8688     0.1993     2.1086 </r>
       <r>    10.9891     0.2229     2.1354 </r>
       <r>    11.1094     0.2663     2.1674 </r>
       <r>    11.2296     0.3109     2.2048 </r>
       <r>    11.3499     0.3317     2.2447 </r>
       <r>    11.4702     0.3126     2.2823 </r>
       <r>    11.5905     0.2568     2.3132 </r>
       <r>    11.7107     0.1868     2.3356 </r>
       <r>    11.8310     0.1345     2.3518 </r>
       <r>    11.9513     0.1277     2.3672 </r>
       <r>    12.0715     0.1745     2.3882 </r>
       <r>    12.1918     0.2548     2.4188 </r>
       <r>    12.3121     0.3280     2.4583 </r>
       <r>    12.4324     0.3568     2.5012 </r>
       <r>    12.5526     0.3324     2.5411 </r>
       <r>    12.6729     0.2778     2.5746 </r>
       <r>    12.7932     0.2293     2.6021 </r>
       <r>    12.9134     0.2114     2.6276 </r>
       <r>    13.0337     0.2254     2.6547 </r>
       <r>    13.1540     0.2557     2.6854 </r>
       <r>    13.2743     0.2834     2.7195 </r>
       <r>    13.3945     0.2974     2.7553 </r>
       <r>    13.5148     0.2975     2.7910 </r>
       <r>    13.6351     0.2896     2.8259 </r>
       <r>    13.7553     0.2806     2.8596 </r>
       <r>    13.8756     0.2742     2.8926 </r>
       <r>    13.9959     0.2717     2.9253 </r>
       <r>    14.1162     0.2735     2.9582 </r>
       <r>    14.2364     0.2783     2.9917 </r>
       <r>    14.3567     0.2819     3.0256 </r>
       <r>    14.4770     0.2781     3.0590 </r>
       <r>    14.5972     0.2649     3.0909 </r>
       <r>    14.7175     0.2501     3.1209 </r>
       <r>    14.8378     0.2487     3.1509 </r>
       <r>    14.9581     0.2706     3.1834 </r>
       <r>    15.0783     0.3092     3.2206 </r>
       <r>    15.1986     0.3421     3.2617 </r>
       <r>    15.3189     0.3448     3.3032 </r>
       <r>    15.4391     0.3085     3.3403 </r>
       <r>    15.5594     0.2482     3.3702 </r>
       <r>    15.6797     0.1953     3.3936 </r>
       <r>    15.8000     0.1774     3.4150 </r>
       <r>    15.9202     0.2011     3.4392 </r>
       <r>    16.0405     0.2484     3.4690 </r>
       <r>    16.1608     0.2916     3.5041 </r>
       <r>    16.2810     0.3124     3.5417 </r>
       <r>    16.4013     0.3105     3.5790 </r>
       <r>    16.5216     0.2968     3.6147 </r>
       <r>    16.6419     0.2801     3.6484 </r>
       <r>    16.7621     0.2620     3.6799 </r>
       <r>    16.8824     0.2420     3.7090 </r>
       <r>    17.0027     0.2250     3.7361 </r>
       <r>    17.1229     0.2194     3.7625 </r>
       <r>    17.2432     0.2286     3.7899 </r>
       <r>    17.3635     0.2457     3.8195 </r>
       <r>    17.4838     0.2607     3.8509 </r>
       <r>    17.6040     0.2719     3.8836 </r>
       <r>    17.7243     0.2865     3.9180 </r>
       <r>    17.8446     0.3080     3.9550 </r>
       <r>    17.9648     0.3259     3.9942 </r>
       <r>    18.0851     0.3216     4.0329 </r>
       <r>    18.2054     0.2861     4.0673 </r>
       <r>    18.3257     0.2318     4.0952 </r>
       <r>    18.4459     0.1861     4.1176 </r>
       <r>    18.5662     0.1740     4.1385 </r>
       <r>    18.6865     0.2020     4.1628 </r>
       <r>    18.8067     0.2548     4.1935 </r>
       <r>    18.9270     0.3051     4.2302 </r>
       <r>    19.0473     0.3316     4.2700 </r>
       <r>    19.1676     0.3299     4.3097 </r>
       <r>    19.2878     0.3105     4.3471 </r>
       <r>    19.4081     0.2869     4.3816 </r>
       <r>    19.5284     0.2668     4.4137 </r>
       <r>    19.6486     0.2528     4.4441 </r>
       <r>    19.7689     0.2468     4.4738 </r>
       <r>    19.8892     0.2493     4.5037 </r>
       <r>    20.0095     0.2572     4.5347 </r>
       <r>    20.1297     0.2652     4.5666 </r>
       <r>    20.2500     0.2716     4.5992 </r>
       <r>    20.3703     0.2785     4.6327 </r>
       <r>    20.4905     0.2843     4.6669 </r>
       <r>    20.6108     0.2794     4.7005 </r>
       <r>    20.7311     0.2542     4.7311 </r>
       <r>    20.8514     0.2121     4.7566 </r>
       <r>    20.9716     0.1712     4.7772 </r>
       <r>    21.0919     0.1511     4.7954 </r>
       <r>    21.2122     0.1566     4.8142 </r>
       <r>    21.3324     0.1765     4.8354 </r>
       <r>    21.4527     0.1936     4.8587 </r>
       <r>    21.5730     0.1963     4.8823 </r>
       <r>    21.6933     0.1814     4.9042 </r>
       <r>    21.8135     0.1528     4.9225 </r>
       <r>    21.9338     0.1178     4.9367 </r>
       <r>    22.0541     0.0858     4.9470 </r>
       <r>    22.1743     0.0638     4.9547 </r>
       <r>    22.2946     0.0537     4.9611 </r>
       <r>    22.4149     0.0524     4.9675 </r>
       <r>    22.5352     0.0547     4.9740 </r>
       <r>    22.6554     0.0558     4.9807 </r>
       <r>    22.7757     0.0532     4.9871 </r>
       <r>    22.8960     0.0459     4.9927 </r>
       <r>    23.0162     0.0345     4.9968 </r>
       <r>    23.1365     0.0214     4.9994 </r>
       <r>    23.2568     0.0096     5.0005 </r>
       <r>    23.3771     0.0015     5.0007 </r>
       <r>    23.4973    -0.0021     5.0005 </r>
       <r>    23.6176    -0.0025     5.0001 </r>
       <r>    23.7379    -0.0015     5.0000 </r>
       <r>    23.8581    -0.0005     4.9999 </r>
       <r>    23.9784     0.0001     4.9999 </r>
       <r>    24.0987     0.0003     5.0000 </r>
       <r>    24.2190     0.0002     5.0000 </r>
       <r>    24.3392     0.0001     5.0000 </r>
       <r>    24.4595     0.0000     5.0000 </r>
       <r>    24.5798     0.0000     5.0000 </r>
       <r>    24.7000     0.0000     5.0000 </r>
       <r>    24.8203     0.0000     5.0000 </r>
       <r>    24.9406     0.0000     5.0000 </r>
       <r>    25.0609     0.0000     5.0000 </r>
       <r>    25.1811     0.0000     5.0000 </r>
       <r>    25.3014     0.0000     5.0000 </r>
       <r>    25.4217     0.0000     5.0000 </r>
       <r>    25.5419     0.0000     5.0000 </r>
       <r>    25.6622     0.0000     5.0000 </r>
       <r>    25.7825     0.0000     5.0000 </r>
       <r>    25.9028     0.0000     5.0000 </r>
       <r>    26.0230     0.0000     5.0000 </r>
       <r>    26.1433     0.0000     5.0000 </r>
       <r>    26.2636     0.0000     5.0000 </r>
       <r>    26.3838     0.0000     5.0000 </r>
       <r>    26.5041     0.0000     5.0000 </r>
       <r>    26.6244     0.0000     5.0000 </r>
       <r>    26.7447     0.0000     5.0000 </r>
       <r>    26.8649     0.0000     5.0000 </r>
       <r>    26.9852     0.0000     5.0000 </r>
       <r>    27.1055     0.0000     5.0000 </r>
       <r>    27.2257     0.0000     5.0000 </r>
       <r>    27.3460     0.0000     5.0000 </r>
       <r>    27.4663     0.0000     5.0000 </r>
       <r>    27.5866     0.0000     5.0000 </r>
       <r>    27.7068     0.0000     5.0000 </r>
       <r>    27.8271     0.0000     5.0000 </r>
       <r>    27.9474     0.0000     5.0000 </r>
      </set>
     </set>
    </array>
   </total>
   <partial>
    <array>
     <dimension dim="1">gridpoints</dimension>
     <dimension dim="2">spin</dimension>
     <dimension dim="3">ion</dimension>
     <field>energy</field>
     <field>    s</field>
     <field>   py</field>
     <field>   pz</field>
     <field>   px</field>
     <field>  dxy</field>
     <field>  dyz</field>
     <field>  dz2</field>
     <field>  dxz</field>
     <field>x2-y2</field>
     <set>
      <set comment="ion 1">
       <set comment="spin 1">
        <r>    -8.1341     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -8.0138     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.8935     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.7732     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.6530     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.5327     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.4124     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.2922     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.1719     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.0516     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.9313     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.8111     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.6908     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.5705     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.4503     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.3300     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.2097     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.0894     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.9692     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.8489     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.7286     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.6084     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.4881     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.3678     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.2475     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.1273     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.0070     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.8867     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.7665     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.6462     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.5259     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.4056     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.2854     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.1651     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.0448     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.9246     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.8043    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.6840    -0.0001     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.5637    -0.0003    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.4435    -0.0003    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.3232    -0.0000    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.2029     0.0008    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.0827     0.0023     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.9624     0.0044     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.8421     0.0067     0.0000     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.7218     0.0088     0.0001     0.0001     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.6016     0.0104     0.0001     0.0002     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.4813     0.0113     0.0001     0.0003     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.3610     0.0121     0.0002     0.0003     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.2408     0.0133     0.0003     0.0005     0.0002     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.1205     0.0149     0.0003     0.0006     0.0002     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.0002     0.0162     0.0003     0.0008     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.8799     0.0165     0.0003     0.0010     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.7597     0.0156     0.0003     0.0010     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.6394     0.0148     0.0003     0.0011     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.5191     0.0153     0.0005     0.0011     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.3989     0.0175     0.0007     0.0014     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.2786     0.0204     0.0009     0.0017     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.1583     0.0223     0.0010     0.0021     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.0380     0.0219     0.0008     0.0024     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.9178     0.0197     0.0006     0.0025     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.7975     0.0171     0.0005     0.0022     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.6772     0.0157     0.0006     0.0020     0.0007     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.5570     0.0168     0.0008     0.0018     0.0011     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.4367     0.0200     0.0013     0.0020     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.3164     0.0240     0.0017     0.0026     0.0018     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.1961     0.0268     0.0019     0.0034     0.0017     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.0759     0.0269     0.0019     0.0039     0.0014     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.0444     0.0242     0.0016     0.0041     0.0009     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.1647     0.0203     0.0012     0.0041     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.2849     0.0175     0.0009     0.0042     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.4052     0.0172     0.0008     0.0042     0.0006     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.5255     0.0191     0.0011     0.0041     0.0014     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.6458     0.0217     0.0017     0.0036     0.0024     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.7660     0.0237     0.0024     0.0031     0.0033     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.8863     0.0249     0.0028     0.0031     0.0037     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.0066     0.0260     0.0030     0.0041     0.0033     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.1268     0.0271     0.0030     0.0060     0.0025     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.2471     0.0277     0.0028     0.0077     0.0016     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.3674     0.0265     0.0025     0.0084     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.4877     0.0231     0.0021     0.0080     0.0006     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.6079     0.0187     0.0016     0.0068     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.7282     0.0155     0.0013     0.0056     0.0008     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.8485     0.0154     0.0016     0.0051     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.9687     0.0185     0.0023     0.0052     0.0028     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.0890     0.0235     0.0033     0.0059     0.0044     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.2093     0.0282     0.0042     0.0069     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.3296     0.0304     0.0047     0.0078     0.0060     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.4498     0.0294     0.0046     0.0083     0.0053     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.5701     0.0262     0.0041     0.0087     0.0040     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.6904     0.0229     0.0036     0.0092     0.0026     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.8106     0.0211     0.0034     0.0102     0.0016     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.9309     0.0209     0.0033     0.0115     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.0512     0.0206     0.0032     0.0123     0.0009     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.1715     0.0190     0.0030     0.0119     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.2917     0.0165     0.0027     0.0105     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.4120     0.0149     0.0029     0.0090     0.0025     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.5323     0.0161     0.0037     0.0085     0.0042     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.6525     0.0199     0.0049     0.0091     0.0064     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.7728     0.0243     0.0061     0.0101     0.0084     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.8931     0.0270     0.0069     0.0108     0.0095     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.0134     0.0270     0.0070     0.0110     0.0095     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.1336     0.0247     0.0064     0.0112     0.0082     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.2539     0.0216     0.0055     0.0118     0.0060     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.3742     0.0191     0.0044     0.0130     0.0038     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.4944     0.0182     0.0034     0.0145     0.0021     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.6147     0.0189     0.0029     0.0159     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.7350     0.0202     0.0029     0.0169     0.0017     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.8553     0.0210     0.0032     0.0172     0.0025     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.9755     0.0200     0.0036     0.0172     0.0034     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.0958     0.0178     0.0041     0.0168     0.0042     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.2161     0.0155     0.0048     0.0162     0.0053     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.3363     0.0147     0.0057     0.0153     0.0070     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.4566     0.0160     0.0066     0.0142     0.0092     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.5769     0.0183     0.0072     0.0131     0.0111     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.6972     0.0203     0.0072     0.0126     0.0121     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.8174     0.0213     0.0067     0.0130     0.0120     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.9377     0.0215     0.0061     0.0143     0.0114     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.0580     0.0215     0.0056     0.0158     0.0108     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.1782     0.0211     0.0052     0.0168     0.0102     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.2985     0.0197     0.0047     0.0167     0.0089     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.4188     0.0170     0.0038     0.0155     0.0068     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.5391     0.0140     0.0030     0.0142     0.0041     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.6593     0.0123     0.0025     0.0139     0.0018     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.7796     0.0125     0.0028     0.0147     0.0009     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.8999     0.0140     0.0038     0.0154     0.0016     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.0201     0.0155     0.0054     0.0144     0.0042     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.1404     0.0161     0.0071     0.0116     0.0082     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.2607     0.0164     0.0084     0.0081     0.0129     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.3810     0.0172     0.0090     0.0056     0.0174     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.5012     0.0189     0.0090     0.0054     0.0205     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.6215     0.0211     0.0085     0.0070     0.0215     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.7418     0.0229     0.0079     0.0096     0.0200     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.8620     0.0238     0.0073     0.0118     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.9823     0.0242     0.0066     0.0128     0.0126     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.1026     0.0244     0.0060     0.0121     0.0097     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.2229     0.0244     0.0055     0.0103     0.0087     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.3431     0.0235     0.0054     0.0087     0.0090     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.4634     0.0209     0.0057     0.0082     0.0091     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.5837     0.0170     0.0059     0.0089     0.0083     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.7039     0.0132     0.0058     0.0102     0.0067     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.8242     0.0112     0.0057     0.0113     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.9445     0.0118     0.0059     0.0120     0.0060     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.0648     0.0153     0.0068     0.0125     0.0081     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.1850     0.0210     0.0084     0.0131     0.0113     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.3053     0.0273     0.0099     0.0135     0.0145     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.4256     0.0319     0.0107     0.0133     0.0169     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.5458     0.0329     0.0106     0.0126     0.0178     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.6661     0.0301     0.0100     0.0118     0.0169     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.7864     0.0246     0.0092     0.0113     0.0148     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.9067     0.0190     0.0085     0.0108     0.0126     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.0269     0.0151     0.0079     0.0103     0.0112     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.1472     0.0140     0.0076     0.0101     0.0109     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.2675     0.0154     0.0078     0.0104     0.0108     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.3877     0.0183     0.0088     0.0112     0.0101     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.5080     0.0216     0.0100     0.0119     0.0086     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.6283     0.0243     0.0107     0.0117     0.0070     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.7486     0.0259     0.0106     0.0107     0.0068     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.8688     0.0268     0.0098     0.0099     0.0089     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.9891     0.0281     0.0094     0.0103     0.0132     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.1094     0.0300     0.0101     0.0119     0.0182     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.2296     0.0315     0.0119     0.0139     0.0221     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.3499     0.0307     0.0137     0.0147     0.0232     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.4702     0.0264     0.0143     0.0137     0.0212     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.5905     0.0194     0.0130     0.0111     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.7107     0.0124     0.0102     0.0082     0.0111     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.8310     0.0085     0.0074     0.0067     0.0071     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.9513     0.0094     0.0063     0.0077     0.0061     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.0715     0.0151     0.0084     0.0113     0.0086     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.1918     0.0231     0.0133     0.0161     0.0131     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.3121     0.0298     0.0190     0.0198     0.0171     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.4324     0.0323     0.0229     0.0204     0.0185     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.5526     0.0301     0.0235     0.0182     0.0169     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.6729     0.0257     0.0218     0.0148     0.0136     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.7932     0.0222     0.0196     0.0123     0.0106     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.9134     0.0213     0.0183     0.0120     0.0095     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.0337     0.0222     0.0181     0.0133     0.0107     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.1540     0.0232     0.0182     0.0151     0.0135     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.2743     0.0224     0.0181     0.0161     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.3945     0.0200     0.0180     0.0159     0.0193     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.5148     0.0170     0.0185     0.0146     0.0211     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.6351     0.0146     0.0200     0.0132     0.0220     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.7553     0.0135     0.0222     0.0123     0.0225     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.8756     0.0137     0.0238     0.0122     0.0230     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.9959     0.0145     0.0240     0.0127     0.0240     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.1162     0.0156     0.0229     0.0137     0.0258     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.2364     0.0164     0.0215     0.0148     0.0276     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.3567     0.0166     0.0207     0.0155     0.0284     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.4770     0.0162     0.0204     0.0153     0.0274     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.5972     0.0152     0.0199     0.0139     0.0250     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.7175     0.0137     0.0194     0.0116     0.0229     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.8378     0.0122     0.0196     0.0092     0.0233     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.9581     0.0112     0.0217     0.0076     0.0271     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.0783     0.0113     0.0250     0.0070     0.0331     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.1986     0.0121     0.0274     0.0069     0.0379     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.3189     0.0130     0.0267     0.0066     0.0385     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.4391     0.0131     0.0221     0.0062     0.0335     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.5594     0.0124     0.0154     0.0065     0.0247     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.6797     0.0115     0.0098     0.0083     0.0158     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.8000     0.0112     0.0082     0.0119     0.0104     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.9202     0.0118     0.0110     0.0160     0.0103     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.0405     0.0125     0.0160     0.0189     0.0143     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.1608     0.0127     0.0203     0.0193     0.0198     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.2810     0.0123     0.0217     0.0174     0.0245     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.4013     0.0118     0.0201     0.0146     0.0271     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.5216     0.0120     0.0168     0.0125     0.0272     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.6419     0.0130     0.0134     0.0124     0.0246     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.7621     0.0140     0.0110     0.0139     0.0201     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.8824     0.0141     0.0097     0.0155     0.0153     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.0027     0.0127     0.0098     0.0156     0.0128     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.1229     0.0103     0.0112     0.0136     0.0138     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.2432     0.0078     0.0134     0.0108     0.0175     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.3635     0.0062     0.0150     0.0088     0.0212     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.4838     0.0065     0.0153     0.0085     0.0227     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.6040     0.0091     0.0142     0.0099     0.0217     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.7243     0.0134     0.0130     0.0122     0.0200     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.8446     0.0177     0.0129     0.0146     0.0194     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.9648     0.0198     0.0140     0.0159     0.0202     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.0851     0.0186     0.0153     0.0153     0.0208     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.2054     0.0147     0.0154     0.0127     0.0196     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.3257     0.0102     0.0139     0.0090     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.4459     0.0072     0.0117     0.0060     0.0131     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.5662     0.0068     0.0108     0.0053     0.0112     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.6865     0.0088     0.0119     0.0073     0.0118     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.8067     0.0119     0.0144     0.0112     0.0139     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.9270     0.0147     0.0165     0.0155     0.0160     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.0473     0.0159     0.0170     0.0187     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.1676     0.0153     0.0161     0.0202     0.0154     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.2878     0.0134     0.0149     0.0203     0.0136     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.4081     0.0112     0.0144     0.0195     0.0122     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.5284     0.0093     0.0148     0.0181     0.0118     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.6486     0.0081     0.0157     0.0163     0.0123     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.7689     0.0075     0.0165     0.0143     0.0132     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.8892     0.0076     0.0171     0.0125     0.0137     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.0095     0.0080     0.0175     0.0112     0.0136     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.1297     0.0085     0.0179     0.0104     0.0132     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.2500     0.0088     0.0190     0.0098     0.0131     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.3703     0.0088     0.0211     0.0091     0.0139     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.4905     0.0086     0.0239     0.0079     0.0152     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.6108     0.0079     0.0257     0.0062     0.0159     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.7311     0.0068     0.0247     0.0046     0.0150     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.8514     0.0052     0.0207     0.0038     0.0126     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.9716     0.0037     0.0154     0.0045     0.0099     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.0919     0.0029     0.0111     0.0066     0.0083     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.2122     0.0030     0.0096     0.0091     0.0081     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.3324     0.0037     0.0110     0.0109     0.0085     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.4527     0.0045     0.0141     0.0114     0.0087     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.5730     0.0048     0.0172     0.0113     0.0080     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.6933     0.0044     0.0189     0.0110     0.0066     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.8135     0.0037     0.0183     0.0107     0.0051     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.9338     0.0029     0.0155     0.0099     0.0037     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.0541     0.0022     0.0118     0.0081     0.0031     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.1743     0.0017     0.0089     0.0055     0.0034     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.2946     0.0014     0.0074     0.0028     0.0044     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.4149     0.0013     0.0070     0.0009     0.0057     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.5352     0.0013     0.0068     0.0001     0.0068     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.6554     0.0011     0.0061     0.0001     0.0072     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.7757     0.0008     0.0047     0.0004     0.0068     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.8960     0.0005     0.0029     0.0005     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.0162     0.0002     0.0012     0.0004     0.0040     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.1365     0.0000     0.0002     0.0002     0.0024     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.2568    -0.0000    -0.0003     0.0001     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.3771    -0.0001    -0.0003    -0.0000     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.4973    -0.0000    -0.0002    -0.0000    -0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.6176    -0.0000    -0.0001    -0.0000    -0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.7379     0.0000     0.0000    -0.0000    -0.0002     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.8581     0.0000     0.0000     0.0000    -0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.9784     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.0987     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.2190     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.3392     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.4595     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.5798     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.7000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.8203     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.9406     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.0609     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.1811     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.3014     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.4217     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.5419     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.6622     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.7825     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.9028     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.0230     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.1433     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.2636     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.3838     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.5041     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.6244     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.7447     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.8649     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.9852     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.1055     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.2257     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.3460     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.4663     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.5866     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.7068     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.8271     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.9474     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
       </set>
       <set comment="spin 2">
        <r>    -8.1341     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -8.0138     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.8935     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.7732     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.6530     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.5327     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.4124     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.2922     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.1719     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -7.0516     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.9313     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.8111     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.6908     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.5705     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.4503     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.3300     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.2097     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -6.0894     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.9692     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.8489     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.7286     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.6084     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.4881     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.3678     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.2475     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.1273     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -5.0070     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.8867     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.7665     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.6462     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.5259     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.4056     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.2854     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.1651     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -4.0448     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.9246     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.8043    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.6840    -0.0001     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.5637    -0.0003    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.4435    -0.0003    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.3232    -0.0000    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.2029     0.0008    -0.0000    -0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -3.0827     0.0023     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.9624     0.0044     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.8421     0.0067     0.0000     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.7218     0.0088     0.0001     0.0001     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.6016     0.0104     0.0001     0.0002     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.4813     0.0113     0.0001     0.0003     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.3610     0.0121     0.0002     0.0003     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.2408     0.0133     0.0003     0.0005     0.0002     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.1205     0.0149     0.0003     0.0006     0.0002     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -2.0002     0.0162     0.0003     0.0008     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.8799     0.0165     0.0003     0.0010     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.7597     0.0156     0.0003     0.0010     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.6394     0.0148     0.0003     0.0011     0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.5191     0.0153     0.0005     0.0011     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.3989     0.0175     0.0007     0.0014     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.2786     0.0204     0.0009     0.0017     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.1583     0.0223     0.0010     0.0021     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -1.0380     0.0219     0.0008     0.0024     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.9178     0.0197     0.0006     0.0025     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.7975     0.0171     0.0005     0.0022     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.6772     0.0157     0.0006     0.0020     0.0007     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.5570     0.0168     0.0008     0.0018     0.0011     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.4367     0.0200     0.0013     0.0020     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.3164     0.0240     0.0017     0.0026     0.0018     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.1961     0.0268     0.0019     0.0034     0.0017     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    -0.0759     0.0269     0.0019     0.0039     0.0014     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.0444     0.0242     0.0016     0.0041     0.0009     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.1647     0.0203     0.0012     0.0041     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.2849     0.0175     0.0009     0.0042     0.0004     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.4052     0.0172     0.0008     0.0042     0.0006     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.5255     0.0191     0.0011     0.0041     0.0014     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.6458     0.0217     0.0017     0.0036     0.0024     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.7660     0.0237     0.0024     0.0031     0.0033     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     0.8863     0.0249     0.0028     0.0031     0.0037     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.0066     0.0260     0.0030     0.0041     0.0033     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.1268     0.0271     0.0030     0.0060     0.0025     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.2471     0.0277     0.0028     0.0077     0.0016     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.3674     0.0265     0.0025     0.0084     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.4877     0.0231     0.0021     0.0080     0.0006     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.6079     0.0187     0.0016     0.0068     0.0005     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.7282     0.0155     0.0013     0.0056     0.0008     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.8485     0.0154     0.0016     0.0051     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     1.9687     0.0185     0.0023     0.0052     0.0028     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.0890     0.0235     0.0033     0.0059     0.0044     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.2093     0.0282     0.0042     0.0069     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.3296     0.0304     0.0047     0.0078     0.0060     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.4498     0.0294     0.0046     0.0083     0.0053     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.5701     0.0262     0.0041     0.0087     0.0040     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.6904     0.0229     0.0036     0.0092     0.0026     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.8106     0.0211     0.0034     0.0102     0.0016     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     2.9309     0.0209     0.0033     0.0115     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.0512     0.0206     0.0032     0.0123     0.0009     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.1715     0.0190     0.0030     0.0119     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.2917     0.0165     0.0027     0.0105     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.4120     0.0149     0.0029     0.0090     0.0025     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.5323     0.0161     0.0037     0.0085     0.0042     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.6525     0.0199     0.0049     0.0091     0.0064     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.7728     0.0243     0.0061     0.0101     0.0084     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     3.8931     0.0270     0.0069     0.0108     0.0095     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.0134     0.0270     0.0070     0.0110     0.0095     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.1336     0.0247     0.0064     0.0112     0.0082     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.2539     0.0216     0.0055     0.0118     0.0060     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.3742     0.0191     0.0044     0.0130     0.0038     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.4944     0.0182     0.0034     0.0145     0.0021     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.6147     0.0189     0.0029     0.0159     0.0015     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.7350     0.0202     0.0029     0.0169     0.0017     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.8553     0.0210     0.0032     0.0172     0.0025     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     4.9755     0.0200     0.0036     0.0172     0.0034     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.0958     0.0178     0.0041     0.0168     0.0042     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.2161     0.0155     0.0048     0.0162     0.0053     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.3363     0.0147     0.0057     0.0153     0.0070     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.4566     0.0160     0.0066     0.0142     0.0092     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.5769     0.0183     0.0072     0.0131     0.0111     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.6972     0.0203     0.0072     0.0126     0.0121     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.8174     0.0213     0.0067     0.0130     0.0120     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     5.9377     0.0215     0.0061     0.0143     0.0114     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.0580     0.0215     0.0056     0.0158     0.0108     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.1782     0.0211     0.0052     0.0168     0.0102     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.2985     0.0197     0.0047     0.0167     0.0089     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.4188     0.0170     0.0038     0.0155     0.0068     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.5391     0.0140     0.0030     0.0142     0.0041     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.6593     0.0123     0.0025     0.0139     0.0018     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.7796     0.0125     0.0028     0.0147     0.0009     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     6.8999     0.0140     0.0038     0.0154     0.0016     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.0201     0.0155     0.0054     0.0144     0.0042     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.1404     0.0161     0.0071     0.0116     0.0082     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.2607     0.0164     0.0084     0.0081     0.0129     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.3810     0.0172     0.0090     0.0056     0.0174     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.5012     0.0189     0.0090     0.0054     0.0205     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.6215     0.0211     0.0085     0.0070     0.0215     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.7418     0.0229     0.0079     0.0096     0.0200     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.8620     0.0238     0.0073     0.0118     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     7.9823     0.0242     0.0066     0.0128     0.0126     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.1026     0.0244     0.0060     0.0121     0.0097     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.2229     0.0244     0.0055     0.0103     0.0087     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.3431     0.0235     0.0054     0.0087     0.0090     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.4634     0.0209     0.0057     0.0082     0.0091     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.5837     0.0170     0.0059     0.0089     0.0083     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.7039     0.0132     0.0058     0.0102     0.0067     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.8242     0.0112     0.0057     0.0113     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     8.9445     0.0118     0.0059     0.0120     0.0060     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.0648     0.0153     0.0068     0.0125     0.0081     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.1850     0.0210     0.0084     0.0131     0.0113     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.3053     0.0273     0.0099     0.0135     0.0145     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.4256     0.0319     0.0107     0.0133     0.0169     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.5458     0.0329     0.0106     0.0126     0.0178     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.6661     0.0301     0.0100     0.0118     0.0169     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.7864     0.0246     0.0092     0.0113     0.0148     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>     9.9067     0.0190     0.0085     0.0108     0.0126     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.0269     0.0151     0.0079     0.0103     0.0112     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.1472     0.0140     0.0076     0.0101     0.0109     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.2675     0.0154     0.0078     0.0104     0.0108     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.3877     0.0183     0.0088     0.0112     0.0101     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.5080     0.0216     0.0100     0.0119     0.0086     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.6283     0.0243     0.0107     0.0117     0.0070     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.7486     0.0259     0.0106     0.0107     0.0068     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.8688     0.0268     0.0098     0.0099     0.0089     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    10.9891     0.0281     0.0094     0.0103     0.0132     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.1094     0.0300     0.0101     0.0119     0.0182     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.2296     0.0315     0.0119     0.0139     0.0221     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.3499     0.0307     0.0137     0.0147     0.0232     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.4702     0.0264     0.0143     0.0137     0.0212     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.5905     0.0194     0.0130     0.0111     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.7107     0.0124     0.0102     0.0082     0.0111     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.8310     0.0085     0.0074     0.0067     0.0071     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    11.9513     0.0094     0.0063     0.0077     0.0061     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.0715     0.0151     0.0084     0.0113     0.0086     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.1918     0.0231     0.0133     0.0161     0.0131     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.3121     0.0298     0.0190     0.0198     0.0171     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.4324     0.0323     0.0229     0.0204     0.0185     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.5526     0.0301     0.0235     0.0182     0.0169     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.6729     0.0257     0.0218     0.0148     0.0136     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.7932     0.0222     0.0196     0.0123     0.0106     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    12.9134     0.0213     0.0183     0.0120     0.0095     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.0337     0.0222     0.0181     0.0133     0.0107     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.1540     0.0232     0.0182     0.0151     0.0135     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.2743     0.0224     0.0181     0.0161     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.3945     0.0200     0.0180     0.0159     0.0193     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.5148     0.0170     0.0185     0.0146     0.0211     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.6351     0.0146     0.0200     0.0132     0.0220     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.7553     0.0135     0.0222     0.0123     0.0225     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.8756     0.0137     0.0238     0.0122     0.0230     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    13.9959     0.0145     0.0240     0.0127     0.0240     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.1162     0.0156     0.0229     0.0137     0.0258     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.2364     0.0164     0.0215     0.0148     0.0276     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.3567     0.0166     0.0207     0.0155     0.0284     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.4770     0.0162     0.0204     0.0153     0.0274     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.5972     0.0152     0.0199     0.0139     0.0250     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.7175     0.0137     0.0194     0.0116     0.0229     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.8378     0.0122     0.0196     0.0092     0.0233     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    14.9581     0.0112     0.0217     0.0076     0.0271     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.0783     0.0113     0.0250     0.0070     0.0331     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.1986     0.0121     0.0274     0.0069     0.0379     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.3189     0.0130     0.0267     0.0066     0.0385     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.4391     0.0131     0.0221     0.0062     0.0335     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.5594     0.0124     0.0154     0.0065     0.0247     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.6797     0.0115     0.0098     0.0083     0.0158     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.8000     0.0112     0.0082     0.0119     0.0104     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    15.9202     0.0118     0.0110     0.0160     0.0103     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.0405     0.0125     0.0160     0.0189     0.0143     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.1608     0.0127     0.0203     0.0193     0.0198     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.2810     0.0123     0.0217     0.0174     0.0245     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.4013     0.0118     0.0201     0.0146     0.0271     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.5216     0.0120     0.0168     0.0125     0.0272     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.6419     0.0130     0.0134     0.0124     0.0246     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.7621     0.0140     0.0110     0.0139     0.0201     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    16.8824     0.0141     0.0097     0.0155     0.0153     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.0027     0.0127     0.0098     0.0156     0.0128     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.1229     0.0103     0.0112     0.0136     0.0138     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.2432     0.0078     0.0134     0.0108     0.0175     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.3635     0.0062     0.0150     0.0088     0.0212     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.4838     0.0065     0.0153     0.0085     0.0227     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.6040     0.0091     0.0142     0.0099     0.0217     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.7243     0.0134     0.0130     0.0122     0.0200     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.8446     0.0177     0.0129     0.0146     0.0194     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    17.9648     0.0198     0.0140     0.0159     0.0202     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.0851     0.0186     0.0153     0.0153     0.0208     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.2054     0.0147     0.0154     0.0127     0.0196     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.3257     0.0102     0.0139     0.0090     0.0166     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.4459     0.0072     0.0118     0.0060     0.0131     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.5662     0.0068     0.0108     0.0053     0.0112     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.6865     0.0088     0.0119     0.0073     0.0117     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.8067     0.0119     0.0144     0.0112     0.0139     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    18.9270     0.0147     0.0165     0.0155     0.0160     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.0473     0.0159     0.0170     0.0187     0.0165     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.1676     0.0153     0.0160     0.0203     0.0154     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.2878     0.0134     0.0148     0.0204     0.0136     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.4081     0.0112     0.0143     0.0196     0.0121     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.5284     0.0093     0.0147     0.0184     0.0116     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.6486     0.0081     0.0154     0.0168     0.0120     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.7689     0.0075     0.0162     0.0150     0.0128     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    19.8892     0.0076     0.0167     0.0132     0.0133     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.0095     0.0080     0.0172     0.0118     0.0133     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.1297     0.0085     0.0178     0.0107     0.0130     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.2500     0.0088     0.0191     0.0099     0.0131     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.3703     0.0087     0.0217     0.0089     0.0140     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.4905     0.0084     0.0250     0.0074     0.0153     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.6108     0.0078     0.0270     0.0056     0.0159     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.7311     0.0066     0.0260     0.0040     0.0149     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.8514     0.0051     0.0217     0.0033     0.0124     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    20.9716     0.0037     0.0159     0.0043     0.0098     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.0919     0.0029     0.0113     0.0066     0.0083     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.2122     0.0031     0.0096     0.0092     0.0081     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.3324     0.0038     0.0110     0.0110     0.0085     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.4527     0.0044     0.0141     0.0118     0.0086     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.5730     0.0044     0.0173     0.0119     0.0079     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.6933     0.0039     0.0189     0.0120     0.0065     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.8135     0.0031     0.0182     0.0119     0.0050     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    21.9338     0.0024     0.0152     0.0110     0.0037     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.0541     0.0018     0.0115     0.0089     0.0031     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.1743     0.0015     0.0085     0.0059     0.0034     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.2946     0.0014     0.0070     0.0029     0.0043     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.4149     0.0014     0.0067     0.0010     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.5352     0.0013     0.0065     0.0003     0.0066     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.6554     0.0012     0.0059     0.0005     0.0070     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.7757     0.0009     0.0044     0.0008     0.0067     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    22.8960     0.0005     0.0027     0.0010     0.0056     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.0162     0.0002     0.0011     0.0007     0.0040     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.1365     0.0000     0.0001     0.0004     0.0024     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.2568    -0.0001    -0.0003     0.0001     0.0010     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.3771    -0.0001    -0.0003    -0.0001     0.0001     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.4973    -0.0000    -0.0002    -0.0001    -0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.6176    -0.0000    -0.0000    -0.0000    -0.0003     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.7379     0.0000     0.0000    -0.0000    -0.0002     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.8581     0.0000     0.0000     0.0000    -0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    23.9784     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.0987     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.2190     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.3392     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.4595     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.5798     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.7000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.8203     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    24.9406     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.0609     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.1811     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.3014     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.4217     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.5419     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.6622     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.7825     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    25.9028     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.0230     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.1433     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.2636     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.3838     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.5041     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.6244     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.7447     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.8649     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    26.9852     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.1055     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.2257     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.3460     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.4663     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.5866     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.7068     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.8271     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
        <r>    27.9474     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000     0.0000 </r>
       </set>
      </set>
     </set>
    </array>
   </partial>
  </dos>
  <projected>
   <eigenvalues>
    <array>
     <dimension dim="1">band</dimension>
     <dimension dim="2">kpoint</dimension>
     <dimension dim="3">spin</dimension>
     <field>eigene</field>
     <field>occ</field>
     <set>
      <set comment="spin 1">
       <set comment="kpoint 1">
        <r>   -3.1341    1.0000 </r>
        <r>   20.6208    0.0000 </r>
        <r>   20.6208    0.0000 </r>
        <r>   20.6208    0.0000 </r>
        <r>   22.1131    0.0000 </r>
       </set>
       <set comment="kpoint 2">
        <r>   -2.9732    1.0000 </r>
        <r>   18.2625    0.0000 </r>
        <r>   20.1120    0.0000 </r>
        <r>   20.1120    0.0000 </r>
        <r>   22.7761    0.0000 </r>
       </set>
       <set comment="kpoint 3">
        <r>   -2.4918    1.0000 </r>
        <r>   15.1545    0.0000 </r>
        <r>   19.4447    0.0000 </r>
        <r>   19.4447    0.0000 </r>
        <r>   22.6316    0.0000 </r>
       </set>
       <set comment="kpoint 4">
        <r>   -1.6930    1.0000 </r>
        <r>   12.1586    0.0000 </r>
        <r>   18.9854    0.0000 </r>
        <r>   18.9854    0.0000 </r>
        <r>   22.3991    0.0000 </r>
       </set>
       <set comment="kpoint 5">
        <r>   -0.5822    1.0000 </r>
        <r>    9.3722    0.0005 </r>
        <r>   18.7803    0.0000 </r>
        <r>   18.7805    0.0000 </r>
        <r>   21.7676    0.0000 </r>
       </set>
       <set comment="kpoint 6">
        <r>    0.8335    1.0000 </r>
        <r>    6.8308    0.9959 </r>
        <r>   18.8406    0.0000 </r>
        <r>   18.8412    0.0000 </r>
        <r>   20.6402    0.0000 </r>
       </set>
       <set comment="kpoint 7">
        <r>    2.5434    1.0000 </r>
        <r>    4.5561    1.0000 </r>
        <r>   19.1696    0.0000 </r>
        <r>   19.1756    0.0000 </r>
        <r>   19.7719    0.0000 </r>
       </set>
       <set comment="kpoint 8">
        <r>   -2.9196    1.0000 </r>
        <r>   19.0625    0.0000 </r>
        <r>   19.0625    0.0000 </r>
        <r>   19.3530    0.0000 </r>
        <r>   21.6182    0.0000 </r>
       </set>
       <set comment="kpoint 9">
        <r>   -2.5452    1.0000 </r>
        <r>   16.2579    0.0000 </r>
        <r>   18.2305    0.0000 </r>
        <r>   19.4493    0.0000 </r>
        <r>   21.9393    0.0000 </r>
       </set>
       <set comment="kpoint 10">
        <r>   -1.8524    1.0000 </r>
        <r>   13.2453    0.0000 </r>
        <r>   17.6567    0.0000 </r>
        <r>   19.1235    0.0000 </r>
        <r>   21.3429    0.0000 </r>
       </set>
       <set comment="kpoint 11">
        <r>   -0.8459    1.0000 </r>
        <r>   10.3980    0.0000 </r>
        <r>   17.3504    0.0000 </r>
        <r>   18.8807    0.0000 </r>
        <r>   19.8568    0.0000 </r>
       </set>
       <set comment="kpoint 12">
        <r>    0.4677    1.0000 </r>
        <r>    7.7792    0.8110 </r>
        <r>   17.3148    0.0000 </r>
        <r>   18.5785    0.0000 </r>
        <r>   18.9007    0.0000 </r>
       </set>
       <set comment="kpoint 13">
        <r>    2.0796    1.0000 </r>
        <r>    5.4172    1.0000 </r>
        <r>   17.5514    0.0000 </r>
        <r>   17.6063    0.0000 </r>
        <r>   19.1519    0.0000 </r>
       </set>
       <set comment="kpoint 14">
        <r>    3.3035    1.0000 </r>
        <r>    4.0052    1.0000 </r>
        <r>   16.8983    0.0000 </r>
        <r>   18.0612    0.0000 </r>
        <r>   19.6932    0.0000 </r>
       </set>
       <set comment="kpoint 15">
        <r>    1.5100    1.0000 </r>
        <r>    6.1739    1.0000 </r>
        <r>   16.4703    0.0000 </r>
        <r>   18.8428    0.0000 </r>
        <r>   20.5110    0.0000 </r>
       </set>
       <set comment="kpoint 16">
        <r>   -0.0038    1.0000 </r>
        <r>    8.6224   -0.0218 </r>
        <r>   16.3249    0.0000 </r>
        <r>   19.8771    0.0000 </r>
        <r>   20.5167    0.0000 </r>
       </set>
       <set comment="kpoint 17">
        <r>   -1.2159    1.0000 </r>
        <r>   11.3163    0.0000 </r>
        <r>   16.4641    0.0000 </r>
        <r>   20.5551    0.0000 </r>
        <r>   21.2450    0.0000 </r>
       </set>
       <set comment="kpoint 18">
        <r>   -2.1185    1.0000 </r>
        <r>   14.2062    0.0000 </r>
        <r>   16.8998    0.0000 </r>
        <r>   20.8873    0.0000 </r>
        <r>   22.6024    0.0000 </r>
       </set>
       <set comment="kpoint 19">
        <r>   -2.7056    1.0000 </r>
        <r>   16.9661    0.0000 </r>
        <r>   17.8588    0.0000 </r>
        <r>   21.0378    0.0000 </r>
        <r>   22.3633    0.0000 </r>
       </set>
       <set comment="kpoint 20">
        <r>   -2.2783    1.0000 </r>
        <r>   17.2814    0.0000 </r>
        <r>   17.2814    0.0000 </r>
        <r>   17.4619    0.0000 </r>
        <r>   20.5141    0.0000 </r>
       </set>
       <set comment="kpoint 21">
        <r>   -1.6930    1.0000 </r>
        <r>   14.5332    0.0000 </r>
        <r>   16.6041    0.0000 </r>
        <r>   17.8754    0.0000 </r>
        <r>   19.2376    0.0000 </r>
       </set>
       <set comment="kpoint 22">
        <r>   -0.7932    1.0000 </r>
        <r>   11.6813    0.0000 </r>
        <r>   16.2000    0.0000 </r>
        <r>   17.5061    0.0000 </r>
        <r>   17.8750    0.0000 </r>
       </set>
       <set comment="kpoint 23">
        <r>    0.4150    1.0000 </r>
        <r>    9.0083    0.0046 </r>
        <r>   16.0694    0.0000 </r>
        <r>   16.1799    0.0000 </r>
        <r>   17.7810    0.0000 </r>
       </set>
       <set comment="kpoint 24">
        <r>    1.9234    1.0000 </r>
        <r>    6.5725    0.9991 </r>
        <r>   15.0818    0.0000 </r>
        <r>   16.2127    0.0000 </r>
        <r>   17.9735    0.0000 </r>
       </set>
       <set comment="kpoint 25">
        <r>    3.7068    1.0000 </r>
        <r>    4.4152    1.0000 </r>
        <r>   14.2658    0.0000 </r>
        <r>   16.6303    0.0000 </r>
        <r>   18.4467    0.0000 </r>
       </set>
       <set comment="kpoint 26">
        <r>    2.4919    1.0000 </r>
        <r>    5.8125    1.0000 </r>
        <r>   13.7414    0.0000 </r>
        <r>   17.3224    0.0000 </r>
        <r>   19.2038    0.0000 </r>
       </set>
       <set comment="kpoint 27">
        <r>    0.8854    1.0000 </r>
        <r>    8.1494    0.1093 </r>
        <r>   13.5200    0.0000 </r>
        <r>   18.2870    0.0000 </r>
        <r>   20.2350    0.0000 </r>
       </set>
       <set comment="kpoint 28">
        <r>   -0.4243    1.0000 </r>
        <r>   10.6937    0.0000 </r>
        <r>   13.6386    0.0000 </r>
        <r>   19.5157    0.0000 </r>
        <r>   21.5178    0.0000 </r>
       </set>
       <set comment="kpoint 29">
        <r>   -1.4278    1.0000 </r>
        <r>   13.0656    0.0000 </r>
        <r>   14.4433    0.0000 </r>
        <r>   20.9757    0.0000 </r>
        <r>   22.2850    0.0000 </r>
       </set>
       <set comment="kpoint 30">
        <r>   -1.2160    1.0000 </r>
        <r>   15.8284    0.0000 </r>
        <r>   15.8284    0.0000 </r>
        <r>   15.8293    0.0000 </r>
        <r>   17.2063    0.0000 </r>
       </set>
       <set comment="kpoint 31">
        <r>   -0.4246    1.0000 </r>
        <r>   13.1660    0.0000 </r>
        <r>   15.3287    0.0000 </r>
        <r>   15.4253    0.0000 </r>
        <r>   16.5538    0.0000 </r>
       </set>
       <set comment="kpoint 32">
        <r>    0.6755    1.0000 </r>
        <r>   10.4935    0.0000 </r>
        <r>   13.9343    0.0000 </r>
        <r>   15.1040    0.0000 </r>
        <r>   16.6700    0.0000 </r>
       </set>
       <set comment="kpoint 33">
        <r>    2.0762    1.0000 </r>
        <r>    8.0069    0.3565 </r>
        <r>   12.7214    0.0000 </r>
        <r>   15.1539    0.0000 </r>
        <r>   16.8304    0.0000 </r>
       </set>
       <set comment="kpoint 34">
        <r>    3.7624    1.0000 </r>
        <r>    5.7667    1.0000 </r>
        <r>   11.8000    0.0000 </r>
        <r>   15.4786    0.0000 </r>
        <r>   17.2330    0.0000 </r>
       </set>
       <set comment="kpoint 35">
        <r>    3.7703    1.0000 </r>
        <r>    5.7430    1.0000 </r>
        <r>   11.1866    0.0000 </r>
        <r>   16.0784    0.0000 </r>
        <r>   17.9154    0.0000 </r>
       </set>
       <set comment="kpoint 36">
        <r>    2.0791    1.0000 </r>
        <r>    7.9166    0.5452 </r>
        <r>   10.9288    0.0000 </r>
        <r>   16.9535    0.0000 </r>
        <r>   18.8771    0.0000 </r>
       </set>
       <set comment="kpoint 37">
        <r>    0.6766    1.0000 </r>
        <r>    9.9061    0.0000 </r>
        <r>   11.4181    0.0000 </r>
        <r>   18.1019    0.0000 </r>
        <r>   20.0968    0.0000 </r>
       </set>
       <set comment="kpoint 38">
        <r>    0.2572    1.0000 </r>
        <r>   13.5342    0.0000 </r>
        <r>   14.5772    0.0000 </r>
        <r>   14.7348    0.0000 </r>
        <r>   14.7348    0.0000 </r>
       </set>
       <set comment="kpoint 39">
        <r>    1.2467    1.0000 </r>
        <r>   11.8647    0.0000 </r>
        <r>   12.2444    0.0000 </r>
        <r>   14.4167    0.0000 </r>
        <r>   15.5526    0.0000 </r>
       </set>
       <set comment="kpoint 40">
        <r>    2.5349    1.0000 </r>
        <r>    9.6698    0.0000 </r>
        <r>   10.6384    0.0000 </r>
        <r>   14.3736    0.0000 </r>
        <r>   15.8675    0.0000 </r>
       </set>
       <set comment="kpoint 41">
        <r>    4.1038    1.0000 </r>
        <r>    7.3956    1.0471 </r>
        <r>    9.6030    0.0001 </r>
        <r>   14.6052    0.0000 </r>
        <r>   16.2338    0.0000 </r>
       </set>
       <set comment="kpoint 42">
        <r>    5.3203    1.0000 </r>
        <r>    5.9384    1.0000 </r>
        <r>    8.9330    0.0049 </r>
        <r>   15.1119    0.0000 </r>
        <r>   16.8439    0.0000 </r>
       </set>
       <set comment="kpoint 43">
        <r>    3.5681    1.0000 </r>
        <r>    7.5369    1.0415 </r>
        <r>    9.0263    0.0044 </r>
        <r>   15.8946    0.0000 </r>
        <r>   17.7290    0.0000 </r>
       </set>
       <set comment="kpoint 44">
        <r>    2.1225    1.0000 </r>
        <r>   10.1899    0.0000 </r>
        <r>   13.7349    0.0000 </r>
        <r>   14.0057    0.0000 </r>
        <r>   14.0057    0.0000 </r>
       </set>
       <set comment="kpoint 45">
        <r>    3.2895    1.0000 </r>
        <r>    8.7585   -0.0025 </r>
        <r>   11.6066    0.0000 </r>
        <r>   13.8699    0.0000 </r>
        <r>   14.9490    0.0000 </r>
       </set>
       <set comment="kpoint 46">
        <r>    4.7029    1.0000 </r>
        <r>    7.6607    0.9672 </r>
        <r>    9.3425    0.0007 </r>
        <r>   14.0086    0.0000 </r>
        <r>   15.4595    0.0000 </r>
       </set>
       <set comment="kpoint 47">
        <r>    5.9615    1.0000 </r>
        <r>    7.0099    0.9961 </r>
        <r>    7.5610    1.0334 </r>
        <r>   14.4224    0.0000 </r>
        <r>   16.0315    0.0000 </r>
       </set>
       <set comment="kpoint 48">
        <r>    4.3007    1.0000 </r>
        <r>    7.2695    1.0247 </r>
        <r>   13.3121    0.0000 </r>
        <r>   13.6415    0.0000 </r>
        <r>   13.6415    0.0000 </r>
       </set>
       <set comment="kpoint 49">
        <r>    5.1754    1.0000 </r>
        <r>    6.5036    0.9995 </r>
        <r>   11.4071    0.0000 </r>
        <r>   13.6877    0.0000 </r>
        <r>   14.7469    0.0000 </r>
       </set>
       <set comment="kpoint 50">
        <r>   -2.0652    1.0000 </r>
        <r>   15.0315    0.0000 </r>
        <r>   16.1512    0.0000 </r>
        <r>   20.0199    0.0000 </r>
        <r>   21.1683    0.0000 </r>
       </set>
       <set comment="kpoint 51">
        <r>   -1.2689    1.0000 </r>
        <r>   12.4574    0.0000 </r>
        <r>   15.2797    0.0000 </r>
        <r>   19.5419    0.0000 </r>
        <r>   20.1691    0.0000 </r>
       </set>
       <set comment="kpoint 52">
        <r>   -0.1614    1.0000 </r>
        <r>    9.7382    0.0000 </r>
        <r>   15.0008    0.0000 </r>
        <r>   18.0988    0.0000 </r>
        <r>   20.0874    0.0000 </r>
       </set>
       <set comment="kpoint 53">
        <r>    1.2498    1.0000 </r>
        <r>    7.2214    1.0161 </r>
        <r>   15.0354    0.0000 </r>
        <r>   16.9219    0.0000 </r>
        <r>   20.1896    0.0000 </r>
       </set>
       <set comment="kpoint 54">
        <r>    2.9531    1.0000 </r>
        <r>    4.9602    1.0000 </r>
        <r>   15.3444    0.0000 </r>
        <r>   16.0412    0.0000 </r>
        <r>   20.5287    0.0000 </r>
       </set>
       <set comment="kpoint 55">
        <r>   -1.0045    1.0000 </r>
        <r>   13.4633    0.0000 </r>
        <r>   14.6801    0.0000 </r>
        <r>   17.4148    0.0000 </r>
        <r>   18.5677    0.0000 </r>
       </set>
       <set comment="kpoint 56">
        <r>   -0.0042    1.0000 </r>
        <r>   11.0792    0.0000 </r>
        <r>   13.9496    0.0000 </r>
        <r>   15.8519    0.0000 </r>
        <r>   18.9514    0.0000 </r>
       </set>
       <set comment="kpoint 57">
        <r>    1.3008    1.0000 </r>
        <r>    8.5425   -0.0369 </r>
        <r>   13.8129    0.0000 </r>
        <r>   14.5909    0.0000 </r>
        <r>   19.0846    0.0000 </r>
       </set>
       <set comment="kpoint 58">
        <r>    2.9004    1.0000 </r>
        <r>    6.2143    1.0000 </r>
        <r>   13.4410    0.0000 </r>
        <r>   14.1743    0.0000 </r>
        <r>   19.3910    0.0000 </r>
       </set>
       <set comment="kpoint 59">
        <r>    4.1120    1.0000 </r>
        <r>    4.8161    1.0000 </r>
        <r>   12.7617    0.0000 </r>
        <r>   14.6560    0.0000 </r>
        <r>   19.9480    0.0000 </r>
       </set>
       <set comment="kpoint 60">
        <r>    2.3367    1.0000 </r>
        <r>    6.9472    0.9951 </r>
        <r>   12.3544    0.0000 </r>
        <r>   15.4615    0.0000 </r>
        <r>   20.7617    0.0000 </r>
       </set>
       <set comment="kpoint 61">
        <r>    0.8329    1.0000 </r>
        <r>    9.3076    0.0009 </r>
        <r>   12.2910    0.0000 </r>
        <r>   16.5547    0.0000 </r>
        <r>   21.7689    0.0000 </r>
       </set>
       <set comment="kpoint 62">
        <r>   -0.3717    1.0000 </r>
        <r>   11.4872    0.0000 </r>
        <r>   12.9443    0.0000 </r>
        <r>   17.9221    0.0000 </r>
        <r>   21.6206    0.0000 </r>
       </set>
       <set comment="kpoint 63">
        <r>    0.4665    1.0000 </r>
        <r>   12.2802    0.0000 </r>
        <r>   13.4000    0.0000 </r>
        <r>   13.9191    0.0000 </r>
        <r>   17.4206    0.0000 </r>
       </set>
       <set comment="kpoint 64">
        <r>    1.6620    1.0000 </r>
        <r>   10.0835    0.0000 </r>
        <r>   12.2632    0.0000 </r>
        <r>   13.1030    0.0000 </r>
        <r>   18.0118    0.0000 </r>
       </set>
       <set comment="kpoint 65">
        <r>    3.1512    1.0000 </r>
        <r>    7.7385    0.8741 </r>
        <r>   11.1916    0.0000 </r>
        <r>   13.1430    0.0000 </r>
        <r>   18.3471    0.0000 </r>
       </set>
       <set comment="kpoint 66">
        <r>    4.8972    1.0000 </r>
        <r>    5.6246    1.0000 </r>
        <r>   10.3862    0.0000 </r>
        <r>   13.5521    0.0000 </r>
        <r>   18.8516    0.0000 </r>
       </set>
       <set comment="kpoint 67">
        <r>    3.7194    1.0000 </r>
        <r>    6.9194    0.9951 </r>
        <r>    9.9299    0.0000 </r>
        <r>   14.2668    0.0000 </r>
        <r>   19.6050    0.0000 </r>
       </set>
       <set comment="kpoint 68">
        <r>    2.1301    1.0000 </r>
        <r>    8.7209   -0.0065 </r>
        <r>   10.2275    0.0000 </r>
        <r>   15.2714    0.0000 </r>
        <r>   20.6107    0.0000 </r>
       </set>
       <set comment="kpoint 69">
        <r>    2.3288    1.0000 </r>
        <r>   10.3743    0.0000 </r>
        <r>   11.4891    0.0000 </r>
        <r>   12.8494    0.0000 </r>
        <r>   16.6458    0.0000 </r>
       </set>
       <set comment="kpoint 70">
        <r>    3.6976    1.0000 </r>
        <r>    9.0878    0.0035 </r>
        <r>    9.5783    0.0001 </r>
        <r>   12.4683    0.0000 </r>
        <r>   17.4349    0.0000 </r>
       </set>
       <set comment="kpoint 71">
        <r>    5.3026    1.0000 </r>
        <r>    7.3024    1.0310 </r>
        <r>    8.3514   -0.0447 </r>
        <r>   12.7441    0.0000 </r>
        <r>   17.9677    0.0000 </r>
       </set>
       <set comment="kpoint 72">
        <r>    5.3813    1.0000 </r>
        <r>    6.7487    0.9971 </r>
        <r>    8.2121    0.0359 </r>
        <r>   13.3561    0.0000 </r>
        <r>   18.6664    0.0000 </r>
       </set>
       <set comment="kpoint 73">
        <r>    4.5020    1.0000 </r>
        <r>    7.4705    1.0514 </r>
        <r>   11.0930    0.0000 </r>
        <r>   12.4574    0.0000 </r>
        <r>   16.2568    0.0000 </r>
       </set>
       <set comment="kpoint 74">
        <r>    5.5683    1.0000 </r>
        <r>    6.9102    0.9951 </r>
        <r>    9.3209    0.0008 </r>
        <r>   12.2701    0.0000 </r>
        <r>   17.2411    0.0000 </r>
       </set>
       <set comment="kpoint 75">
        <r>    1.0927    1.0000 </r>
        <r>   10.3020    0.0000 </r>
        <r>   11.7654    0.0000 </r>
        <r>   14.3713    0.0000 </r>
        <r>   20.5530    0.0000 </r>
       </set>
       <set comment="kpoint 76">
        <r>    2.4891    1.0000 </r>
        <r>    8.3117   -0.0321 </r>
        <r>   11.2845    0.0000 </r>
        <r>   13.1729    0.0000 </r>
        <r>   21.2332    0.0000 </r>
       </set>
       <set comment="kpoint 77">
        <r>    4.1696    1.0000 </r>
        <r>    6.1448    1.0000 </r>
        <r>   11.4850    0.0000 </r>
        <r>   12.3146    0.0000 </r>
        <r>   21.6515    0.0000 </r>
       </set>
       <set comment="kpoint 78">
        <r>    2.9462    1.0000 </r>
        <r>    9.5115    0.0001 </r>
        <r>   10.6699    0.0000 </r>
        <r>   11.3690    0.0000 </r>
        <r>   19.7986    0.0000 </r>
       </set>
       <set comment="kpoint 79">
        <r>    4.5096    1.0000 </r>
        <r>    7.7074    0.9158 </r>
        <r>    9.8927    0.0000 </r>
        <r>   10.8643    0.0000 </r>
        <r>   20.6961    0.0000 </r>
       </set>
       <set comment="kpoint 80">
        <r>    5.7016    1.0000 </r>
        <r>    6.3354    0.9999 </r>
        <r>    9.2983    0.0010 </r>
        <r>   11.2554    0.0000 </r>
        <r>   21.3526    0.0000 </r>
       </set>
       <set comment="kpoint 81">
        <r>    3.9686    1.0000 </r>
        <r>    7.9315    0.5138 </r>
        <r>    9.4006    0.0004 </r>
        <r>   12.0547    0.0000 </r>
        <r>   22.1590    0.0000 </r>
       </set>
       <set comment="kpoint 82">
        <r>    5.1044    1.0000 </r>
        <r>    8.0578    0.2576 </r>
        <r>    9.1163    0.0031 </r>
        <r>   10.6891    0.0000 </r>
        <r>   19.4155    0.0000 </r>
       </set>
       <set comment="kpoint 83">
        <r>    6.3550    0.9999 </r>
        <r>    7.3825    1.0453 </r>
        <r>    7.8767    0.6283 </r>
        <r>   10.5937    0.0000 </r>
        <r>   20.5118    0.0000 </r>
       </set>
       <set comment="kpoint 84">
        <r>    6.1012    1.0000 </r>
        <r>    7.5371    1.0415 </r>
        <r>    8.6857   -0.0113 </r>
        <r>    9.5044    0.0002 </r>
        <r>   22.9474    0.0000 </r>
       </set>
      </set>
      <set comment="spin 2">
       <set comment="kpoint 1">
        <r>   -3.1341    1.0000 </r>
        <r>   20.6208    0.0000 </r>
        <r>   20.6208    0.0000 </r>
        <r>   20.6208    0.0000 </r>
        <r>   22.1131    0.0000 </r>
       </set>
       <set comment="kpoint 2">
        <r>   -2.9732    1.0000 </r>
        <r>   18.2625    0.0000 </r>
        <r>   20.1120    0.0000 </r>
        <r>   20.1120    0.0000 </r>
        <r>   22.7714    0.0000 </r>
       </set>
       <set comment="kpoint 3">
        <r>   -2.4918    1.0000 </r>
        <r>   15.1545    0.0000 </r>
        <r>   19.4447    0.0000 </r>
        <r>   19.4447    0.0000 </r>
        <r>   22.6316    0.0000 </r>
       </set>
       <set comment="kpoint 4">
        <r>   -1.6930    1.0000 </r>
        <r>   12.1586    0.0000 </r>
        <r>   18.9854    0.0000 </r>
        <r>   18.9854    0.0000 </r>
        <r>   22.3980    0.0000 </r>
       </set>
       <set comment="kpoint 5">
        <r>   -0.5822    1.0000 </r>
        <r>    9.3722    0.0005 </r>
        <r>   18.7803    0.0000 </r>
        <r>   18.7805    0.0000 </r>
        <r>   21.7676    0.0000 </r>
       </set>
       <set comment="kpoint 6">
        <r>    0.8335    1.0000 </r>
        <r>    6.8308    0.9959 </r>
        <r>   18.8406    0.0000 </r>
        <r>   18.8412    0.0000 </r>
        <r>   20.6402    0.0000 </r>
       </set>
       <set comment="kpoint 7">
        <r>    2.5434    1.0000 </r>
        <r>    4.5561    1.0000 </r>
        <r>   19.1696    0.0000 </r>
        <r>   19.1741    0.0000 </r>
        <r>   19.7719    0.0000 </r>
       </set>
       <set comment="kpoint 8">
        <r>   -2.9196    1.0000 </r>
        <r>   19.0625    0.0000 </r>
        <r>   19.0625    0.0000 </r>
        <r>   19.3530    0.0000 </r>
        <r>   21.6182    0.0000 </r>
       </set>
       <set comment="kpoint 9">
        <r>   -2.5452    1.0000 </r>
        <r>   16.2579    0.0000 </r>
        <r>   18.2305    0.0000 </r>
        <r>   19.4493    0.0000 </r>
        <r>   21.9398    0.0000 </r>
       </set>
       <set comment="kpoint 10">
        <r>   -1.8524    1.0000 </r>
        <r>   13.2453    0.0000 </r>
        <r>   17.6567    0.0000 </r>
        <r>   19.1235    0.0000 </r>
        <r>   21.3428    0.0000 </r>
       </set>
       <set comment="kpoint 11">
        <r>   -0.8459    1.0000 </r>
        <r>   10.3980    0.0000 </r>
        <r>   17.3504    0.0000 </r>
        <r>   18.8807    0.0000 </r>
        <r>   19.8568    0.0000 </r>
       </set>
       <set comment="kpoint 12">
        <r>    0.4677    1.0000 </r>
        <r>    7.7792    0.8110 </r>
        <r>   17.3148    0.0000 </r>
        <r>   18.5785    0.0000 </r>
        <r>   18.9007    0.0000 </r>
       </set>
       <set comment="kpoint 13">
        <r>    2.0796    1.0000 </r>
        <r>    5.4172    1.0000 </r>
        <r>   17.5514    0.0000 </r>
        <r>   17.6063    0.0000 </r>
        <r>   19.1519    0.0000 </r>
       </set>
       <set comment="kpoint 14">
        <r>    3.3035    1.0000 </r>
        <r>    4.0052    1.0000 </r>
        <r>   16.8983    0.0000 </r>
        <r>   18.0612    0.0000 </r>
        <r>   19.6931    0.0000 </r>
       </set>
       <set comment="kpoint 15">
        <r>    1.5100    1.0000 </r>
        <r>    6.1739    1.0000 </r>
        <r>   16.4703    0.0000 </r>
        <r>   18.8428    0.0000 </r>
        <r>   20.5099    0.0000 </r>
       </set>
       <set comment="kpoint 16">
        <r>   -0.0038    1.0000 </r>
        <r>    8.6224   -0.0218 </r>
        <r>   16.3249    0.0000 </r>
        <r>   19.8771    0.0000 </r>
        <r>   20.5167    0.0000 </r>
       </set>
       <set comment="kpoint 17">
        <r>   -1.2159    1.0000 </r>
        <r>   11.3163    0.0000 </r>
        <r>   16.4641    0.0000 </r>
        <r>   20.5551    0.0000 </r>
        <r>   21.2451    0.0000 </r>
       </set>
       <set comment="kpoint 18">
        <r>   -2.1185    1.0000 </r>
        <r>   14.2062    0.0000 </r>
        <r>   16.8998    0.0000 </r>
        <r>   20.8874    0.0000 </r>
        <r>   22.6175    0.0000 </r>
       </set>
       <set comment="kpoint 19">
        <r>   -2.7056    1.0000 </r>
        <r>   16.9661    0.0000 </r>
        <r>   17.8588    0.0000 </r>
        <r>   21.0378    0.0000 </r>
        <r>   22.3631    0.0000 </r>
       </set>
       <set comment="kpoint 20">
        <r>   -2.2783    1.0000 </r>
        <r>   17.2814    0.0000 </r>
        <r>   17.2814    0.0000 </r>
        <r>   17.4619    0.0000 </r>
        <r>   20.5135    0.0000 </r>
       </set>
       <set comment="kpoint 21">
        <r>   -1.6930    1.0000 </r>
        <r>   14.5332    0.0000 </r>
        <r>   16.6041    0.0000 </r>
        <r>   17.8754    0.0000 </r>
        <r>   19.2377    0.0000 </r>
       </set>
       <set comment="kpoint 22">
        <r>   -0.7932    1.0000 </r>
        <r>   11.6813    0.0000 </r>
        <r>   16.2000    0.0000 </r>
        <r>   17.5061    0.0000 </r>
        <r>   17.8750    0.0000 </r>
       </set>
       <set comment="kpoint 23">
        <r>    0.4150    1.0000 </r>
        <r>    9.0083    0.0046 </r>
        <r>   16.0694    0.0000 </r>
        <r>   16.1799    0.0000 </r>
        <r>   17.7810    0.0000 </r>
       </set>
       <set comment="kpoint 24">
        <r>    1.9234    1.0000 </r>
        <r>    6.5725    0.9991 </r>
        <r>   15.0818    0.0000 </r>
        <r>   16.2127    0.0000 </r>
        <r>   17.9735    0.0000 </r>
       </set>
       <set comment="kpoint 25">
        <r>    3.7068    1.0000 </r>
        <r>    4.4152    1.0000 </r>
        <r>   14.2658    0.0000 </r>
        <r>   16.6303    0.0000 </r>
        <r>   18.4467    0.0000 </r>
       </set>
       <set comment="kpoint 26">
        <r>    2.4919    1.0000 </r>
        <r>    5.8125    1.0000 </r>
        <r>   13.7414    0.0000 </r>
        <r>   17.3224    0.0000 </r>
        <r>   19.2038    0.0000 </r>
       </set>
       <set comment="kpoint 27">
        <r>    0.8854    1.0000 </r>
        <r>    8.1494    0.1093 </r>
        <r>   13.5200    0.0000 </r>
        <r>   18.2870    0.0000 </r>
        <r>   20.2350    0.0000 </r>
       </set>
       <set comment="kpoint 28">
        <r>   -0.4243    1.0000 </r>
        <r>   10.6937    0.0000 </r>
        <r>   13.6386    0.0000 </r>
        <r>   19.5157    0.0000 </r>
        <r>   21.5173    0.0000 </r>
       </set>
       <set comment="kpoint 29">
        <r>   -1.4278    1.0000 </r>
        <r>   13.0656    0.0000 </r>
        <r>   14.4433    0.0000 </r>
        <r>   20.9757    0.0000 </r>
        <r>   22.2850    0.0000 </r>
       </set>
       <set comment="kpoint 30">
        <r>   -1.2160    1.0000 </r>
        <r>   15.8284    0.0000 </r>
        <r>   15.8284    0.0000 </r>
        <r>   15.8293    0.0000 </r>
        <r>   17.2063    0.0000 </r>
       </set>
       <set comment="kpoint 31">
        <r>   -0.4246    1.0000 </r>
        <r>   13.1660    0.0000 </r>
        <r>   15.3287    0.0000 </r>
        <r>   15.4253    0.0000 </r>
        <r>   16.5538    0.0000 </r>
       </set>
       <set comment="kpoint 32">
        <r>    0.6755    1.0000 </r>
        <r>   10.4935    0.0000 </r>
        <r>   13.9343    0.0000 </r>
        <r>   15.1040    0.0000 </r>
        <r>   16.6700    0.0000 </r>
       </set>
       <set comment="kpoint 33">
        <r>    2.0762    1.0000 </r>
        <r>    8.0069    0.3565 </r>
        <r>   12.7214    0.0000 </r>
        <r>   15.1539    0.0000 </r>
        <r>   16.8304    0.0000 </r>
       </set>
       <set comment="kpoint 34">
        <r>    3.7624    1.0000 </r>
        <r>    5.7667    1.0000 </r>
        <r>   11.8000    0.0000 </r>
        <r>   15.4786    0.0000 </r>
        <r>   17.2330    0.0000 </r>
       </set>
       <set comment="kpoint 35">
        <r>    3.7703    1.0000 </r>
        <r>    5.7430    1.0000 </r>
        <r>   11.1866    0.0000 </r>
        <r>   16.0784    0.0000 </r>
        <r>   17.9154    0.0000 </r>
       </set>
       <set comment="kpoint 36">
        <r>    2.0791    1.0000 </r>
        <r>    7.9166    0.5452 </r>
        <r>   10.9288    0.0000 </r>
        <r>   16.9535    0.0000 </r>
        <r>   18.8771    0.0000 </r>
       </set>
       <set comment="kpoint 37">
        <r>    0.6766    1.0000 </r>
        <r>    9.9061    0.0000 </r>
        <r>   11.4181    0.0000 </r>
        <r>   18.1019    0.0000 </r>
        <r>   20.0968    0.0000 </r>
       </set>
       <set comment="kpoint 38">
        <r>    0.2572    1.0000 </r>
        <r>   13.5342    0.0000 </r>
        <r>   14.5772    0.0000 </r>
        <r>   14.7348    0.0000 </r>
        <r>   14.7348    0.0000 </r>
       </set>
       <set comment="kpoint 39">
        <r>    1.2467    1.0000 </r>
        <r>   11.8647    0.0000 </r>
        <r>   12.2444    0.0000 </r>
        <r>   14.4167    0.0000 </r>
        <r>   15.5526    0.0000 </r>
       </set>
       <set comment="kpoint 40">
        <r>    2.5349    1.0000 </r>
        <r>    9.6698    0.0000 </r>
        <r>   10.6384    0.0000 </r>
        <r>   14.3736    0.0000 </r>
        <r>   15.8675    0.0000 </r>
       </set>
       <set comment="kpoint 41">
        <r>    4.1038    1.0000 </r>
        <r>    7.3956    1.0471 </r>
        <r>    9.6030    0.0001 </r>
        <r>   14.6052    0.0000 </r>
        <r>   16.2338    0.0000 </r>
       </set>
       <set comment="kpoint 42">
        <r>    5.3203    1.0000 </r>
        <r>    5.9384    1.0000 </r>
        <r>    8.9330    0.0049 </r>
        <r>   15.1119    0.0000 </r>
        <r>   16.8439    0.0000 </r>
       </set>
       <set comment="kpoint 43">
        <r>    3.5681    1.0000 </r>
        <r>    7.5369    1.0415 </r>
        <r>    9.0263    0.0044 </r>
        <r>   15.8946    0.0000 </r>
        <r>   17.7290    0.0000 </r>
       </set>
       <set comment="kpoint 44">
        <r>    2.1225    1.0000 </r>
        <r>   10.1899    0.0000 </r>
        <r>   13.7349    0.0000 </r>
        <r>   14.0057    0.0000 </r>
        <r>   14.0057    0.0000 </r>
       </set>
       <set comment="kpoint 45">
        <r>    3.2895    1.0000 </r>
        <r>    8.7585   -0.0025 </r>
        <r>   11.6066    0.0000 </r>
        <r>   13.8699    0.0000 </r>
        <r>   14.9490    0.0000 </r>
       </set>
       <set comment="kpoint 46">
        <r>    4.7029    1.0000 </r>
        <r>    7.6607    0.9672 </r>
        <r>    9.3425    0.0007 </r>
        <r>   14.0086    0.0000 </r>
        <r>   15.4595    0.0000 </r>
       </set>
       <set comment="kpoint 47">
        <r>    5.9615    1.0000 </r>
        <r>    7.0099    0.9961 </r>
        <r>    7.5610    1.0334 </r>
        <r>   14.4224    0.0000 </r>
        <r>   16.0315    0.0000 </r>
       </set>
       <set comment="kpoint 48">
        <r>    4.3007    1.0000 </r>
        <r>    7.2695    1.0247 </r>
        <r>   13.3121    0.0000 </r>
        <r>   13.6415    0.0000 </r>
        <r>   13.6415    0.0000 </r>
       </set>
       <set comment="kpoint 49">
        <r>    5.1754    1.0000 </r>
        <r>    6.5036    0.9995 </r>
        <r>   11.4071    0.0000 </r>
        <r>   13.6877    0.0000 </r>
        <r>   14.7469    0.0000 </r>
       </set>
       <set comment="kpoint 50">
        <r>   -2.0652    1.0000 </r>
        <r>   15.0315    0.0000 </r>
        <r>   16.1512    0.0000 </r>
        <r>   20.0199    0.0000 </r>
        <r>   21.1664    0.0000 </r>
       </set>
       <set comment="kpoint 51">
        <r>   -1.2689    1.0000 </r>
        <r>   12.4574    0.0000 </r>
        <r>   15.2797    0.0000 </r>
        <r>   19.5420    0.0000 </r>
        <r>   20.1691    0.0000 </r>
       </set>
       <set comment="kpoint 52">
        <r>   -0.1614    1.0000 </r>
        <r>    9.7382    0.0000 </r>
        <r>   15.0008    0.0000 </r>
        <r>   18.0988    0.0000 </r>
        <r>   20.0874    0.0000 </r>
       </set>
       <set comment="kpoint 53">
        <r>    1.2498    1.0000 </r>
        <r>    7.2214    1.0161 </r>
        <r>   15.0354    0.0000 </r>
        <r>   16.9219    0.0000 </r>
        <r>   20.1896    0.0000 </r>
       </set>
       <set comment="kpoint 54">
        <r>    2.9531    1.0000 </r>
        <r>    4.9602    1.0000 </r>
        <r>   15.3444    0.0000 </r>
        <r>   16.0412    0.0000 </r>
        <r>   20.5273    0.0000 </r>
       </set>
       <set comment="kpoint 55">
        <r>   -1.0045    1.0000 </r>
        <r>   13.4633    0.0000 </r>
        <r>   14.6801    0.0000 </r>
        <r>   17.4148    0.0000 </r>
        <r>   18.5677    0.0000 </r>
       </set>
       <set comment="kpoint 56">
        <r>   -0.0042    1.0000 </r>
        <r>   11.0792    0.0000 </r>
        <r>   13.9496    0.0000 </r>
        <r>   15.8519    0.0000 </r>
        <r>   18.9514    0.0000 </r>
       </set>
       <set comment="kpoint 57">
        <r>    1.3008    1.0000 </r>
        <r>    8.5425   -0.0369 </r>
        <r>   13.8129    0.0000 </r>
        <r>   14.5909    0.0000 </r>
        <r>   19.0846    0.0000 </r>
       </set>
       <set comment="kpoint 58">
        <r>    2.9004    1.0000 </r>
        <r>    6.2143    1.0000 </r>
        <r>   13.4410    0.0000 </r>
        <r>   14.1743    0.0000 </r>
        <r>   19.3910    0.0000 </r>
       </set>
       <set comment="kpoint 59">
        <r>    4.1120    1.0000 </r>
        <r>    4.8161    1.0000 </r>
        <r>   12.7617    0.0000 </r>
        <r>   14.6560    0.0000 </r>
        <r>   19.9480    0.0000 </r>
       </set>
       <set comment="kpoint 60">
        <r>    2.3367    1.0000 </r>
        <r>    6.9472    0.9951 </r>
        <r>   12.3544    0.0000 </r>
        <r>   15.4615    0.0000 </r>
        <r>   20.7617    0.0000 </r>
       </set>
       <set comment="kpoint 61">
        <r>    0.8329    1.0000 </r>
        <r>    9.3076    0.0009 </r>
        <r>   12.2910    0.0000 </r>
        <r>   16.5547    0.0000 </r>
        <r>   21.7678    0.0000 </r>
       </set>
       <set comment="kpoint 62">
        <r>   -0.3717    1.0000 </r>
        <r>   11.4872    0.0000 </r>
        <r>   12.9443    0.0000 </r>
        <r>   17.9221    0.0000 </r>
        <r>   21.6206    0.0000 </r>
       </set>
       <set comment="kpoint 63">
        <r>    0.4665    1.0000 </r>
        <r>   12.2802    0.0000 </r>
        <r>   13.4000    0.0000 </r>
        <r>   13.9191    0.0000 </r>
        <r>   17.4206    0.0000 </r>
       </set>
       <set comment="kpoint 64">
        <r>    1.6620    1.0000 </r>
        <r>   10.0835    0.0000 </r>
        <r>   12.2632    0.0000 </r>
        <r>   13.1030    0.0000 </r>
        <r>   18.0118    0.0000 </r>
       </set>
       <set comment="kpoint 65">
        <r>    3.1512    1.0000 </r>
        <r>    7.7385    0.8741 </r>
        <r>   11.1916    0.0000 </r>
        <r>   13.1430    0.0000 </r>
        <r>   18.3471    0.0000 </r>
       </set>
       <set comment="kpoint 66">
        <r>    4.8972    1.0000 </r>
        <r>    5.6246    1.0000 </r>
        <r>   10.3862    0.0000 </r>
        <r>   13.5521    0.0000 </r>
        <r>   18.8516    0.0000 </r>
       </set>
       <set comment="kpoint 67">
        <r>    3.7194    1.0000 </r>
        <r>    6.9194    0.9951 </r>
        <r>    9.9299    0.0000 </r>
        <r>   14.2668    0.0000 </r>
        <r>   19.6050    0.0000 </r>
       </set>
       <set comment="kpoint 68">
        <r>    2.1301    1.0000 </r>
        <r>    8.7209   -0.0065 </r>
        <r>   10.2275    0.0000 </r>
        <r>   15.2714    0.0000 </r>
        <r>   20.6107    0.0000 </r>
       </set>
       <set comment="kpoint 69">
        <r>    2.3288    1.0000 </r>
        <r>   10.3743    0.0000 </r>
        <r>   11.4891    0.0000 </r>
        <r>   12.8494    0.0000 </r>
        <r>   16.6458    0.0000 </r>
       </set>
       <set comment="kpoint 70">
        <r>    3.6976    1.0000 </r>
        <r>    9.0878    0.0035 </r>
        <r>    9.5783    0.0001 </r>
        <r>   12.4683    0.0000 </r>
        <r>   17.4349    0.0000 </r>
       </set>
       <set comment="kpoint 71">
        <r>    5.3026    1.0000 </r>
        <r>    7.3024    1.0310 </r>
        <r>    8.3514   -0.0447 </r>
        <r>   12.7441    0.0000 </r>
        <r>   17.9677    0.0000 </r>
       </set>
       <set comment="kpoint 72">
        <r>    5.3813    1.0000 </r>
        <r>    6.7487    0.9971 </r>
        <r>    8.2121    0.0359 </r>
        <r>   13.3561    0.0000 </r>
        <r>   18.6664    0.0000 </r>
       </set>
       <set comment="kpoint 73">
        <r>    4.5020    1.0000 </r>
        <r>    7.4705    1.0514 </r>
        <r>   11.0930    0.0000 </r>
        <r>   12.4574    0.0000 </r>
        <r>   16.2568    0.0000 </r>
       </set>
       <set comment="kpoint 74">
        <r>    5.5683    1.0000 </r>
        <r>    6.9102    0.9951 </r>
        <r>    9.3209    0.0008 </r>
        <r>   12.2701    0.0000 </r>
        <r>   17.2411    0.0000 </r>
       </set>
       <set comment="kpoint 75">
        <r>    1.0927    1.0000 </r>
        <r>   10.3020    0.0000 </r>
        <r>   11.7654    0.0000 </r>
        <r>   14.3713    0.0000 </r>
        <r>   20.5530    0.0000 </r>
       </set>
       <set comment="kpoint 76">
        <r>    2.4891    1.0000 </r>
        <r>    8.3117   -0.0321 </r>
        <r>   11.2845    0.0000 </r>
        <r>   13.1729    0.0000 </r>
        <r>   21.2332    0.0000 </r>
       </set>
       <set comment="kpoint 77">
        <r>    4.1696    1.0000 </r>
        <r>    6.1448    1.0000 </r>
        <r>   11.4850    0.0000 </r>
        <r>   12.3146    0.0000 </r>
        <r>   21.6506    0.0000 </r>
       </set>
       <set comment="kpoint 78">
        <r>    2.9462    1.0000 </r>
        <r>    9.5115    0.0001 </r>
        <r>   10.6699    0.0000 </r>
        <r>   11.3690    0.0000 </r>
        <r>   19.7985    0.0000 </r>
       </set>
       <set comment="kpoint 79">
        <r>    4.5096    1.0000 </r>
        <r>    7.7074    0.9158 </r>
        <r>    9.8927    0.0000 </r>
        <r>   10.8643    0.0000 </r>
        <r>   20.6961    0.0000 </r>
       </set>
       <set comment="kpoint 80">
        <r>    5.7016    1.0000 </r>
        <r>    6.3354    0.9999 </r>
        <r>    9.2983    0.0010 </r>
        <r>   11.2554    0.0000 </r>
        <r>   21.3526    0.0000 </r>
       </set>
       <set comment="kpoint 81">
        <r>    3.9686    1.0000 </r>
        <r>    7.9315    0.5138 </r>
        <r>    9.4006    0.0004 </r>
        <r>   12.0547    0.0000 </r>
        <r>   22.1590    0.0000 </r>
       </set>
       <set comment="kpoint 82">
        <r>    5.1044    1.0000 </r>
        <r>    8.0578    0.2576 </r>
        <r>    9.1163    0.0031 </r>
        <r>   10.6891    0.0000 </r>
        <r>   19.4155    0.0000 </r>
       </set>
       <set comment="kpoint 83">
        <r>    6.3550    0.9999 </r>
        <r>    7.3825    1.0453 </r>
        <r>    7.8767    0.6283 </r>
        <r>   10.5937    0.0000 </r>
        <r>   20.5118    0.0000 </r>
       </set>
       <set comment="kpoint 84">
        <r>    6.1012    1.0000 </r>
        <r>    7.5371    1.0415 </r>
        <r>    8.6857   -0.0113 </r>
        <r>    9.5044    0.0002 </r>
        <r>   22.9401    0.0000 </r>
       </set>
      </set>
     </set>
    </array>
   </eigenvalues>
   <array>
    <dimension dim="1">ion</dimension>
    <dimension dim="2">band</dimension>
    <dimension dim="3">kpoint</dimension>
    <dimension dim="4">spin</dimension>
    <field>    s</field>
    <field>   py</field>
    <field>   pz</field>
    <field>   px</field>
    <field>  dxy</field>
    <field>  dyz</field>
    <field>  dz2</field>
    <field>  dxz</field>
    <field>x2-y2</field>
    <set>
     <set comment="spin1">
      <set comment="kpoint 1">
       <set comment="band 1">
        <r>  0.2701  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.3415  0.0196  0.1553  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 2">
       <set comment="band 1">
        <r>  0.2669  0.0011  0.0011  0.0011  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0073  0.0501  0.0500  0.0501  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0013  0.0690  0.0574  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0837  0.0159  0.0278  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0233  0.1008  0.0484  0.1429  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 3">
       <set comment="band 1">
        <r>  0.2574  0.0043  0.0043  0.0043  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0266  0.0565  0.0565  0.0565  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1324  0.0413  0.0284  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0024  0.0934  0.1063  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0676  0.0527  0.0518  0.0533  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 4">
       <set comment="band 1">
        <r>  0.2418  0.0094  0.0094  0.0094  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0521  0.0559  0.0559  0.0559  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0615  0.0230  0.1501  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0948  0.1334  0.0063  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1033  0.0262  0.0288  0.0213  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 5">
       <set comment="band 1">
        <r>  0.2205  0.0162  0.0162  0.0162  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0818  0.0520  0.0520  0.0520  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0496  0.0349  0.1662  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1165  0.1310  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0472  0.0221  0.1315  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 6">
       <set comment="band 1">
        <r>  0.1935  0.0244  0.0244  0.0244  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1147  0.0456  0.0456  0.0456  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0032  0.1280  0.1249  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1628  0.0415  0.0445  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0027  0.1070  0.1033  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 7">
       <set comment="band 1">
        <r>  0.1543  0.0352  0.0352  0.0352  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1570  0.0357  0.0357  0.0357  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0983  0.0076  0.1569  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0607  0.1411  0.0173  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0453  0.0295  0.1413  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 8">
       <set comment="band 1">
        <r>  0.2659  0.0000  0.0044  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.1037  0.0000  0.0963  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0963  0.0000  0.1037  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0449  0.0000  0.4510  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 9">
       <set comment="band 1">
        <r>  0.2584  0.0011  0.0097  0.0011  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0141  0.0712  0.0227  0.0712  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1182  0.0000  0.1182  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0214  0.0086  0.0751  0.0086  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0429  0.0469  0.3301  0.0469  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 10">
       <set comment="band 1">
        <r>  0.2448  0.0042  0.0170  0.0042  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0388  0.0691  0.0322  0.0691  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1288  0.0000  0.1288  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0460  0.0001  0.1196  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0129  0.0126  0.2021  0.0126  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 11">
       <set comment="band 1">
        <r>  0.2254  0.0092  0.0258  0.0092  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0680  0.0646  0.0332  0.0646  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1345  0.0000  0.1345  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0544  0.0045  0.1049  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0204  0.0132  0.1992  0.0132  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 12">
       <set comment="band 1">
        <r>  0.2004  0.0157  0.0359  0.0157  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1002  0.0579  0.0300  0.0579  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1365  0.0000  0.1365  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0065  0.0082  0.1029  0.0082  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0785  0.0224  0.1821  0.0224  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 13">
       <set comment="band 1">
        <r>  0.1676  0.0239  0.0473  0.0239  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1373  0.0488  0.0238  0.0488  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1355  0.0000  0.1355  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0234  0.0187  0.1540  0.0187  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0657  0.0284  0.1062  0.0283  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 14">
       <set comment="band 1">
        <r>  0.1098  0.0530  0.0314  0.0530  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1961  0.0181  0.0440  0.0181  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0299  0.0264  0.1504  0.0264  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1327  0.0000  0.1327  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0573  0.0392  0.0813  0.0392  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 15">
       <set comment="band 1">
        <r>  0.1799  0.0348  0.0153  0.0348  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1237  0.0336  0.0638  0.0336  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0339  0.0353  0.1383  0.0353  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1306  0.0000  0.1301  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0453  0.0639  0.0621  0.0366  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 16">
       <set comment="band 1">
        <r>  0.2094  0.0251  0.0089  0.0251  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0885  0.0393  0.0740  0.0393  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0357  0.0456  0.1205  0.0457  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1568  0.0000  0.1573  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0651  0.0000  0.0664  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 17">
       <set comment="band 1">
        <r>  0.2326  0.0166  0.0041  0.0166  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0565  0.0410  0.0845  0.0410  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0355  0.0581  0.0967  0.0581  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1420  0.0000  0.1418  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0667  0.0000  0.0665  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 18">
       <set comment="band 1">
        <r>  0.2501  0.0096  0.0011  0.0096  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0262  0.0341  0.1010  0.0341  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0342  0.0761  0.0619  0.0761  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0794  0.0000  0.0794  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.1529  0.0001  0.1529  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 19">
       <set comment="band 1">
        <r>  0.2616  0.0043  0.0000  0.0043  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1429  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0300  0.1127  0.0000  0.1127  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.2705  0.0001  0.2743  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 20">
       <set comment="band 1">
        <r>  0.2531  0.0000  0.0172  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.1622  0.0000  0.0989  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0989  0.0000  0.1622  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1289  0.0000  0.3083  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 21">
       <set comment="band 1">
        <r>  0.2416  0.0010  0.0266  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0180  0.0793  0.0116  0.0793  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1393  0.0000  0.1393  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0325  0.0238  0.0293  0.0238  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0301  0.0003  0.1870  0.0003  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 22">
       <set comment="band 1">
        <r>  0.2241  0.0041  0.0374  0.0041  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0472  0.0763  0.0172  0.0763  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1448  0.0000  0.1448  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0017  0.0045  0.0647  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1218  0.0001  0.1836  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 23">
       <set comment="band 1">
        <r>  0.2008  0.0088  0.0494  0.0088  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0797  0.0702  0.0170  0.0702  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1472  0.0000  0.1472  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0349  0.0058  0.1479  0.0058  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1204  0.0019  0.0884  0.0019  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 24">
       <set comment="band 1">
        <r>  0.1710  0.0151  0.0625  0.0151  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1153  0.0616  0.0135  0.0616  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0506  0.0104  0.1490  0.0104  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1470  0.0000  0.1470  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1207  0.0099  0.0614  0.0099  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 25">
       <set comment="band 1">
        <r>  0.0996  0.0314  0.0815  0.0314  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1887  0.0421  0.0032  0.0421  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0622  0.0174  0.1381  0.0174  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1446  0.0000  0.1446  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1139  0.0215  0.0411  0.0215  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 26">
       <set comment="band 1">
        <r>  0.1592  0.0466  0.0087  0.0466  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1267  0.0223  0.0851  0.0223  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0711  0.0269  0.1203  0.0269  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1409  0.0000  0.1409  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1007  0.0349  0.0251  0.0348  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 27">
       <set comment="band 1">
        <r>  0.1922  0.0354  0.0038  0.0354  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0859  0.0266  0.1018  0.0266  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0785  0.0400  0.0956  0.0400  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1372  0.0000  0.1372  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0824  0.0481  0.0131  0.0489  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 28">
       <set comment="band 1">
        <r>  0.2172  0.0254  0.0010  0.0254  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0438  0.0224  0.1274  0.0224  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0879  0.0609  0.0578  0.0609  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1363  0.0000  0.1364  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0609  0.0567  0.0045  0.0649  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 29">
       <set comment="band 1">
        <r>  0.2366  0.0167  0.0000  0.0167  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1712  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0975  0.0991  0.0000  0.0991  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1489  0.0000  0.1489  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0001  0.0000  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 30">
       <set comment="band 1">
        <r>  0.2320  0.0000  0.0379  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0606  0.0000  0.2346  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.2346  0.0000  0.0606  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0263  0.0000  0.1465  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 31">
       <set comment="band 1">
        <r>  0.2164  0.0010  0.0508  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0190  0.0840  0.0040  0.0840  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1533  0.0000  0.1533  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0384  0.0008  0.1466  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0617  0.0304  0.0400  0.0304  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 32">
       <set comment="band 1">
        <r>  0.1949  0.0038  0.0648  0.0038  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0512  0.0814  0.0062  0.0814  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0609  0.0017  0.1521  0.0017  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1562  0.0000  0.1562  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1233  0.0030  0.0468  0.0030  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 33">
       <set comment="band 1">
        <r>  0.1669  0.0083  0.0800  0.0083  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0869  0.0742  0.0054  0.0742  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0797  0.0049  0.1441  0.0049  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1566  0.0000  0.1566  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1442  0.0007  0.0346  0.0007  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 34">
       <set comment="band 1">
        <r>  0.1252  0.0155  0.0974  0.0155  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1311  0.0625  0.0022  0.0625  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0962  0.0107  0.1278  0.0107  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1545  0.0000  0.1545  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1455  0.0079  0.0207  0.0079  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 35">
       <set comment="band 1">
        <r>  0.1287  0.0599  0.0048  0.0599  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1228  0.0116  0.1117  0.0116  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1123  0.0200  0.1033  0.0200  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1504  0.0000  0.1504  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1355  0.0196  0.0097  0.0196  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 36">
       <set comment="band 1">
        <r>  0.1688  0.0467  0.0010  0.0467  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0629  0.0131  0.1451  0.0131  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1354  0.0360  0.0629  0.0360  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1449  0.0000  0.1449  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1173  0.0336  0.0026  0.0336  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 37">
       <set comment="band 1">
        <r>  0.1960  0.0356  0.0000  0.0356  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1967  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1663  0.0661  0.0000  0.0661  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1391  0.0000  0.1391  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0924  0.0481  0.0000  0.0481  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 38">
       <set comment="band 1">
        <r>  0.2025  0.0000  0.0657  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0576  0.0000  0.1496  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0293  0.0000  0.2896  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.2896  0.0000  0.0293  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 39">
       <set comment="band 1">
        <r>  0.1824  0.0009  0.0820  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0387  0.0191  0.1091  0.0191  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0627  0.0676  0.0374  0.0676  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1630  0.0000  0.1630  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0658  0.0361  0.0137  0.0361  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 40">
       <set comment="band 1">
        <r>  0.1551  0.0035  0.1001  0.0035  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0328  0.0859  0.0004  0.0859  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1295  0.0000  0.1343  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1640  0.0000  0.1640  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1373  0.0041  0.0158  0.0041  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 41">
       <set comment="band 1">
        <r>  0.1167  0.0076  0.1221  0.0076  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0817  0.0778  0.0000  0.0778  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1450  0.0037  0.1116  0.0037  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1623  0.0000  0.1623  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1586  0.0005  0.0084  0.0005  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 42">
       <set comment="band 1">
        <r>  0.0719  0.0758  0.0079  0.0758  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1113  0.0008  0.1496  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1773  0.0132  0.0701  0.0132  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1583  0.0000  0.1583  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1566  0.0073  0.0022  0.0073  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 43">
       <set comment="band 1">
        <r>  0.1386  0.0590  0.0000  0.0590  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2187  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2297  0.0329  0.0000  0.0329  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1523  0.0000  0.1523  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1418  0.0192  0.0000  0.0191  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 44">
       <set comment="band 1">
        <r>  0.1626  0.0000  0.1014  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0992  0.0000  0.1327  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0846  0.0000  0.2500  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.2500  0.0000  0.0846  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 45">
       <set comment="band 1">
        <r>  0.1331  0.0008  0.1249  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1223  0.0020  0.1131  0.0020  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0327  0.0861  0.0017  0.0861  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1690  0.0000  0.1689  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0715  0.0388  0.0033  0.0388  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 46">
       <set comment="band 1">
        <r>  0.0830  0.0024  0.1628  0.0024  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1438  0.0103  0.0746  0.0103  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0958  0.0780  0.0021  0.0780  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1679  0.0000  0.1679  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1448  0.0045  0.0017  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 47">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2348  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0029  0.0843  0.0000  0.0843  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.3440  0.0055  0.0000  0.0055  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1643  0.0000  0.1643  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1631  0.0004  0.0000  0.0004  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 48">
       <set comment="band 1">
        <r>  0.0894  0.0000  0.1647  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1738  0.0000  0.0785  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.2349  0.0000  0.1077  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.1077  0.0000  0.2349  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 49">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2432  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.2575  0.0024  0.0000  0.0024  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0320  0.0869  0.0000  0.0869  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1709  0.0000  0.1709  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0738  0.0397  0.0000  0.0397  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 50">
       <set comment="band 1">
        <r>  0.2490  0.0042  0.0171  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1571  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0528  0.1698  0.0608  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0878  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0296  0.0262  0.2683  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 51">
       <set comment="band 1">
        <r>  0.2334  0.0093  0.0262  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0337  0.0388  0.0189  0.1136  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0542  0.1073  0.0406  0.0602  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0267  0.0139  0.1992  0.0014  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0184  0.0322  0.0374  0.0569  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 52">
       <set comment="band 1">
        <r>  0.2121  0.0160  0.0366  0.0040  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0693  0.0462  0.0230  0.0930  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0512  0.0794  0.0262  0.0959  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0345  0.0152  0.1871  0.0050  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0321  0.0737  0.0564  0.0325  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 53">
       <set comment="band 1">
        <r>  0.1850  0.0241  0.0481  0.0086  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1048  0.0430  0.0206  0.0798  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0461  0.0640  0.0144  0.1179  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0456  0.0187  0.1806  0.0116  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0349  0.1079  0.0536  0.0201  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 54">
       <set comment="band 1">
        <r>  0.1451  0.0347  0.0618  0.0162  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1486  0.0342  0.0140  0.0659  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0283  0.0610  0.0010  0.1183  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0667  0.0164  0.1746  0.0357  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0266  0.1546  0.0466  0.0155  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 55">
       <set comment="band 1">
        <r>  0.2280  0.0041  0.0377  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1668  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0654  0.1983  0.0277  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0441  0.0005  0.1822  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1241  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 56">
       <set comment="band 1">
        <r>  0.2086  0.0089  0.0501  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0388  0.0423  0.0088  0.1218  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0596  0.1244  0.0130  0.0604  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0667  0.0012  0.1836  0.0027  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0337  0.0268  0.0236  0.0828  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 57">
       <set comment="band 1">
        <r>  0.1833  0.0152  0.0636  0.0038  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0785  0.0500  0.0100  0.0987  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0319  0.0953  0.0001  0.0843  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1047  0.0002  0.1783  0.0198  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0558  0.0646  0.0282  0.0479  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 58">
       <set comment="band 1">
        <r>  0.1495  0.0231  0.0784  0.0086  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1188  0.0450  0.0069  0.0831  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0355  0.0421  0.1112  0.0011  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1096  0.0400  0.0458  0.1333  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0616  0.0961  0.0222  0.0282  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 59">
       <set comment="band 1">
        <r>  0.0915  0.0518  0.0170  0.0803  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1774  0.0153  0.0799  0.0026  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0727  0.0372  0.1194  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0775  0.0401  0.0129  0.1568  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0584  0.1253  0.0149  0.0161  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 60">
       <set comment="band 1">
        <r>  0.1631  0.0337  0.0038  0.0616  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0996  0.0288  0.1074  0.0119  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0895  0.0490  0.0980  0.0138  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0631  0.0314  0.0047  0.1725  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0472  0.1589  0.0093  0.0097  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 61">
       <set comment="band 1">
        <r>  0.1929  0.0243  0.0010  0.0487  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0518  0.0250  0.1362  0.0114  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1069  0.0727  0.0596  0.0280  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0507  0.0235  0.0011  0.1844  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0224  0.2119  0.0101  0.0114  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 62">
       <set comment="band 1">
        <r>  0.2160  0.0161  0.0000  0.0369  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1831  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1258  0.1197  0.0000  0.0520  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0386  0.0175  0.0000  0.1926  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0874  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 63">
       <set comment="band 1">
        <r>  0.1987  0.0039  0.0653  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1732  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0086  0.1609  0.0136  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1466  0.0538  0.1641  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1436  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 64">
       <set comment="band 1">
        <r>  0.1748  0.0083  0.0810  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0384  0.0461  0.0016  0.1268  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0400  0.0285  0.1135  0.0056  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1418  0.1033  0.0499  0.0591  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0407  0.0254  0.0099  0.0963  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 65">
       <set comment="band 1">
        <r>  0.1432  0.0141  0.0984  0.0036  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0818  0.0531  0.0016  0.1024  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0903  0.0177  0.1257  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1045  0.0813  0.0146  0.1065  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0674  0.0630  0.0098  0.0555  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 66">
       <set comment="band 1">
        <r>  0.0727  0.0286  0.1187  0.0158  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1524  0.0383  0.0008  0.0774  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1179  0.0249  0.1058  0.0029  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0891  0.0613  0.0051  0.1352  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0733  0.0941  0.0053  0.0320  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 67">
       <set comment="band 1">
        <r>  0.1340  0.0443  0.0011  0.0761  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0729  0.0139  0.1507  0.0043  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1492  0.0415  0.0660  0.0109  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0769  0.0458  0.0011  0.1568  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0688  0.1216  0.0016  0.0174  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 68">
       <set comment="band 1">
        <r>  0.1676  0.0338  0.0000  0.0618  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2072  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1912  0.0757  0.0000  0.0273  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0644  0.0335  0.0000  0.1732  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0565  0.1490  0.0000  0.0088  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 69">
       <set comment="band 1">
        <r>  0.1589  0.0035  0.1007  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0843  0.0097  0.1268  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.1772  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1240  0.2150  0.0149  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1559  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 70">
       <set comment="band 1">
        <r>  0.1258  0.0072  0.1235  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0507  0.0343  0.0865  0.0236  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1158  0.0216  0.0249  0.1075  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1166  0.1287  0.0059  0.0644  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0445  0.0255  0.0024  0.1036  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 71">
       <set comment="band 1">
        <r>  0.0707  0.0105  0.1595  0.0032  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0485  0.0622  0.0037  0.1045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2110  0.0086  0.0700  0.0005  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1008  0.0877  0.0013  0.1074  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0724  0.0632  0.0011  0.0589  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 72">
       <set comment="band 1">
        <r>  0.0934  0.0558  0.0000  0.0920  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2265  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2561  0.0347  0.0000  0.0054  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0889  0.0636  0.0000  0.1355  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0766  0.0940  0.0000  0.0333  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 73">
       <set comment="band 1">
        <r>  0.0863  0.0024  0.1638  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1649  0.0072  0.0770  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.1791  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1235  0.2233  0.0014  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1621  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 74">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2390  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.2280  0.0228  0.0000  0.0035  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0697  0.0393  0.0000  0.1293  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1161  0.1314  0.0000  0.0649  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0457  0.0256  0.0000  0.1059  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 75">
       <set comment="band 1">
        <r>  0.1871  0.0153  0.0640  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1923  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1334  0.1345  0.0159  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0930  0.0047  0.1774  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1172  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 76">
       <set comment="band 1">
        <r>  0.1593  0.0228  0.0789  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0572  0.0273  0.0037  0.1418  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0997  0.0839  0.0038  0.0604  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1205  0.0089  0.1623  0.0031  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0201  0.0635  0.0109  0.0769  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 77">
       <set comment="band 1">
        <r>  0.1179  0.0326  0.0960  0.0048  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1144  0.0287  0.0018  0.1097  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0403  0.0692  0.0029  0.0807  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1755  0.0070  0.1354  0.0263  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0184  0.1470  0.0161  0.0493  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 78">
       <set comment="band 1">
        <r>  0.1476  0.0140  0.0989  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1986  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0041  0.1122  0.0362  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2837  0.0330  0.1083  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1298  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 79">
       <set comment="band 1">
        <r>  0.1098  0.0201  0.1205  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0552  0.0303  0.0000  0.1444  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0483  0.0511  0.0859  0.0081  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2428  0.0433  0.0286  0.0579  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0260  0.0548  0.0025  0.0854  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 80">
       <set comment="band 1">
        <r>  0.0628  0.0450  0.0072  0.1146  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0991  0.0086  0.1477  0.0013  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1299  0.0506  0.0662  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1711  0.0334  0.0037  0.1093  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0352  0.1109  0.0013  0.0524  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 81">
       <set comment="band 1">
        <r>  0.1278  0.0311  0.0000  0.0958  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2145  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1882  0.0855  0.0000  0.0029  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1405  0.0202  0.0000  0.1394  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0267  0.1640  0.0000  0.0357  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 82">
       <set comment="band 1">
        <r>  0.0771  0.0094  0.1608  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1266  0.0329  0.0720  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.2018  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2430  0.1194  0.0038  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1352  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 83">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2305  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0005  0.0686  0.0000  0.1096  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2621  0.0127  0.0000  0.0390  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2000  0.0643  0.0000  0.0646  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0275  0.0534  0.0000  0.0879  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 84">
       <set comment="band 1">
        <r>  0.0590  0.0220  0.1549  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.2181  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0093  0.1139  0.0296  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.4033  0.0009  0.0423  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0021  0.0004  0.1100  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
     </set>
     <set comment="spin2">
      <set comment="kpoint 1">
       <set comment="band 1">
        <r>  0.2701  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0460  0.2380  0.2323  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 2">
       <set comment="band 1">
        <r>  0.2669  0.0011  0.0011  0.0011  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0073  0.0501  0.0501  0.0501  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0045  0.0779  0.0451  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0804  0.0071  0.0401  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0244  0.0552  0.1027  0.1278  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 3">
       <set comment="band 1">
        <r>  0.2574  0.0043  0.0043  0.0043  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0266  0.0565  0.0565  0.0565  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1193  0.0098  0.0730  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0154  0.1249  0.0617  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0677  0.0526  0.0542  0.0512  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 4">
       <set comment="band 1">
        <r>  0.2418  0.0094  0.0094  0.0094  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0521  0.0559  0.0559  0.0559  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1443  0.0666  0.0236  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0120  0.0897  0.1327  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1035  0.0286  0.0244  0.0231  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 5">
       <set comment="band 1">
        <r>  0.2205  0.0162  0.0162  0.0162  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0818  0.0520  0.0520  0.0520  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0342  0.1176  0.0989  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1317  0.0491  0.0676  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0250  0.0948  0.0810  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 6">
       <set comment="band 1">
        <r>  0.1935  0.0244  0.0244  0.0244  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1147  0.0456  0.0456  0.0456  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1619  0.0251  0.0691  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0086  0.1416  0.0988  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.1348  0.0212  0.0570  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 7">
       <set comment="band 1">
        <r>  0.1543  0.0352  0.0352  0.0352  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1570  0.0357  0.0357  0.0357  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0075  0.1463  0.1091  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1431  0.0246  0.0564  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0066  0.1193  0.0901  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 8">
       <set comment="band 1">
        <r>  0.2659  0.0000  0.0044  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0918  0.0000  0.1081  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1081  0.0000  0.0918  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0449  0.0000  0.4510  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 9">
       <set comment="band 1">
        <r>  0.2584  0.0011  0.0097  0.0011  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0141  0.0712  0.0227  0.0712  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1182  0.0000  0.1182  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0214  0.0086  0.0750  0.0086  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0430  0.0466  0.3280  0.0467  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 10">
       <set comment="band 1">
        <r>  0.2448  0.0042  0.0170  0.0042  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0388  0.0691  0.0322  0.0691  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1288  0.0000  0.1288  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0460  0.0001  0.1196  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0129  0.0127  0.2024  0.0127  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 11">
       <set comment="band 1">
        <r>  0.2254  0.0092  0.0258  0.0092  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0680  0.0646  0.0332  0.0646  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1345  0.0000  0.1345  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0544  0.0045  0.1049  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0204  0.0132  0.1992  0.0133  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 12">
       <set comment="band 1">
        <r>  0.2004  0.0157  0.0359  0.0157  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1002  0.0579  0.0300  0.0579  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1365  0.0000  0.1365  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0065  0.0082  0.1029  0.0082  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0785  0.0224  0.1821  0.0224  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 13">
       <set comment="band 1">
        <r>  0.1676  0.0239  0.0473  0.0239  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1373  0.0488  0.0238  0.0488  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1355  0.0000  0.1355  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0234  0.0187  0.1540  0.0187  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0657  0.0284  0.1063  0.0284  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 14">
       <set comment="band 1">
        <r>  0.1098  0.0530  0.0314  0.0530  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1961  0.0181  0.0440  0.0181  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0299  0.0264  0.1504  0.0264  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1327  0.0000  0.1327  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0573  0.0393  0.0813  0.0391  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 15">
       <set comment="band 1">
        <r>  0.1799  0.0348  0.0153  0.0348  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1237  0.0336  0.0638  0.0336  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0339  0.0353  0.1383  0.0353  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1304  0.0000  0.1302  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0456  0.0548  0.0627  0.0451  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 16">
       <set comment="band 1">
        <r>  0.2094  0.0251  0.0089  0.0251  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0885  0.0393  0.0740  0.0393  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0357  0.0457  0.1205  0.0457  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1572  0.0000  0.1570  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0661  0.0000  0.0654  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 17">
       <set comment="band 1">
        <r>  0.2326  0.0166  0.0041  0.0166  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0565  0.0410  0.0845  0.0410  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0355  0.0581  0.0967  0.0581  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1419  0.0000  0.1419  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0663  0.0000  0.0668  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 18">
       <set comment="band 1">
        <r>  0.2501  0.0096  0.0011  0.0096  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0262  0.0341  0.1010  0.0341  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0342  0.0761  0.0619  0.0761  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0791  0.0000  0.0790  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0005  0.1508  0.0035  0.1466  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 19">
       <set comment="band 1">
        <r>  0.2616  0.0043  0.0000  0.0043  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1429  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0300  0.1127  0.0000  0.1127  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.2712  0.0001  0.2737  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 20">
       <set comment="band 1">
        <r>  0.2531  0.0000  0.0172  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.1914  0.0000  0.0697  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0697  0.0000  0.1914  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1292  0.0000  0.3102  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 21">
       <set comment="band 1">
        <r>  0.2416  0.0010  0.0266  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0180  0.0793  0.0116  0.0793  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1393  0.0000  0.1393  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0325  0.0238  0.0293  0.0238  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0301  0.0003  0.1869  0.0003  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 22">
       <set comment="band 1">
        <r>  0.2241  0.0041  0.0374  0.0041  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0472  0.0763  0.0172  0.0763  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1448  0.0000  0.1448  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0017  0.0045  0.0647  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1218  0.0001  0.1835  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 23">
       <set comment="band 1">
        <r>  0.2008  0.0088  0.0494  0.0088  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0797  0.0702  0.0170  0.0702  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1472  0.0000  0.1472  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0349  0.0058  0.1479  0.0058  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1204  0.0019  0.0884  0.0019  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 24">
       <set comment="band 1">
        <r>  0.1710  0.0151  0.0625  0.0151  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1153  0.0616  0.0135  0.0616  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0506  0.0104  0.1490  0.0104  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1470  0.0000  0.1470  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1207  0.0099  0.0614  0.0099  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 25">
       <set comment="band 1">
        <r>  0.0996  0.0314  0.0815  0.0314  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1887  0.0421  0.0032  0.0421  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0622  0.0174  0.1381  0.0174  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1446  0.0000  0.1446  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1139  0.0215  0.0411  0.0215  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 26">
       <set comment="band 1">
        <r>  0.1592  0.0466  0.0087  0.0466  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1267  0.0223  0.0851  0.0223  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0711  0.0269  0.1203  0.0269  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1409  0.0000  0.1409  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1007  0.0349  0.0251  0.0349  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 27">
       <set comment="band 1">
        <r>  0.1922  0.0354  0.0038  0.0354  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0859  0.0266  0.1018  0.0266  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0785  0.0400  0.0956  0.0400  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1372  0.0000  0.1372  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0824  0.0485  0.0131  0.0486  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 28">
       <set comment="band 1">
        <r>  0.2172  0.0254  0.0010  0.0254  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0438  0.0224  0.1274  0.0224  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0879  0.0609  0.0578  0.0609  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1363  0.0000  0.1364  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0609  0.0583  0.0046  0.0633  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 29">
       <set comment="band 1">
        <r>  0.2366  0.0167  0.0000  0.0167  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1712  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0975  0.0991  0.0000  0.0991  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1490  0.0000  0.1488  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0001  0.0000  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 30">
       <set comment="band 1">
        <r>  0.2320  0.0000  0.0379  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.1214  0.0000  0.1738  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1738  0.0000  0.1214  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0262  0.0000  0.1465  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 31">
       <set comment="band 1">
        <r>  0.2164  0.0010  0.0508  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0190  0.0840  0.0040  0.0840  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.1533  0.0000  0.1533  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0384  0.0008  0.1466  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0617  0.0304  0.0400  0.0304  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 32">
       <set comment="band 1">
        <r>  0.1949  0.0038  0.0648  0.0038  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0512  0.0814  0.0062  0.0814  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0609  0.0017  0.1521  0.0017  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1563  0.0000  0.1562  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1233  0.0030  0.0468  0.0030  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 33">
       <set comment="band 1">
        <r>  0.1669  0.0083  0.0800  0.0083  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0869  0.0742  0.0054  0.0742  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0797  0.0049  0.1441  0.0049  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1566  0.0000  0.1566  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1442  0.0007  0.0346  0.0007  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 34">
       <set comment="band 1">
        <r>  0.1252  0.0155  0.0974  0.0155  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1311  0.0625  0.0022  0.0625  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0962  0.0107  0.1278  0.0107  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1545  0.0000  0.1545  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1455  0.0079  0.0207  0.0079  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 35">
       <set comment="band 1">
        <r>  0.1287  0.0599  0.0048  0.0599  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1228  0.0116  0.1117  0.0116  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1123  0.0200  0.1033  0.0200  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1504  0.0000  0.1504  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1355  0.0196  0.0097  0.0196  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 36">
       <set comment="band 1">
        <r>  0.1688  0.0467  0.0010  0.0467  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0629  0.0131  0.1451  0.0131  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1354  0.0360  0.0629  0.0360  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1449  0.0000  0.1449  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1173  0.0336  0.0026  0.0336  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 37">
       <set comment="band 1">
        <r>  0.1960  0.0356  0.0000  0.0356  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1967  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1663  0.0661  0.0000  0.0661  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1391  0.0000  0.1391  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0924  0.0481  0.0000  0.0481  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 38">
       <set comment="band 1">
        <r>  0.2025  0.0000  0.0657  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0576  0.0000  0.1496  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0055  0.0000  0.3134  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.3134  0.0000  0.0055  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 39">
       <set comment="band 1">
        <r>  0.1824  0.0009  0.0820  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0387  0.0191  0.1091  0.0191  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0627  0.0676  0.0374  0.0676  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1630  0.0000  0.1630  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0658  0.0361  0.0137  0.0360  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 40">
       <set comment="band 1">
        <r>  0.1551  0.0035  0.1001  0.0035  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0328  0.0859  0.0004  0.0859  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1295  0.0000  0.1343  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1640  0.0000  0.1640  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1373  0.0041  0.0158  0.0041  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 41">
       <set comment="band 1">
        <r>  0.1167  0.0076  0.1221  0.0076  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0817  0.0778  0.0000  0.0778  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1450  0.0037  0.1116  0.0037  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1623  0.0000  0.1623  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1586  0.0005  0.0084  0.0005  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 42">
       <set comment="band 1">
        <r>  0.0719  0.0758  0.0079  0.0758  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1113  0.0008  0.1496  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1773  0.0132  0.0701  0.0132  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1583  0.0000  0.1583  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1566  0.0073  0.0022  0.0073  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 43">
       <set comment="band 1">
        <r>  0.1386  0.0590  0.0000  0.0590  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2187  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2297  0.0329  0.0000  0.0329  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1523  0.0000  0.1523  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1418  0.0192  0.0000  0.0192  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 44">
       <set comment="band 1">
        <r>  0.1626  0.0000  0.1014  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0992  0.0000  0.1327  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0343  0.0000  0.3004  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.3004  0.0000  0.0343  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 45">
       <set comment="band 1">
        <r>  0.1331  0.0008  0.1249  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1223  0.0020  0.1131  0.0020  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0327  0.0861  0.0017  0.0861  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1689  0.0000  0.1689  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0716  0.0388  0.0033  0.0388  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 46">
       <set comment="band 1">
        <r>  0.0830  0.0024  0.1628  0.0024  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1438  0.0103  0.0746  0.0103  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0958  0.0780  0.0021  0.0780  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1679  0.0000  0.1679  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1448  0.0045  0.0017  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 47">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2348  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0029  0.0843  0.0000  0.0843  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.3440  0.0055  0.0000  0.0055  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1643  0.0000  0.1643  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.1631  0.0004  0.0000  0.0004  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 48">
       <set comment="band 1">
        <r>  0.0894  0.0000  0.1647  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1738  0.0000  0.0785  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0288  0.0000  0.3138  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.3138  0.0000  0.0288  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 49">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2432  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.2575  0.0024  0.0000  0.0024  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0320  0.0869  0.0000  0.0869  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.1709  0.0000  0.1709  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0738  0.0397  0.0000  0.0397  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 50">
       <set comment="band 1">
        <r>  0.2490  0.0042  0.0171  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1571  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0528  0.1698  0.0608  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0000  0.0000  0.0000  0.0878  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0307  0.0280  0.2741  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 51">
       <set comment="band 1">
        <r>  0.2334  0.0093  0.0262  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0337  0.0388  0.0189  0.1136  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0542  0.1073  0.0406  0.0602  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0267  0.0139  0.1991  0.0014  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0184  0.0323  0.0375  0.0569  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 52">
       <set comment="band 1">
        <r>  0.2121  0.0160  0.0366  0.0040  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0693  0.0462  0.0230  0.0930  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0512  0.0794  0.0262  0.0959  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0345  0.0152  0.1871  0.0050  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0321  0.0737  0.0564  0.0326  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 53">
       <set comment="band 1">
        <r>  0.1850  0.0241  0.0481  0.0086  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1048  0.0430  0.0206  0.0798  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0461  0.0640  0.0144  0.1179  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0456  0.0187  0.1806  0.0116  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0349  0.1078  0.0536  0.0201  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 54">
       <set comment="band 1">
        <r>  0.1451  0.0347  0.0618  0.0162  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1486  0.0342  0.0140  0.0659  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0283  0.0610  0.0010  0.1183  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0667  0.0164  0.1746  0.0357  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0222  0.1681  0.0469  0.0179  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 55">
       <set comment="band 1">
        <r>  0.2280  0.0041  0.0377  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1668  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0654  0.1983  0.0277  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0441  0.0005  0.1822  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1241  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 56">
       <set comment="band 1">
        <r>  0.2086  0.0089  0.0501  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0388  0.0423  0.0088  0.1218  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0596  0.1244  0.0130  0.0604  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0667  0.0012  0.1836  0.0027  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0336  0.0267  0.0235  0.0829  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 57">
       <set comment="band 1">
        <r>  0.1833  0.0152  0.0636  0.0038  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0785  0.0500  0.0100  0.0987  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0319  0.0953  0.0001  0.0843  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1047  0.0002  0.1783  0.0198  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0558  0.0647  0.0282  0.0479  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 58">
       <set comment="band 1">
        <r>  0.1495  0.0231  0.0784  0.0086  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1188  0.0450  0.0069  0.0831  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0355  0.0421  0.1112  0.0011  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1096  0.0400  0.0458  0.1333  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0616  0.0961  0.0222  0.0282  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 59">
       <set comment="band 1">
        <r>  0.0915  0.0518  0.0170  0.0803  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1774  0.0153  0.0799  0.0026  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0727  0.0372  0.1194  0.0045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0775  0.0401  0.0129  0.1568  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0585  0.1252  0.0149  0.0161  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 60">
       <set comment="band 1">
        <r>  0.1631  0.0337  0.0038  0.0616  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0996  0.0288  0.1074  0.0119  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0895  0.0490  0.0980  0.0138  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0631  0.0314  0.0047  0.1725  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0473  0.1586  0.0093  0.0096  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 61">
       <set comment="band 1">
        <r>  0.1929  0.0243  0.0010  0.0487  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0518  0.0250  0.1362  0.0114  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1069  0.0727  0.0596  0.0280  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0507  0.0235  0.0011  0.1845  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0103  0.2111  0.0236  0.0175  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 62">
       <set comment="band 1">
        <r>  0.2160  0.0161  0.0000  0.0369  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.1831  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1258  0.1197  0.0000  0.0520  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0386  0.0175  0.0000  0.1926  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0874  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 63">
       <set comment="band 1">
        <r>  0.1987  0.0039  0.0653  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1732  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0087  0.1609  0.0136  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1466  0.0538  0.1641  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1436  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 64">
       <set comment="band 1">
        <r>  0.1748  0.0083  0.0810  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0384  0.0461  0.0016  0.1268  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0400  0.0285  0.1135  0.0056  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1418  0.1033  0.0499  0.0591  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0407  0.0254  0.0099  0.0963  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 65">
       <set comment="band 1">
        <r>  0.1432  0.0141  0.0984  0.0036  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0818  0.0531  0.0016  0.1024  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0903  0.0177  0.1257  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1045  0.0813  0.0146  0.1065  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0674  0.0631  0.0098  0.0555  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 66">
       <set comment="band 1">
        <r>  0.0727  0.0286  0.1187  0.0158  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1524  0.0383  0.0008  0.0774  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1179  0.0249  0.1058  0.0029  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0891  0.0613  0.0051  0.1352  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0733  0.0941  0.0053  0.0320  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 67">
       <set comment="band 1">
        <r>  0.1340  0.0443  0.0011  0.0761  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0729  0.0139  0.1507  0.0043  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1492  0.0415  0.0660  0.0109  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0769  0.0458  0.0011  0.1568  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0688  0.1216  0.0016  0.0174  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 68">
       <set comment="band 1">
        <r>  0.1676  0.0338  0.0000  0.0618  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2072  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1912  0.0757  0.0000  0.0273  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0644  0.0335  0.0000  0.1732  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0565  0.1490  0.0000  0.0088  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 69">
       <set comment="band 1">
        <r>  0.1589  0.0035  0.1007  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0843  0.0097  0.1268  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.1772  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1240  0.2150  0.0149  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1559  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 70">
       <set comment="band 1">
        <r>  0.1258  0.0072  0.1235  0.0008  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0507  0.0343  0.0865  0.0236  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1158  0.0216  0.0249  0.1075  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1166  0.1287  0.0059  0.0644  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0445  0.0255  0.0024  0.1036  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 71">
       <set comment="band 1">
        <r>  0.0707  0.0105  0.1595  0.0032  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0485  0.0622  0.0037  0.1045  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2110  0.0086  0.0700  0.0005  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1008  0.0877  0.0013  0.1074  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0724  0.0632  0.0011  0.0589  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 72">
       <set comment="band 1">
        <r>  0.0934  0.0558  0.0000  0.0920  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2265  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2561  0.0347  0.0000  0.0054  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0889  0.0636  0.0000  0.1355  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0766  0.0940  0.0000  0.0333  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 73">
       <set comment="band 1">
        <r>  0.0863  0.0024  0.1638  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1649  0.0072  0.0770  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.1791  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1235  0.2233  0.0014  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1621  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 74">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2390  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.2280  0.0228  0.0000  0.0035  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0697  0.0393  0.0000  0.1293  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1161  0.1314  0.0000  0.0649  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0456  0.0256  0.0000  0.1060  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 75">
       <set comment="band 1">
        <r>  0.1871  0.0153  0.0640  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1923  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1334  0.1345  0.0159  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.0930  0.0047  0.1774  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1172  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 76">
       <set comment="band 1">
        <r>  0.1593  0.0228  0.0789  0.0009  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0572  0.0273  0.0037  0.1418  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0997  0.0839  0.0038  0.0604  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1205  0.0089  0.1623  0.0031  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0201  0.0632  0.0109  0.0770  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 77">
       <set comment="band 1">
        <r>  0.1179  0.0326  0.0960  0.0048  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1144  0.0287  0.0018  0.1097  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0403  0.0692  0.0029  0.0807  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1755  0.0070  0.1354  0.0263  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0175  0.1498  0.0166  0.0492  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 78">
       <set comment="band 1">
        <r>  0.1476  0.0140  0.0989  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.1986  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0041  0.1122  0.0362  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2837  0.0330  0.1083  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1298  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 79">
       <set comment="band 1">
        <r>  0.1098  0.0201  0.1205  0.0010  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0552  0.0303  0.0000  0.1444  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0483  0.0511  0.0859  0.0081  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2428  0.0433  0.0286  0.0579  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0262  0.0553  0.0025  0.0852  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 80">
       <set comment="band 1">
        <r>  0.0628  0.0450  0.0072  0.1146  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0991  0.0086  0.1477  0.0013  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1299  0.0506  0.0662  0.0001  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1711  0.0334  0.0037  0.1093  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0352  0.1109  0.0013  0.0524  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 81">
       <set comment="band 1">
        <r>  0.1278  0.0311  0.0000  0.0958  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.2145  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.1882  0.0855  0.0000  0.0029  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.1405  0.0202  0.0000  0.1394  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0266  0.1647  0.0000  0.0357  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 82">
       <set comment="band 1">
        <r>  0.0771  0.0094  0.1608  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.1266  0.0329  0.0720  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0000  0.0000  0.0000  0.2018  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2430  0.1194  0.0038  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0000  0.0000  0.1352  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 83">
       <set comment="band 1">
        <r>  0.0000  0.0000  0.2305  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0005  0.0686  0.0000  0.1096  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.2621  0.0127  0.0000  0.0390  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.2000  0.0643  0.0000  0.0646  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0276  0.0535  0.0000  0.0878  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
      <set comment="kpoint 84">
       <set comment="band 1">
        <r>  0.0590  0.0220  0.1549  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 2">
        <r>  0.0000  0.0000  0.0000  0.2181  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 3">
        <r>  0.0093  0.1139  0.0296  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 4">
        <r>  0.4033  0.0009  0.0423  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
       <set comment="band 5">
        <r>  0.0000  0.0002  0.0000  0.1109  0.0000  0.0000  0.0000  0.0000  0.0000 </r>
       </set>
      </set>
     </set>
    </set>
   </array>
  </projected>
 </calculation>
 <structure name="finalpos" >
  <crystal>
   <varray name="basis" >
    <v>       0.00000000       2.02055508       2.02055508 </v>
    <v>       2.02055508       0.00000000       2.02055508 </v>
    <v>       2.02055508       2.02055508      -0.00000000 </v>
   </varray>
   <i name="volume">     16.49840937 </i>
   <varray name="rec_basis" >
    <v>      -0.24745675       0.24745675       0.24745675 </v>
    <v>       0.24745675      -0.24745675       0.24745675 </v>
    <v>       0.24745675       0.24745675      -0.24745675 </v>
   </varray>
  </crystal>
  <varray name="positions" >
   <v>       0.00000000       0.00000000       0.00000000 </v>
  </varray>
 </structure>
</modeling>
