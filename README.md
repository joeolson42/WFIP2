# WFIP2
This is the development version of the RAP/HRRR physics in WRF-ARW in support of WFIP2

Note that this repository houses the WRF-ARW code for two difference code bases: WRFV3.7.1 and WRFV3.9. The v3.7.1 code base is actually much closer to WRFV3.8.1, with new modifications from NOAA/GSD. The v3.9 code also has many bug fixes included in the official release of WRF-ARWv3.9.1, plus many new code developments from NOAA/GSD.

The most important physics development for WFIP2 can be found in:

phys/module_bl_mynn.F
phys/module_sf_mynn.F
phys/module_bl_gwdo.F

As of 22 August 2017, this code has been frozen for final WFIP2 testing.
