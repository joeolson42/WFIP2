# WFIP2
development version of the RAP/HRRR physics in WRF-ARW

Note that this repository houses the WRF-ARW code under development for WFIP2, which is in the
WRFV3.7.1 directory. The code base is actually much closer to WRFV3.8, with new modifications from NOAA/GSD.
The following directories are exactly V3.8:

WRFV3.7.1/frame
WRFV3.7.1/external
WRFV3.7.1/arch

Many of the other directories are close to V3.8.

The most important HRRR-related physics for WFIP2 are:
phys/module_bl_mynn.F
phys/module_sf_mynn.F
phys/module_sf_ruclsm.F

All of these modules have been updated since v3.8, with some stochastic physics options.
