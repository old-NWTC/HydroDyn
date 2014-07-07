**********************************************************************************************
*** There is a known issue with the included HydroDynDriver.vfproj and Visual Studio 2010. ***
**********************************************************************************************

- When you first open the project in Visual Studio and build the solution, you may encounter build errors like this:
      error #7002: Error in opening the compiled module file.  Check INCLUDE paths.   [NWTC_LIBRARY]

There is a problem with the build order this first time you use the project.  At this time the reason is unkown.

If you receive these build errors, try shutting down Visual Studio, and then reopen the HydroDynDriver.vfproj project, and 
try building the solution again.  This usually takes care of the build errors.  

