@ECHO OFF
@ECHO.


REM  Set up environment variables.  You will probably have to change these.


@SET Compare=FC /T
@SET CRUNCH=..\bin\crunch_win32.exe
rem @SET CRUNCH=Call Crunch

@SET MATLAB=matlab
@SET MBC_SOURCE=C:\Users\bjonkman\Documents\DATA\Downloads\MBC\Source
rem @SET MBC_SOURCE=C:\Users\bjonkman\Data\DesignCodes\MBC\Source
@SET DateTime=DateTime.exe
@SET Editor=NotePad.EXE
@SET CompareFile=Samples.out

::=======================================================================================================
IF /I "%1"=="-DEBUG" GOTO debugVer
IF /I "%1"=="-GFORTRAN" GOTO gfortran
IF /I "%1"=="-IFORT" GOTO ifort
IF /I "%1"=="-DEVBUILD" GOTO devBuild
IF /I "%1"=="-DEVDEBUG" GOTO devDebugBuild

:releaseVer
@SET EXE_VER=Using released version of HydroDynDriver (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_win32.exe
goto Samples

:debugVer
@SET EXE_VER=Using HydroDynDriver compiled in debug mode (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_debug_win32.exe
goto Samples

:gfortran
@SET EXE_VER=Using HydroDynDriver compiled with makefile (gfortran)
@SET HydroDyn=..\compiling\HydroDynDriver_gwin32.exe
goto Samples

:ifort
@SET EXE_VER=Using HydroDynDriver compiled with Compile_HydroDynDriver.bat (IVF)
@SET HydroDyn=..\compiling\HydroDynDriver_iwin32.exe
goto Samples

:devBuild
@SET EXE_VER=Using HydroDynDriver compiled with Visual Studio Project, release mode (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_dev_win32.exe
goto Samples

:devDebugBuild
@SET EXE_VER=Using HydroDynDriver compiled with Visual Studio Project, debug mode (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_dev_debug_win32.exe
goto Samples

::=======================================================================================================


:Samples


REM  HydroDyn test sequence definition:

@SET  TEST01=Test #01: Version of the OC3 Monopile.
@SET  TEST02=Test #02: OC3 Tripod.
@SET  TEST03=Test #03: OC4 Jacket under static displacement.
@SET  TEST04=Test #04: OC3 Hywind.
@SET  TEST05=Test #05: MIT/NREL TLP.
@SET  TEST06=Test #06: ITI Barge4.

@SET  DASHES=---------------------------------------------------------------------------------------------
@SET  POUNDS=#############################################################################################

@IF EXIST Samples.out  DEL Samples.out

ECHO.                                               >> Samples.out
ECHO           ************************************ >> Samples.out
ECHO           **  HydroDyn Stand-alone Samples  ** >> Samples.out
ECHO           ************************************ >> Samples.out

ECHO.                                                                             >> Samples.out
ECHO ############################################################################ >> Samples.out
ECHO # Inspect this file for any differences between your results and the saved # >> Samples.out
ECHO # results.  Any differing lines and the two lines surrounding them will be # >> Samples.out
ECHO # listed.  The only differences should be the time stamps at the start of  # >> Samples.out
ECHO # each file.                                                               # >> Samples.out
ECHO #                                                                          # >> Samples.out
ECHO # If you are running on something other than a PC, you may see differences # >> Samples.out
ECHO # in the last significant digit of many of the numbers.                    # >> Samples.out
ECHO ############################################################################ >> Samples.out

rem ECHO.                                            >> Samples.out
rem ECHO Date and time this acceptance test was run: >> Samples.out
rem %DateTime%                                       >> Samples.out
rem ECHO.                                            >> Samples.out


ECHO.                                            >> Samples.out
ECHO %EXE_VER%                                   >> Samples.out
ECHO HydroDyn called with this command:              >> Samples.out
ECHO %HydroDyn%                                      >> Samples.out
ECHO.                                            >> Samples.out


echo %DASHES%
echo %EXE_VER%
echo %HydroDyn%
echo %DASHES%


rem *******************************************************
:Test1
@CALL :GenTestHeader %Test01%
@CALL :RunHydroDyn NRELOffshrBsline5MW_OC3Monopile_HydroDyn.dvr

rem *******************************************************
:Test2
@CALL :GenTestHeader %Test02%
@CALL :RunHydroDyn NRELOffshrBsline5MW_OC3Tripod_HydroDyn.dvr

rem *******************************************************
:Test3
@CALL :GenTestHeader %Test03%
@CALL :RunHydroDyn NRELOffshrBsline5MW_OC4Jacket_HydroDyn.dvr

rem *******************************************************
:Test4
@CALL :GenTestHeader %Test04%
@CALL :RunHydroDyn NRELOffshrBsline5MW_OC3Hywind_HydroDyn.dvr

rem *******************************************************
:Test5
@CALL :GenTestHeader %Test05%
@CALL :RunHydroDyn NRELOffshrBsline5MW_MIT_NREL_TLP_HydroDyn.dvr

rem *******************************************************
:Test6
@CALL :GenTestHeader %Test06%
@CALL :RunHydroDyn NRELOffshrBsline5MW_ITIBarge4_HydroDyn.dvr

rem ******************************************************
rem  Let's look at the comparisons.
:MatlabComparisons

rem %MATLAB% /r PlotSamplesResults('.','.\TstFiles');exit;


rem %Editor% Samples.out
goto END

rem ******************************************************
:GenTestHeader
echo %POUNDS%
@echo HydroDyn %*
echo %POUNDS%

echo.                                    >> %CompareFile%
echo %POUNDS%                            >> %CompareFile%
echo.                                    >> %CompareFile%
echo %*                                  >> %CompareFile%
EXIT /B

rem ******************************************************
:RunHydroDyn
:: Run HydroDyn.
@SET TEST=%1

%HydroDyn% %1

IF ERRORLEVEL 1  GOTO ERROR

echo %DASHES%



EXIT /B




:ERROR
:: Sets clears memory and stops the batch immediately
@echo ** An error has occurred in Test #%TEST% **
@echo ** An error has occurred in Test #%TEST% ** >> %CompareFile%

@call :end
@call :__ErrorExit 2> nul
EXIT /B

:__ErrorExit
rem Creates a syntax error, stops immediately
()
EXIT /B


:END

@SET CRUNCH=
@SET MATLAB=
@SET MBC_SOURCE=
@SET Compare=
@SET CompareFile=
@SET DASHES=
@SET DateTime=
@SET Editor=
@SET HydroDyn=
@SET POUNDS=
@SET TEST=
@SET TEST01=
@SET TEST02=
@SET TEST03=
@SET TEST04=


SET EXE_VER=

rem type Bell.txt
@echo Processing complete.

EXIT /B