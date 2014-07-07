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
@SET CompareFile=CertTest.out

::=======================================================================================================
IF /I "%1"=="-DEBUG" GOTO debugVer
IF /I "%1"=="-GFORTRAN" GOTO gfortran
IF /I "%1"=="-IFORT" GOTO ifort
IF /I "%1"=="-DEVBUILD" GOTO devBuild
IF /I "%1"=="-DEVDEBUG" GOTO devDebugBuild

:releaseVer
@SET EXE_VER=Using released version of HydroDynDriver (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_win32.exe
goto CertTest

:debugVer
@SET EXE_VER=Using HydroDynDriver compiled in debug mode (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_debug_win32.exe
goto CertTest

:gfortran
@SET EXE_VER=Using HydroDynDriver compiled with makefile (gfortran)
@SET HydroDyn=..\compiling\HydroDynDriver_gwin32.exe
goto CertTest

:ifort
@SET EXE_VER=Using HydroDynDriver compiled with Compile_HydroDynDriver.bat (IVF)
@SET HydroDyn=..\compiling\HydroDynDriver_iwin32.exe
goto CertTest

:devBuild
@SET EXE_VER=Using HydroDynDriver compiled with Visual Studio Project, release mode (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_dev_win32.exe
goto CertTest

:devDebugBuild
@SET EXE_VER=Using HydroDynDriver compiled with Visual Studio Project, debug mode (IVF/VS)
@SET HydroDyn=..\bin\HydroDynDriver_dev_debug_win32.exe
goto CertTest

::=======================================================================================================


:CertTest


REM  HydroDyn test sequence definition:

@SET  TEST01=Test #01: Simple tripod.
@SET  TEST02=Test #02: A 2D house frame aligned along the X-axis.


@SET  DASHES=---------------------------------------------------------------------------------------------
@SET  POUNDS=#############################################################################################

@IF EXIST CertTest.out  DEL CertTest.out

ECHO.                                               >> CertTest.out
ECHO           ************************************ >> CertTest.out
ECHO           **  HydroDyn Stand-alone CertTest ** >> CertTest.out
ECHO           ************************************ >> CertTest.out

ECHO.                                                                             >> CertTest.out
ECHO ############################################################################ >> CertTest.out
ECHO # Inspect this file for any differences between your results and the saved # >> CertTest.out
ECHO # results.  Any differing lines and the two lines surrounding them will be # >> CertTest.out
ECHO # listed.  The only differences should be the time stamps at the start of  # >> CertTest.out
ECHO # each file.                                                               # >> CertTest.out
ECHO #                                                                          # >> CertTest.out
ECHO # If you are running on something other than a PC, you may see differences # >> CertTest.out
ECHO # in the last significant digit of many of the numbers.                    # >> CertTest.out
ECHO ############################################################################ >> CertTest.out

rem ECHO.                                            >> CertTest.out
rem ECHO Date and time this acceptance test was run: >> CertTest.out
rem %DateTime%                                       >> CertTest.out
rem ECHO.                                            >> CertTest.out


ECHO.                                            >> CertTest.out
ECHO %EXE_VER%                                   >> CertTest.out
ECHO HydroDyn called with this command:              >> CertTest.out
ECHO %HydroDyn%                                      >> CertTest.out
ECHO.                                            >> CertTest.out


echo %DASHES%
echo %EXE_VER%
echo %HydroDyn%
echo %DASHES%


rem *******************************************************
:Test1
@CALL :GenTestHeader %Test01%
@CALL :RunHydroDyn Test_001.dvr

rem *******************************************************
:Test2
@CALL :GenTestHeader %Test02%
@CALL :RunHydroDyn Test_002.dvr


rem ******************************************************
rem  Let's look at the comparisons.
:MatlabComparisons

rem %MATLAB% /r PlotCertTestResults('.','.\TstFiles');exit;


rem %Editor% CertTest.out
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