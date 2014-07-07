@ECHO OFF

set lines=======================================================================
echo %lines%
IF "%1"=="" (
ECHO.
ECHO   The calling syntax for this script is
ECHO             RunRegistry ModuleName
ECHO.
GOTO Done
)


REM ----------------------------------------------------------------------------
REM ------------------------- LOCAL PATHS --------------------------------------
REM ----------------------------------------------------------------------------
REM -- USERS MAY EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL MACHINES. -
REM -- NOTE: do not use quotation marks around the path names!!!! --------------
REM ----------------------------------------------------------------------------
REM ----------------------------------------------------------------------------

SET Registry=..\..\bin\Registry_win32.exe
SET Source_Loc=..\..\Source


SET HD_Reg_Loc=%Source_Loc%\RegistryFiles
SET NWTC_Lib_Loc=%HD_Reg_Loc%\include




IF /I "%2"=="bjonkman" CALL ..\Set_FAST_paths.bat

SET ModuleName=%1


REM ----------------------------------------------------------------------------
REM ---------------- RUN THE REGISTRY TO AUTO-GENERATE FILES -------------------
REM ----------------------------------------------------------------------------




%REGISTRY% "%HD_Reg_Loc%\%ModuleName%.txt" -I %NWTC_Lib_Loc% -I %HD_Reg_Loc%
GOTO checkError


:checkError
ECHO.
IF %ERRORLEVEL% NEQ 0 (
ECHO Error running  Registry for SubDyn.
) ELSE (
ECHO %ModuleName%_Types.f90 was created.
COPY /Y "%ModuleName%_Types.f90" "%Source_Loc%"
)




:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ----------------------------------------------------------------------------
ECHO. 


SET REGISTRY=

SET NWTC_Lib_Loc=
SET Source_Loc=

SET CURR_LOC=
:Done
echo %lines%
set lines=