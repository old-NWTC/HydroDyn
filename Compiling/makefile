#================================================================================#
# This makefile created by B. Jonkman on 2-Apr-2013,                             #
# adapted from Crunch (M. Buhl on 25-Jan-2013).                                  #
# (c) 2013 National Renewable Energy Laboratory                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
                                                                                #
# Older versions of the source code may not work with this makefile.             #
#================================================================================#

   # 32-bit or 64-bit?
BITS = 32
#BITS = 64


   # Location of source files for FAST, AeroDyn, InflowWind, and the NWTC Library.
   # You will probably need to change these for your system.


ifeq ($(OS),Windows_NT)
   NWTC_Lib_DIR = C:/Users/bjonkman/Documents/DATA/DesignCodes/miscellaneous/nwtc_subs/SVNdirectory/trunk/source
   HD_DIR       = C:/Users/bjonkman/Documents/DATA/DesignCodes/simulators/HydroDyn/SVNdirectory/branches/HydroDyn_Modularization/Source
   HD_Reg_DIR   = C:/Users/bjonkman/Documents/DATA/DesignCodes/simulators/HydroDyn/SVNdirectory/branches/HydroDyn_Modularization/Source/RegistryFiles

   Registry     = REGISTRY
else
   NWTC_Lib_DIR = $(HOME)/PC/CAEtools/miscellaneous/nwtc_subs/SVNdirectory/trunk/source
   HD_DIR       = $(HOME)/PC/CAEtools/simulators/HydroDyn/SVNdirectory/branches/HydroDyn_Modularization/Source
   HD_Reg_DIR   = $(HOME)/PC/CAEtools/simulators/HydroDyn/SVNdirectory/branches/HydroDyn_Modularization/Source/RegistryFiles

   Registry     = registry.exe

endif


   # Name of compiler to use and flags to use.
FC      = gfortran
FFLAGS  = -O2 -m$(BITS) -fbacktrace -finit-local-zero -fno-automatic -ffree-line-length-none -x f95-cpp-input
LDFLAGS = -O2 -m$(BITS) -fbacktrace

# other useful gfortran options:
#  -fdefault-real-8 -fcheck=bounds -std=f2003 -O0 -v


   # Precision.

# Use "SingPrec" for single precision and "DoubPrec" for double precision. You may also need to change an option switch to make constants DP.
PREC = SingPrec



   # Destination and RootName for executable

OUTPUT_NAME = HydroDyn
DEST_DIR    = .

   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)
      # Windows
   DEL_CMD   = del
   EXE_EXT   = _gwin$(BITS).exe
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else
      # Linux
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS)
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux
endif


   # Source files (by module)

LIB_SOURCES =           \
	$(PREC).f90          \
	NWTC_Base.f90        \
	$(SYS_FILE).f90      \
	NWTC_IO.f90          \
	NWTC_Num.f90         \
	ModMesh_Types.f90    \
	ModMesh.f90          \
	ModMesh_Mapping.f90  \
	NWTC_Library_Types.f90 \
	NWTC_Library.f90

HD_SOURCES =                 \
	FFTMod.f90                \
	SS_Radiation_Types.f90    \
	SS_Radiation.f90          \
	Waves_Types.f90           \
	Waves.f90                 \
	Current_Types.f90         \
	Current.f90               \
	Morison_Types.f90         \
	Morison_Output.f90        \
	Morison.f90               \
	Conv_Radiation_Types.f90  \
	Conv_Radiation.f90        \
	WAMIT_Types.f90           \
	WAMIT_Output.f90          \
	WAMIT.f90                 \
	WAMIT_Interp.f90			  \
	HydroDyn_Output.f90       \
	HydroDyn_Types.f90        \
	HydroDyn_Input.f90        \
	HydroDyn.f90              \
	HydroDyn_DriverCode.f90


vpath %.f90 $(NWTC_Lib_DIR) $(HD_DIR) $(HD_Reg_DIR)
vpath %.f   $(HD_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)
vpath %.txt $(HD_Reg_DIR)

LIB_OBJS      = $(LIB_SOURCES:.f90=.obj)
HD_OBJS       = $(HD_SOURCES:.f90=.obj) fftpack.obj
ALL_OBJS      = $(LIB_OBJS) $(HD_OBJS)

   # Rule to do everything.

all:     default
default:  $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)


   # General rules for compliling the files.

%.obj: %.f90
	$(FC)  -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)

fftpack.obj: fftpack.f
	$(F77) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)




   #  Dependency rules.
#NWTC Library dependency rules:
NWTC_Base.obj:           $(PREC).obj
$(SYS_FILE).obj:         NWTC_Base.obj
NWTC_Library_Types.obj:  $(SYS_FILE).obj
NWTC_IO.obj:             NWTC_Library_Types.obj
NWTC_Num.obj:            NWTC_IO.obj
ModMesh_Types.obj:       NWTC_Num.obj
ModMesh.obj:             ModMesh_Types.obj
ModMesh_Mapping.obj:     ModMesh.obj
NWTC_Library.obj:        ModMesh.obj  ModMesh_Mapping.obj


#HydroDyn dependency rules:
fftpack.obj:                NWTC_Library.obj
FFTMod.obj:                 NWTC_Library.obj  fftpack.obj

SS_Radiation_Types.obj:     NWTC_Library.obj
SS_Radiation.obj:           SS_Radiation_Types.obj

Waves_Types.obj:            NWTC_Library.obj
Waves.obj:                  Waves_Types.obj  FFTMod.obj

Current_Types.obj:          NWTC_Library.obj
Current.obj:                Current_Types.obj

Morison_Types.obj:          NWTC_Library.obj
Morison_Output.obj:         Morison_Types.obj Waves.obj
Morison.obj:                Morison_Types.obj Morison_Output.obj

Conv_Radiation_Types.obj:   NWTC_Library.obj
Conv_Radiation.obj:         Conv_Radiation_Types.obj

WAMIT_Types.obj:            NWTC_Library.obj Conv_Radiation_Types.obj SS_Radiation_Types.obj Waves_Types.obj
WAMIT.obj:                  WAMIT_Types.obj WAMIT_Output.obj Waves_Types.obj Conv_Radiation.obj SS_Radiation.obj FFTMod.obj WAMIT_Interp.obj

WAMIT_Interp.obj:				 NWTC_Library.obj

HydroDyn_Output_Types.obj:  NWTC_Library.obj
WAMIT_Output.obj:           WAMIT_Types.obj Waves.obj
HydroDyn_Output.obj:        HydroDyn_Types.obj Waves.obj

HydroDyn_Types.obj:         NWTC_Library.obj  HydroDyn_Output_Types.obj  Current_Types.obj  Waves_Types.obj  Conv_Radiation_Types.obj SS_Radiation_Types.obj WAMIT_Types.obj Morison_Types.obj
HydroDyn_Input.obj:         HydroDyn_Types.obj  Waves.obj Morison.obj
HydroDyn.obj:               HydroDyn_Types.obj  HydroDyn_Input.obj  HydroDyn_Output.obj  WAMIT.obj  Current.obj  Morison.obj


# HydroDyn DriverCode dependency rules:
HydroDyn_DriverCode.obj:              NWTC_Library.obj  HydroDyn.obj


   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # Rules for creating _Types.f90 files from Registry files.

Current_Types.f90:
	$(Registry) ./RegistryFiles/Current.txt


Waves_Types.f90:
	$(Registry) ./RegistryFiles/Waves.txt


SS_Radiation_Types.f90:
	$(Registry) ./RegistryFiles/SS_Radiation.txt


Conv_Radiation_Types.f90:
	$(Registry) ./RegistryFiles/Conv_Radiation.txt


WAMIT_Types.f90:
	$(Registry) ./RegistryFiles/WAMIT.txt


Morison_Types.f90:
	$(Registry) ./RegistryFiles/Morison.txt


HydroDyn_Types.f90:
	$(Registry) ./RegistryFiles/HydroDyn.txt




   # For linking HydroDyn Driver Code.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj Current_Types.f90 Waves_Types.f90 SS_Radiation_Types.f90 Conv_Radiation_Types.f90 WAMIT_Types.f90 Morison_Types.f90 HydroDyn_Types.f90
