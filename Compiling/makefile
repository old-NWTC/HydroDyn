#================================================================================#
# This makefile created by B. Jonkman on 20-May-2015,                            #
# adapted from FAST v8.                                                          #
# (c) 2015 National Renewable Energy Laboratory                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
#================================================================================#


   # Name of compiler to use and flags to use.
FC      = gfortran

# 32-bit or 64-bit?
BITS = 32
#BITS = 64


FFLAGS  = -O2 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input
LDFLAGS = -O2 -m$(BITS) -fbacktrace


#FFLAGS  = -O0 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input -g -pg
#LDFLAGS = -O0 -m$(BITS) -fbacktrace -Wl,--stack=999999999,--large-address-aware -g -pg

#.........................................
# preprocessor definitions:
#  -DDOUBLE_PRECISION
#
# some useful gfortran options:
#  -DFPE_TRAP_ENABLED
#  -Wl,--large-address-aware   # not necessary when $(BITS) is 64 (64-bit target)
#  -Wl,--stack=999999999       # not necessary when $(BITS) is 64
#  -Wconversion-extra -Wconversion
#  -fdefault-real-8 -fcheck=bounds,do,mem,pointer -std=f2003 -O0 -v -Wall
#  -pg                          # generate debugging info for debugger
# http://gcc.gnu.org/onlinedocs/gfortran/Option-Index.html#Option-Index
#
# makefile options:
#  --warn-undefined-variables
#.........................................


   # Location of source files for FAST, and its modules.

MAIN_DIR     = ../Source
BIN_DIR      = ../bin

NWTC_Lib_DIR = $(MAIN_DIR)/dependencies/NWTC_Library
NETLIB_DIR   = $(MAIN_DIR)/dependencies/NetLib
ED_DIR       = $(MAIN_DIR)/dependencies/ElastoDyn
HD_DIR       = $(MAIN_DIR)

HD_DIR_Reg   = $(HD_DIR)/RegistryFiles

   # Names and locations of the Registry, MAP libraries, and instructions for linking with LAPACK.
   # You will probably need to change these for your system.

ifeq ($(OS),Windows_NT)

   Registry     = $(BIN_DIR)/Registry_win32.exe
   LDFLAGS     := $(LDFLAGS) -llapack -lblas -LC:/LAPACK/win32
   LAPACK_LINK :=

else

   #Registry     = ~/software-development/windsvn/FAST_Registry/registry.exe
   Registry    := registry.exe

   UNAME := $(shell uname -s)

   ifeq ($(UNAME), Darwin)
         # Mac needs some additional flags prepended for linking to get lapack working:
      LDFLAGS     := -framework Accelerate $(LDFLAGS)
      LAPACK_LINK :=
   else
        # ifeq ($(UNAME), Linux)
        # (bjj verify this:)
        # I've been told that the lapack linking needs to occur at the end of the line (after all the object files) on linux:
      LAPACK_LINK  := -llapack -lblas
   endif

endif



   # Destination and RootName for executable

OUTPUT_NAME = HydroDyn
DEST_DIR    = $(BIN_DIR)

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
   CP_CMD    = copy
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else
      # Linux/Mac
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS)
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   CP_CMD    = cp
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux
endif


   # Source files (by module)

LIB_SOURCES =             \
	SingPrec.f90           \
	NWTC_Base.f90          \
	$(SYS_FILE).f90        \
	NWTC_Library_Types.f90 \
	NWTC_IO.f90            \
	NWTC_Num.f90           \
	ModMesh_Types.f90      \
	ModMesh.f90            \
	ModMesh_Mapping.f90    \
	NWTC_Library.f90

NETLIB_SOURCES=            \
	NWTC_FFTPACK.f90        \
	NWTC_LAPACK.f90         \
	fftpack4.1.f


HD_SOURCES =                 \
	SS_Radiation_Types.f90    \
	SS_Radiation.f90          \
	Waves2_Types.f90          \
	Waves2_Output.f90         \
	Waves2.f90                \
	Waves_Types.f90           \
	Waves.f90                 \
	Current_Types.f90         \
	Current.f90               \
	Morison_Types.f90         \
	Morison_Output.f90        \
	Morison.f90               \
	Conv_Radiation_Types.f90  \
	Conv_Radiation.f90        \
	WAMIT2_Types.f90          \
	WAMIT2_Output.f90         \
	WAMIT2.f90                \
	WAMIT_Types.f90           \
	WAMIT_Interp.f90          \
	WAMIT_Output.f90          \
	WAMIT.f90                 \
	HydroDyn_Output.f90       \
	HydroDyn_Types.f90        \
	HydroDyn_Input.f90        \
	HydroDyn.f90              \
	HydroDyn_DriverCode.f90



vpath %.f90       $(NWTC_Lib_DIR) $(NETLIB_DIR) $(HD_DIR)
vpath %.f         $(NETLIB_DIR)
vpath %.mod       $(INTER_DIR)
vpath %$(OBJ_EXT) $(INTER_DIR)
vpath %.txt       $(HD_DIR_Reg)


ALL_SOURCES = $(LIB_SOURCES) $(NETLIB_SOURCES) $(HD_SOURCES)
ALL_OBJS    = $(ALL_SOURCES:.f90=$(OBJ_EXT))
ALL_OBJS   := $(ALL_OBJS:.f=$(OBJ_EXT))


   # Rule to do everything.

all:     default
# use this for "all" for debugging: print variables:
#all:     ; $(info $$HD_DIR_Reg is [${HD_DIR_Reg}]) echo debugging

default:  $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)


   # General rules for compiling the files.

%$(OBJ_EXT): %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)

%$(OBJ_EXT): %.F90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)

%$(OBJ_EXT): %.f
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)



   # General rules for creating _Types.f90 files from Registry files.

   # HydroDyn Registry files:
$(HD_DIR)/%_Types.f90: $(HD_DIR_Reg)/%.txt
	$(Registry) $< -I $(NWTC_Lib_DIR) -I $(HD_DIR_Reg) -O $(HD_DIR)


   #  Dependency rules.
#NWTC Library dependency rules:
NWTC_Base$(OBJ_EXT):              SingPrec$(OBJ_EXT)
$(SYS_FILE)$(OBJ_EXT):            NWTC_Base$(OBJ_EXT)
NWTC_Library_Types$(OBJ_EXT):     $(SYS_FILE)$(OBJ_EXT)
NWTC_IO$(OBJ_EXT):                NWTC_Library_Types$(OBJ_EXT)
NWTC_Num$(OBJ_EXT):               NWTC_IO$(OBJ_EXT)
ModMesh_Types$(OBJ_EXT):          NWTC_Num$(OBJ_EXT)
ModMesh$(OBJ_EXT):                ModMesh_Types$(OBJ_EXT)
ModMesh_Mapping$(OBJ_EXT):        ModMesh$(OBJ_EXT) NWTC_LAPACK$(OBJ_EXT)
NWTC_Library$(OBJ_EXT):           ModMesh$(OBJ_EXT)  ModMesh_Mapping$(OBJ_EXT)

NWTC_LAPACK$(OBJ_EXT):            NWTC_Base$(OBJ_EXT)
NWTC_FFTPACK$(OBJ_EXT):           NWTC_Library$(OBJ_EXT) fftpack4.1$(OBJ_EXT)
fftpack4.1$(OBJ_EXT):             SingPrec$(OBJ_EXT)


#HydroDyn dependency rules:
SS_Radiation_Types$(OBJ_EXT):     NWTC_Library$(OBJ_EXT) $(HD_DIR)/SS_Radiation_Types.f90
SS_Radiation$(OBJ_EXT):           SS_Radiation_Types$(OBJ_EXT)

Waves2_Types$(OBJ_EXT):           NWTC_Library$(OBJ_EXT) $(HD_DIR)/Waves2_Types.f90
Waves2_Output$(OBJ_EXT):          Waves2_Types$(OBJ_EXT)
Waves2$(OBJ_EXT):                 NWTC_Library$(OBJ_EXT)  Waves2_Types$(OBJ_EXT)  NWTC_FFTPACK$(OBJ_EXT)  \
                                     Waves2_Output$(OBJ_EXT)  Waves$(OBJ_EXT)

Waves_Types$(OBJ_EXT):            NWTC_Library$(OBJ_EXT) $(HD_DIR)/Waves_Types.f90
Waves$(OBJ_EXT):                  Waves_Types$(OBJ_EXT)   NWTC_FFTPACK$(OBJ_EXT)

Current_Types$(OBJ_EXT):          NWTC_Library$(OBJ_EXT) $(HD_DIR)/Current_Types.f90
Current$(OBJ_EXT):                Current_Types$(OBJ_EXT)

Morison_Types$(OBJ_EXT):          NWTC_Library$(OBJ_EXT) $(HD_DIR)/Morison_Types.f90
Morison_Output$(OBJ_EXT):         Morison_Types$(OBJ_EXT) Waves$(OBJ_EXT)
Morison$(OBJ_EXT):                Morison_Types$(OBJ_EXT) Morison_Output$(OBJ_EXT)

Conv_Radiation_Types$(OBJ_EXT):   NWTC_Library$(OBJ_EXT) $(HD_DIR)/Conv_Radiation_Types.f90
Conv_Radiation$(OBJ_EXT):         Conv_Radiation_Types$(OBJ_EXT)

WAMIT2_Types$(OBJ_EXT):           NWTC_Library$(OBJ_EXT) $(HD_DIR)/WAMIT2_Types.f90
WAMIT2_Output$(OBJ_EXT):          NWTC_Library$(OBJ_EXT)  WAMIT2_Types$(OBJ_EXT)
WAMIT2$(OBJ_EXT):                 NWTC_Library$(OBJ_EXT)  WAMIT_Interp$(OBJ_EXT)  WAMIT2_Output$(OBJ_EXT)  \
                                    NWTC_FFTPACK$(OBJ_EXT)  Waves$(OBJ_EXT)

WAMIT_Types$(OBJ_EXT):            NWTC_Library$(OBJ_EXT) Conv_Radiation_Types$(OBJ_EXT) SS_Radiation_Types$(OBJ_EXT) Waves_Types$(OBJ_EXT) $(HD_DIR)/WAMIT_Types.f90
WAMIT$(OBJ_EXT):                  WAMIT_Types$(OBJ_EXT) WAMIT_Output$(OBJ_EXT) Waves_Types$(OBJ_EXT) Conv_Radiation$(OBJ_EXT) \
                                    SS_Radiation$(OBJ_EXT) NWTC_FFTPACK$(OBJ_EXT) WAMIT_Interp$(OBJ_EXT)

WAMIT_Interp$(OBJ_EXT):           NWTC_Library$(OBJ_EXT)
WAMIT_Output$(OBJ_EXT):           WAMIT_Types$(OBJ_EXT) Waves$(OBJ_EXT)
HydroDyn_Output$(OBJ_EXT):        HydroDyn_Types$(OBJ_EXT) Waves$(OBJ_EXT)

HydroDyn_Types$(OBJ_EXT):         NWTC_Library$(OBJ_EXT)  Current_Types$(OBJ_EXT)  Waves_Types$(OBJ_EXT)  Waves2_Types$(OBJ_EXT) \
                                     Conv_Radiation_Types$(OBJ_EXT) SS_Radiation_Types$(OBJ_EXT) \
                                     WAMIT_Types$(OBJ_EXT)  WAMIT2_Types$(OBJ_EXT)  Morison_Types$(OBJ_EXT) $(HD_DIR)/HydroDyn_Types.f90
HydroDyn_Input$(OBJ_EXT):         HydroDyn_Types$(OBJ_EXT)  Waves$(OBJ_EXT) Morison$(OBJ_EXT)  HydroDyn_Output$(OBJ_EXT)  \
                                     WAMIT_Output$(OBJ_EXT)  Waves2_Output$(OBJ_EXT)  WAMIT2_Output$(OBJ_EXT)
HydroDyn$(OBJ_EXT):               HydroDyn_Types$(OBJ_EXT)  HydroDyn_Input$(OBJ_EXT)  HydroDyn_Output$(OBJ_EXT)  WAMIT$(OBJ_EXT) \
                                     Current$(OBJ_EXT)  Morison$(OBJ_EXT)  Waves2$(OBJ_EXT)  WAMIT2$(OBJ_EXT)

HydroDyn_DriverCode$(OBJ_EXT):    HydroDyn$(OBJ_EXT)



   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # For linking FAST.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) $(LAPACK_LINK)

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*$(OBJ_EXT)

superclean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*$(OBJ_EXT)
	$(DEL_CMD) "$(DEST_DIR)$(PATH_SEP)$(OUTPUT_NAME)$(EXE_EXT)"
	$(DEL_CMD) "$(HD_DIR)$(PATH_SEP)"*_Types.f90
	$(DEL_CMD) *_Types.f90
