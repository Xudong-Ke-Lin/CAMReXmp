# This should point to the location of the amrex source code.
AMREX_HOME = ../../amrex
TOPTOP := .

USE_EB =FALSE
PRECISION  = DOUBLE
PROFILE    = FALSE

# Debug mode compiles without optimisation, so a debugger returns
# accurate line information, without optimising out variables.
# Additionally, additional bounds checking is implemented in AMReX's
# debug mode
#DEBUG      = TRUE
DEBUG      = FALSE

# Dimensionality is specified at compile time - this example has been written for 2D only
DIM        = 2

# Compiler can be specified here
COMP	   = gnu

# MPI and OpenMP compatibility must be included at compile time.  OpenMP will also require appropriate pragmas.
USE_MPI    = TRUE
USE_OMP    = FALSE

Bpack   := $(TOPTOP)/Make.package 
Blocs   := $(TOPTOP) 

# Include Eigen library
INCLUDE_LOCATIONS += /usr/include/eigen3/

TOP := $(TOPTOP)/..

EBASE := CAMREXmp

TmpBuildDir = tmp_build_dir_$(EBASE)
include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Bdirs 	:= Source Source/Src_nd
Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs   += $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

Pdirs 	:= Base Boundary AmrCore Amr Particle
Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

# to use linear solvers
Pdirs   += LinearSolvers/MLMG

include $(Ppack)

all: $(executable) 
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


