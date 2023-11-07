AMREX_HOME ?= ../../../../amrex
INCLUDE_LOCATIONS += /usr/local/opt/cgal/lib/cmake/CGAL/

PROFILE      = FALSE # Full profiling
TINY_PROFILE = FALSE

DIM	= 3

COMP    = gnu

USE_EB    = TRUE
USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE

include ./Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
