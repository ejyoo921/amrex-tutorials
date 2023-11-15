AMREX_HOME ?= ../../../../amrex

# We advise you to add the locations of your CGAL and boost.
INCLUDE_LOCATIONS += /home/eyoo/spack/opt/spack/linux-rhel7-skylake_avx512/gcc-8.4.0/cgal-5.5.2-774me4awmah35ndfwa6h2cdjdxbi7wip/include/
INCLUDE_LOCATIONS += /home/eyoo/spack/opt/spack/linux-rhel7-skylake_avx512/gcc-8.4.0/boost-1.83.0-lm4clgzwoav4coynyuhxj7zr5pt2jhnb/include/

PROFILE      = TRUE # Full profiling
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
