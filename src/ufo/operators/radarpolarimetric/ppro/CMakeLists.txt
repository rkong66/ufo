# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Note: Core physics modules (dualpol_op_mod.f90, dualpol_op_tlad_mod.f90)
# have been moved to the external ppro-lib library

set ( ppro_src_files
    ObsPPRO.h
    ObsPPRO.cc
    ObsPPROTLAD.h
    ObsPPROTLAD.cc
    ObsPPRO.interface.F90
    ObsPPRO.interface.h
    ObsPPROTLAD.interface.F90
    ObsPPROTLAD.interface.h
    ufo_PPRO_mod.F90
    ufo_PPRO_tlad_mod.F90
PARENT_SCOPE
)
