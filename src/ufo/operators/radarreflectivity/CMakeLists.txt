##########################
add_subdirectory( directZDA )
add_subdirectory( genericreflectivity )
add_subdirectory( ppro )

PREPEND( _p_directZDA_files     "operators/radarreflectivity/directZDA"     ${directZDA_src_files} )
PREPEND( _p_genericreflectivity_files     "operators/radarreflectivity/genericreflectivity"     ${genericreflectivity_src_files} )
PREPEND( _p_ppro_files     "operators/radarreflectivity/ppro"     ${ppro_src_files} )
PREPEND( _p_microphysics_file     "operators/radarreflectivity"     "MicrophysicsOptions.h" ) # K.R. - Shared microphysics options

set ( radarreflectivity_src_files
        ${_p_microphysics_file}
        ${_p_directZDA_files}
        ${_p_genericreflectivity_files}
	${_p_ppro_files}
    PARENT_SCOPE  
)
