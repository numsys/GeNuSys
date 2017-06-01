
if(GeNuSys_INCLUDE_DIRS AND GeNuSys_LIBRARIES)
    # Already in cache, be silent
    set(GeNuSys_FIND_QUIETLY TRUE)
endif()

find_path(GeNuSys_INCLUDE_DIRS GeNuSys/number_system.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GeNuSys DEFAULT_MSG GeNuSys_INCLUDE_DIRS)

mark_as_advanced(GeNuSys_INCLUDE_DIRS)
