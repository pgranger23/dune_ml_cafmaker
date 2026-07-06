cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 20)


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
	    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
	      
	      # Set the possible values for cmake-gui drop-down
	        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
	endif()


# C++ ABI is set per-target/directory to match specific dependency requirements (HDF5 vs ROOT)

project(dune_ml_cafmaker)

# Robust HDF5 detection
if(HDF5_INSTALL)
  set(HDF5_ROOT ${HDF5_INSTALL})
  set(HDF5_DIR ${HDF5_INSTALL})
endif()

# Explicitly request CXX and HL components
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

# Handle the case where HDF5_INCLUDE_DIRS is empty (common in CONFIG mode)
if("${HDF5_INCLUDE_DIRS}" STREQUAL "")
  get_filename_component(_HDF5_BASE "${HDF5_DIR}/../../../../" ABSOLUTE)
  if(NOT EXISTS "${_HDF5_BASE}/include/H5Cpp.h")
    get_filename_component(_HDF5_BASE "${HDF5_DIR}/../../../" ABSOLUTE)
  endif()
  if(EXISTS "${_HDF5_BASE}/include/H5Cpp.h")
    set(HDF5_INCLUDE_DIRS "${_HDF5_BASE}/include")
  endif()
endif()

# Ensure we have the C++ and HL libraries
if(TARGET hdf5::hdf5_cpp-shared)
    set(HDF5_LIBRARIES hdf5::hdf5_cpp-shared hdf5::hdf5_hl_cpp-shared hdf5::hdf5-shared)
elseif(TARGET hdf5::hdf5_cpp-static)
    set(HDF5_LIBRARIES hdf5::hdf5_cpp-static hdf5::hdf5_hl_cpp-static hdf5::hdf5-static)
endif()

# Debugging info
message(STATUS "HDF5 Found: ${HDF5_FOUND}")
message(STATUS "HDF5 Include Dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 Libraries: ${HDF5_LIBRARIES}")

# Build local duneanaobj from source
add_subdirectory(duneanaobj)
if(TARGET duneanaobj_StandardRecord)
    add_library(duneanaobj::StandardRecord ALIAS duneanaobj_StandardRecord)
elseif(TARGET StandardRecord)
    add_library(duneanaobj::StandardRecord ALIAS StandardRecord)
endif()
find_package(ZLIB REQUIRED)
find_package(ROOT REQUIRED)

file(GLOB SPINE_SOURCES src/*.cc)

# Helper function to ensure all includes and libraries are correctly attached to targets
function(link_dlp_dependencies target_name)
    target_include_directories(${target_name} PUBLIC 
        include 
        ${HDF5_INCLUDE_DIRS} 
        ${ROOT_INCLUDE_DIRS}
    )
    target_link_libraries(${target_name} PRIVATE 
        ${HDF5_LIBRARIES} 
        ZLIB::ZLIB 
        duneanaobj::StandardRecord
    )
endfunction()

# Libraries
add_library(dlp_simulation SHARED ${SPINE_SOURCES})
link_dlp_dependencies(dlp_simulation)
target_compile_definitions(dlp_simulation PRIVATE MC_NOT_DATA)

add_library(dlp_data SHARED ${SPINE_SOURCES})
link_dlp_dependencies(dlp_data)

# Helper for executables
function(link_executable_dependencies target_name)
    link_dlp_dependencies(${target_name})
    target_link_libraries(${target_name} PRIVATE ${ROOT_LIBRARIES})
endfunction()

# Merge executables
add_executable(merge_sources_simulation merge_sources.cc)
target_compile_definitions(merge_sources_simulation PRIVATE MC_NOT_DATA)
link_executable_dependencies(merge_sources_simulation)
target_link_libraries(merge_sources_simulation PRIVATE dlp_simulation)

add_executable(merge_sources_data merge_sources.cc)
link_executable_dependencies(merge_sources_data)
target_link_libraries(merge_sources_data PRIVATE dlp_data)

# Standalone executables
add_executable(make_standalone_simulation make_standalone.cc)
target_compile_definitions(make_standalone_simulation PRIVATE MC_NOT_DATA)
link_executable_dependencies(make_standalone_simulation)
target_link_libraries(make_standalone_simulation PRIVATE dlp_simulation)

add_executable(make_standalone_data make_standalone.cc)
link_executable_dependencies(make_standalone_data)
target_link_libraries(make_standalone_data PRIVATE dlp_data)
