cmake_minimum_required(VERSION 3.21)
project(TraceBoxes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(trace_boxes src/main.cpp src/image_io.cpp src/traces.cpp)
# The program is build/<preset>/trace_boxes on every platform (no Release/ subfolder with Visual Studio).
set_target_properties(trace_boxes PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$<1:${CMAKE_BINARY_DIR}>")

# `cmake --workflow --preset run`: configure, build and run in one command (the file dialog opens).
add_custom_target(run COMMAND trace_boxes WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}" USES_TERMINAL)
if(MSVC)
  target_compile_options(trace_boxes PRIVATE /W4 /utf-8)
  target_compile_definitions(trace_boxes PRIVATE NOMINMAX _CRT_SECURE_NO_WARNINGS)
  target_link_libraries(trace_boxes PRIVATE comdlg32)
else()
  target_compile_options(trace_boxes PRIVATE -Wall -Wextra)
endif()

# `ctest`: run on every TIFF in the project folder.
enable_testing()
file(GLOB SAMPLES "${CMAKE_SOURCE_DIR}/*.tif" "${CMAKE_SOURCE_DIR}/*.tiff")
if(SAMPLES)
  add_test(NAME samples COMMAND trace_boxes -o "${CMAKE_BINARY_DIR}/test_output" ${SAMPLES})
endif()
