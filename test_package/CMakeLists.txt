cmake_minimum_required (VERSION 3.20)
project (test_package)

set (CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package (xtal--math CONFIG REQUIRED)
find_package (Catch2    CONFIG REQUIRED)
find_package (Eigen3    CONFIG REQUIRED)

cmake_path (GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH xtal--math_SOURCE_DIR)

if (EXISTS ${xtal--math_SOURCE_DIR}/include/FIXME.cc)
	add_executable (test ${xtal--math_SOURCE_DIR}/include/FIXME.cc)
else ()
	add_executable (test ${xtal--math_SOURCE_DIR}/include/xtal--math/all.cc)
endif ()


if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
#
else ()
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffp-contract=fast")
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffinite-math-only")
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -freciprocal-math")
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-math-errno")
#	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffp-model=fast")
#	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fapprox-func")
#	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
#	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
#	set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -m64")
endif ()

target_link_libraries (test xtal--math::xtal--math Catch2::Catch2WithMain Eigen3::Eigen)

target_precompile_headers (test INTERFACE xtal--math::xtal--math Catch2::Catch2WithMain Eigen3::Eigen)
