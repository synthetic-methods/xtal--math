cmake_minimum_required(VERSION 3.14)

project(xtal--math VERSION 0.0.0)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)


find_package(xtal CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} INTERFACE xtal::xtal)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
