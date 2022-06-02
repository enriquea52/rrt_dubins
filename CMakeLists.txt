cmake_minimum_required(VERSION 3.0.2)
project(rrt_star_planner)

find_package(catkin REQUIRED COMPONENTS
    pybind11_catkin)

catkin_package(
    LIBRARIES
    CATKIN_DEPENDS pybind11_catkin
)

include_directories(${catkin_INCLUDE_DIRS})

pybind_add_module(rrt_star_planner MODULE modules/main.cpp)


catkin_install_python(PROGRAMS scripts/main.py
                      DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

install(TARGETS rrt_star_planner
  LIBRARY DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
