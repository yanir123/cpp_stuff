cmake_minimum_required (VERSION 3.0)
project (myproject)
 
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 
add_executable (dmrg Mpo.cpp Main.cpp)
target_link_libraries (dmrg Eigen3::Eigen)
