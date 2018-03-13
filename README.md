# TrajOpt_DE
Spacecraft trajectory optimization using differential evolution (DE).
At the heart of this code lies a multiple thread capable differential evolution based optimizer. 
Coarse grain parallelized by decomposition of initial population.
Run in distributed systems by manual decomposition of domains.
Wide variety of benchmark problems hard coded in a long switch-case construct. Eliminate this
# Trajectory Optimization
Indirect optimal control theory with the application of Pontryagin's minimum principle.
Resulting two point boundary value problem (TPBVP) solved using the DE optimizer.
Full 3D trajectories have been generated.
# Other related stuff
64bit C++17 standard. Originally compiled with g++ 7.2.0
STL and Boost libraries have been heavily used.
Check out Project Pluto --> https://www.projectpluto.com/jpl_eph.htm#eph_basics is amazing for reading JPL ephemeris data. 
Github page --> https://github.com/Bill-Gray/jpl_eph
Some stupid looking/sounding wrappers have been written to make the ephemeris code C++ friendly and hope to keep up with RAII.
# Additional comments
This code is quite fast..!! 
Multithreading scales well with cores. Superlinear speedups have also been observed when threads are reasonable in number. 
It's yet to be profiled though. Main emphasis till now has been on the accuracy of the numerical results obtained.
# Future work
Make this MPI capable instead of only the C++17 threading framework.
