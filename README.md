# TrajOpt_DE
SOME FILES OF THE CODE ARE RATHER DEVOID OF DOCUMENTATION (THESE ARE RELATED TO THE ACTUAL SPACECRAFT STATE AND COSTATE DYNAMICS)
THE REST SEEMS REASONABLY OK.
READ THE INDIVIDUAL README FILES FOR INPUTS.
EDIT THE MAKEFILE TO POINT TO THE ACTUAL LOCATION OF THE BOOST LIBRARY IN THE SYSTEM. SPECIFICALLY THE $(INC) THINGY.
SORRY FOR SHOUTING..!! THE ACTUAL README FOLLOWS THIS.
# Actual readme stuff - maybe
Spacecraft trajectory optimization using differential evolution (DE).
At the heart of this code lies a multiple thread capable differential evolution based optimizer. 
Coarse grain parallelized by decomposition of initial population.
Run in distributed systems by manual decomposition of domains.
Wide variety of benchmark problems hard coded in a long switch-case construct. Eliminate this
# cpp files
de_solver.cpp has the actual DE solver coded in (Refer Storn and Price for more details). Knuth shuffle is used as a nifty little trick here for making the implementation much cleaner. de_func.cpp generates the initial population and some other stuff. Other files are mostly for support and modularity. cost_fn.cpp and cost_fn_print.cpp are absolutely bloated and can be trimmed out if parts are not necessary.
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
