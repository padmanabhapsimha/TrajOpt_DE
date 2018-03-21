# TrajOpt_DE
This code can optimize continuous low thrust spacecraft trajectories by solving the two point boundary value problem that results from the application of the Pontryagin's minimum principle. At the heart of this code lies a multiple thread capable optimization routine based on Differential Evolution (Refer to Storn and Price for details). This version of the code is coarse grain parallelized by domain decomposition. This has been designed for shared memory parallel systems. Can be run in distributed systems by manual decomposition of the search hyperbox. This is a full 3D code and is quite capable of handling arbitrary thruster performance characteristics. Nuclear electric and inverse square solar electric thrust models have been hard coded in. Other exotic thrusters can be manually put in along with stuff like eclipsing of the Sun etc. This code can as of now handle transfers between arbitrary parking orbits of different planets. (Earth-Mars transfer from 25000km EPO at -23.4 degrees to Earth equator to 17500km MPO at 90 degrees to the ecliptic plane has been verified)
GTO to GSO transfers have also been verified and have been compared with results from literature. Better optimal results have been obtained. 3D to 3D heliocentric transfer portion of this code as of now is the most capable and tested part.
# Actual readme stuff - Important info
Read individual input file readme(s) for supplying the inputs.
Boost library header files are absolutely essential. Version 1.66.0 has been used. 
Wide variety of benchmark problems hard coded in a long switch-case construct. Eliminate this if not needed.
# cpp files
de_solver.cpp has the parallelized version of the DE solver coded. Durstenfeld's version of the Fisher-Yates shuffle is used as a nifty little trick here for making the implementation much cleaner. de_func.cpp generates the initial population and some other stuff. Other files are mostly for support and modularity. cost_fn.cpp and cost_fn_print.cpp are absolutely bloated and can be trimmed out if parts are not necessary.
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
