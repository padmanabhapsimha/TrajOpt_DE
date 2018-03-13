%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Description of the input file for differential evolution code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1>  Number of generations to solve for
2>  Dimensionality of problem
3>  Total number of agents
4>  CR value
5>  F value
6>  Selector string - details below
7>  Lower bound of variables to randomly spawn
8>  Upper bound of variables to randomly spawn
9>  Problem cost function type - details below
10> Number of threads to use to solve problem
11> Seed for random number generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Selector string
carry_toAll -> only provide 1 upper and lower bound and these are used over all dimensions
anything else -> need to provide equal number of upper and lower bounds as per dimension number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Problem cost function listing
1 - Sphere function - verified
2 - 2D Ackerley function - verified
3 - nD Rosenbrock function - verified
4 - Easom function - verified
5 - Schaffer function N.2 - verified
6 - Styblinski-Tang function - verified
7 - Eggholder function - verified
8 - Himmelblau function - verified
9 - nD Rastrigin function - verified
10 - Schwefel function - verified
11 - Cross in tray function - verified
12 - Chichinadze function - verified
13 - Cube function - verified
14 - Damavandi function - verified - difficult
15 - Holder table function - verified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
16 - Optimal time car acceleration profile - verified (3D problem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
17 - 2D solar electric propulsion spacecraft optimal control - verified
18 - 2D nuclear electric propulsion spacecraft optimal control - verified
19 - 2D continuous thrusting minimum effort optimal control - working
20 - 2D continuous thrust tangential to achieve radius only - working
21 - 2D continuous thrust effort optimal ax, ay formulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
22 - 3D continuous thrusting minimum effort optimal control
23 - 2D continuous thrusting fuel optimal no intermediate state formulation
24 - 2D continuous thrusting r based integration termination min fuel
25 - 3D continuous thrusting fuel optimal general code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%