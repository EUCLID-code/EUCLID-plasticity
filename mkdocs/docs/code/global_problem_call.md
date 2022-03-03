# `GLOBAL_PROBLEM_CALL` - executes the finite element analysis
###  Comments
In this function, the finite element analysis is executed. Dependent on
the type of material (e.g., type of hardening), different finite element
solvers are called.

###  Input Arguments
`algorithm` (_struct_) - structural array containing information needed
for the algorithm

`material` (_struct_) - structural array containing material information

`mesh` (_struct_) - structural array containing information about the
finite element mesh

`bc` (_struct_) - structural array containing information about the
boundary conditions

###  Output Arguments
`results` (_struct_) - structural array containing finite element results

