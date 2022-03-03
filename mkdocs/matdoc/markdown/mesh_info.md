# `MESH_INFO` - generates important information for a given finite element mesh
###  Comments
Given a finite element mesh, this function computes information such as
the number of degrees of freedom, the number of elements, the Jacobian of
the mapping between physical and reference element, the shape function
derivatives etc..

###  Input Arguments
`mesh` (_struct_) - structural array containing information about the
finite element mesh

###  Output Arguments
`mesh` (_struct_) - structural array containing information about the
finite element mesh such as for example the determinant of the Jacobian
`mesh.detJ_GP` and the shape function derivatives `mesh.B_GP` at each
Gauss point

