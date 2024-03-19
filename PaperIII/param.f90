! The external parameters of the solver - in general only the point cloud 
! of the reference triangle. Possible choices are CG1 for linear elements
! and CG3 for cubic elements. Please note that while continuous Galerkin
! node distributions are used is it still a discontinuous galerkin solver.
module param
        implicit none
        character(len=*), parameter :: OpDesign = "CG1"
end module param
