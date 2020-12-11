############################################
#
#    dolfin.cpp.la
#
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/cpp/la/index.html
#
############################################

@pyclass dolfin Matrix FeObject FeMatrix
@pyclass dolfin Vector FeObject FeVector
@pyclass dolfin PETScVector FeVector
@pyclass dolfin PETScMatrix FeMatrix
@pyclass dolfin EigenVector FeVector
@pyclass dolfin EigenMatrix FeMatrix
@pyclass dolfin LUSolver
@pyclass dolfin SLEPcEigenSolver
@pyclass dolfin KrylovSolver
@pyclass dolfin NewtonSolver
@pyclass dolfin NonlinearProblem

@pyfunc ufl as_tensor
@pyfunc ufl as_vector
@pyfunc ufl as_matrix

@pyfunc dolfin list_linear_algebra_backends
@pyfunc dolfin has_linear_algebra_backend
@pyfunc dolfin has_slepc
@pyfunc dolfin info

export FeMatrix, FeVector, PETScVector, PETScMatrix, EigenVector, EigenVector,
       LUSolver, SLEPcEigenSolver, NewtonSolver, NonlinearProblem, KrylovSolver, 
       as_tensor, as_vector, as_matrix, list_linear_algebra_backends, 
       has_linear_algebra_backend, has_slepc, info


