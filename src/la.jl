############################################
#
#    dolfin.cpp.la
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/cpp/la/index.html
#
#    dolfin.la
# https://fenicsproject.org/docs/dolfin/latest/python/_autogenerated/dolfin.la.html
#
############################################

@pyclass dolfin Matrix FeObject FeMatrix
@pyclass dolfin Vector FeObject FeVector
@pyclass dolfin PETScVector FeVector
@pyclass dolfin PETScMatrix FeMatrix
@pyclass dolfin EigenVector FeVector
@pyclass dolfin EigenMatrix FeMatrix
@pyclass dolfin VectorSpaceBasis
@pyclass dolfin LUSolver
@pyclass dolfin SLEPcEigenSolver
@pyclass dolfin KrylovSolver
@pyclass dolfin NewtonSolver
@pyclass dolfin PETScKrylovSolver
@pyclass dolfin NonlinearProblem

@pyfunc ufl as_tensor
@pyfunc ufl as_vector
@pyfunc ufl as_matrix

@pyfunc dolfin list_linear_algebra_backends
@pyfunc dolfin has_linear_algebra_backend
@pyfunc dolfin info
@pyfunc dolfin has_krylov_solver_method
@pyfunc dolfin has_krylov_solver_preconditioner

@pyfunc dolfin as_backend_type
@pyfunc dolfin la_index_dtype

*(v::Union{FeMatrix, FeVector}, c::Real) = to_fetype(v.pyobject * c)
*(c::Real, v::Union{FeMatrix, FeVector}) = to_fetype(c * v.pyobject)

+(v1::FeVector, v2::FeVector) = to_fetype(v1.pyobject + v2.pyobject)

-(v1::FeVector, v2::FeVector) = to_fetype(v1.pyobject - v2.pyobject)

/(v::Union{FeMatrix, FeVector}, c::Real) = to_fetype(v.pyobject / c)

export FeMatrix, FeVector, PETScVector, PETScMatrix, EigenVector, EigenMatrix, 
       VectorSpaceBasis, LUSolver, SLEPcEigenSolver, NewtonSolver, 
       PETScKrylovSolver, NonlinearProblem, KrylovSolver, as_tensor, as_vector, 
       as_matrix, list_linear_algebra_backends, has_linear_algebra_backend, 
       info, has_krylov_solver_method,has_krylov_solver_preconditioner,
       as_backend_type, la_index_dtype


