
############################################
#
#    dolfin.cpp.function
# https://fenicsproject.org/docs/dolfin/latest/python/_autogenerated/dolfin.cpp.function.html   
#
#    dolfin.function
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/functions/index.html
#
############################################

@pyclass dolfin Expression
@pyclass ufl Coefficient Expression
@pyclass dolfin Constant Expression
@pyclass dolfin MeshCoordinates Expression
@pyclass dolfin FacetArea Expression

@pyclass dolfin FunctionSpace
@pyclass dolfin Function FeObject FeFunction
@pyclass dolfin MultiMeshFunctionSpace
@pyclass dolfin MultiMeshFunction

@pyfunc dolfin interpolate
@pyfunc dolfin TensorFunctionSpace
@pyfunc dolfin VectorFunctionSpace
@pyfunc dolfin TrialFunction
@pyfunc dolfin TrialFunctions
@pyfunc dolfin TestFunction
@pyfunc dolfin TestFunctions

OpType = Union{Expression, Variable, FeFunction}

export Expression, Coefficient, Constant, FeFunction, FunctionSpace, MultiMeshFunction, 
       FacetArea, MeshCoordinates, interpolate, VectorFunctionSpace, MultiMeshFunctionSpace, 
       TensorFunctionSpace, TrialFunction, TrialFunctions, TestFunction, TestFunctions


############################################
#
#    ufi Expression
#
# https://fenicsproject.org/docs/ufl/1.6.0/ufl.html
#
############################################

# no need to export
# not exported in `dolfin` or `ufl`
@pyclass dolfin Equation Expression

# export 
@pyclass ufl Argument Expression
@pyclass dolfin Form Expression
@pyclass ufl Measure
@pyclass ufl Identity Expression

array(vec::GenericVector) = vec.pyobject.gather_on_zero()
array(vec::GenericMatrix) = vec.pyobject.array()
array(form::Expression) =  array(assemble(form)).gather_on_zero()
function array(solution::FeFunction)
    generic_vector = solution.pyobject.vector()
    instantiated_vector = dolfin.Vector(generic_vector)
    instantiated_vector.gather_on_zero()
end

len(U::OpType) = length(U.pyobject)

export Argument, Form, Measure, Identity, array, len


