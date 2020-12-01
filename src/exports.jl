# must be explicitly imported to be extended
import Base: getproperty, *, +, -, /, ^, ==, <<, sqrt, div, split

@pyclass fenics Matrix FeObject FeMatrix
@pyclass fenics Vector FeObject FeVector
@pyclass fenics PETScVector FeVector
export FeMatrix, FeVector, PETScVector


############################################
#
#    dolfin.cpp.function
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.function.html   
#
#    dolfin.function
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/functions/index.html
#
############################################

@pyclass fenics Expression
@pyclass fenics Constant Expression
@pyclass ufl FacetNormal Expression
@pyclass fenics Function FeObject FeFunction
@pyclass fenics FunctionSpace

@pyfunc fenics interpolate
@pyfunc fenics VectorFunctionSpace
@pyfunc fenics TrialFunction
@pyfunc fenics TestFunction
@pyfunc fenics TestFunctions

export Constant, Expression, FeFunction, FunctionSpace, interpolate, VectorFunctionSpace, TrialFunction, TestFunction, TestFunctions, FacetNormal

############################################
#
#    dolfin.mesh
#    dolfin.cpp.mesh module
#
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.mesh.html
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/cpp/mesh/index.html
#
############################################

@pyclass fenics Mesh
@pyclass fenics MeshFunction

# no need to export
@pyclass fenics MeshDomains

export Mesh, MeshFunction

############################################
#
#    dolfin.cpp.generation
#
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.generation.html
#
############################################

@pyclass fenics BoxMesh
@pyclass fenics RectangleMesh
@pyclass fenics UnitSquareMesh

export BoxMesh, RectangleMesh, UnitSquareMesh

############################################
#
#    dolfin.cpp.fem
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.fem.html
#
#    dolfin.fem module
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/fem/index.html
#
############################################

@pyclass fenics DirichletBC

@pyfunc fenics assemble
@pyfunc fenics errornorm
@pyfunc fenics project
@pyfunc fenics solve

#import Base: split  # must be explicitly imported to be extended

function split(fun::FeFunction)
     vec = fenics.split(fun.pyobject)
     expr_vec = [FeFunction(spl) for spl in vec]
     return expr_vec
end

export DirichletBC, assemble, errornorm, project, solve

############################################
#
#    dolfin.cpp.geometry
#
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.geometry.html
#
############################################

@pyclass fenics Point
export Point

############################################
#
#    dolfin.common
#
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/common/index.html
#
############################################

@pyfunc fenics plot
export plot


############################################
#
#    ufi      
#
# https://fenicsproject.org/docs/ufl/1.6.0/ufl.html
# https://fenicsproject.org/pub/documents/ufl/ufl-user-manual/ufl-user-manual.pdf
#
############################################
@pyclass ufl Cell
@pyclass fenics FiniteElement
@pyclass fenics MixedElement

export FiniteElement, MixedElement

#import Base: div, sqrt  # must be explicitly imported to be extended

@pyclass fenics Measure
@pyclass fenics Identity Expression

# no need to export
@pyclass fenics Form Expression
@pyclass ufl Argument Expression

# no need to export
# not exported in `feincs` and `ufl` 

@pyclass ufl Equation Expression

# < Operator < Expr
@pyclass ufl ListTensor Expression
@pyclass ufl ComponentTensor Expression
@pyclass ufl Indexed Expression

# < CompoundTensorOperator < Operator < Expr
@pyclass ufl Transposed Expression
@pyclass ufl Inner Expression
@pyclass ufl Dot Expression
@pyclass ufl Sym Expression

# < CompoundDerivative < Derivative < Operator < Expr
@pyclass ufl Grad Expression
@pyclass ufl NablaGrad Expression
@pyclass ufl NablaDiv Expression

@pyfunc fenics grad
@pyfunc fenics nabla_grad
@pyfunc fenics dot
@pyfunc fenics inner
@pyfunc fenics sym
@pyfunc fenics lhs
@pyfunc fenics rhs

@pyfunc ufl nabla_div

OpType = Union{Expression, FeFunction}

tr(u::OpType) = Expression(fenics.tr(u.pyobject))
div(u::OpType) = Expression(fenics.div(u.pyobject))
sqrt(u::OpType) = Expression(fenics.sqrt(u.pyobject))

export Measure, Identity, nabla_div, grad, nabla_grad, dot, tr, inner, sym, lhs, rhs

############################################
#
#    mshr
#
# https://bitbucket.org/fenics-project/mshr/wiki/Home
#
############################################

# no need to export
@pyclass mshr CSGGeometry

@pyclass mshr Rectangle CSGGeometry
@pyclass mshr Circle CSGGeometry

@pyfunc mshr generate_mesh

export Rectangle, Circle, generate_mesh

############################################
#
#    dolfin.cpp.io module
#
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.io.html
#
############################################

@pyclass fenics File
@pyclass fenics XDMFFile

export File, XDMFFile

############################################
#
#    dolfin.cpp.adaptivity
#
# https://fenicsproject.org/docs/dolfin/2018.1.0/python/_autogenerated/dolfin.cpp.adaptivity.html?highlight=timeseries#dolfin.cpp.adaptivity.TimeSeries
#
############################################

@pyclass fenics TimeSeries
export TimeSeries

############################################
#    Operator
############################################

#import Base: +, -, *, /, ^, ==, <<   # must be explicitly imported to be extended

*(expr1::OpType , expr2::Measure) = Expression(expr2.pyobject.__rmul__(expr1.pyobject) )

*(expr1::OpType , expr2::OpType ) = Expression(expr1.pyobject.__mul__(expr2.pyobject) )
*(expr1::Real, expr2::OpType) = Expression(expr2.pyobject.__mul__(expr1) )
*(expr1::OpType , expr2::Real) = Expression(expr1.pyobject.__mul__(expr2) )

+(expr1::OpType, expr2::Real) = Expression(expr1.pyobject.__add__(expr2) )
+(expr1::Real, expr2::OpType) = Expression(expr2.pyobject.__add__(expr1) )
+(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject.__add__(expr2.pyobject) )

-(expr::OpType) = (-1)*expr
-(expr1::OpType, expr2::Real) = Expression(expr1.pyobject.__sub__(expr2) )
-(expr1::Real, expr2::OpType) = (-1)*(Expression(expr2.pyobject.__sub__(expr1) ))
-(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject.__sub__(expr2.pyobject) )

/(expr1::OpType, expr2::Real) = Expression(expr1.pyobject.__div__(expr2) )
/(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject.__div__(expr2.pyobject) )
/(expr1::Real,expr2::OpType) = (expr1 * expr2) / (expr2 * expr2)

^(expr1::OpType, expr2::Real) = Expression(expr1.pyobject.__pow__(expr2) )
^(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject.__pow__(expr2.pyobject) )


==(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject == expr2.pyobject)
==(expr1::OpType, expr2::Real) = Expression(expr1.pyobject == expr2)
==(expr1::Real, expr2::OpType) = Expression(expr2.pyobject == expr1)

<<(file::File, u::Mesh) =  file.pyobject << u.pyobject
<<(file::File, u::FeFunction) =  file.pyobject << u.pyobject
<<(file::File, u::Tuple{FeFunction,Real}) =  file.pyobject <<  (u[1].pyobject, u[2])

+(geo1::CSGGeometry, geo2::CSGGeometry) = CSGGeometry(geo1.pyobject + geo2.pyobject)
-(geo1::CSGGeometry, geo2::CSGGeometry) = CSGGeometry(geo1.pyobject - geo2.pyobject) 

#array(matrix::FeMatrix) = matrix.pyobject.gather_on_zero()
array(form::Expression) =  array(assemble(form)).gather_on_zero()
function array(solution::FeFunction)
    generic_vector = solution.pyobject.vector()
    instantiated_vector = fenics.Vector(generic_vector)
    instantiated_vector.gather_on_zero()
end
export array

len(U::OpType) = length(U.pyobject)
export len

@pyfunc fenics as_vector
export as_vector







