mshrfunc = []
mshrclass = []

fenicsfunc = []
fenicsclass = []

uflfunc = []
uflclass = []

export dx, ds, dS, dP

############################################
#    mshr
############################################

push!(mshrclass, [:Rectangle, :Circle])
push!(mshrfunc, [:generate_mesh])
export Rectangle, Circle, generate_mesh

############################################
#    dolfin.cpp
############################################

push!(fenicsclass, [:TimeSeries])
array(matrix) = matrix.gather_on_zero()
export TimeSeries, array

############################################
#    dolfin.cpp.io module
############################################

push!(fenicsclass, [:XDMFFile, :File])
export XDMFFile, File

############################################
#    dolfin.cpp.mesh module
#
# https://fenicsproject.org/docs/dolfin/1.6.0/python/programmers-reference/cpp/mesh/index.html
#
############################################

push!(fenicsclass, [:Mesh, :BoxMesh, :Point, :RectangleMesh, :UnitSquareMesh])
export Mesh, BoxMesh, Point, RectangleMesh, UnitSquareMesh

############################################
#    dolfin.cpp.fem module
#
# https://fenicsproject.org/docs/dolfin/1.6.0/python/programmers-reference/cpp/fem/index.html
#
############################################

push!(fenicsclass, :FiniteElement, :MixedElement)
push!(fenicsfunc, :assemble)
export FiniteElement, MixedElement, assemble

#############################################
#    dolfin.fem.module
#
# https://fenicsproject.org/docs/dolfin/1.6.0/python/programmers-reference/fem/index.html
#
############################################

push!(fenicsclass, [:DirichletBC, :AutoSubDomain])
push!(fenicsfunc, [:interpolate, :errornorm, :project, :solve])
export DirichletBC, AutoSubDomain, interpolate, errornorm, project, solve

############################################
#    dolfin.functions module
############################################

import Base: Function  # must be explicitly imported to be extended

push!(fenicsclass, [:Constant, :Expression])
push!(fenicsfunc, [:VectorFunctionSpace, :FunctionSpace, :TrialFunction, :TestFunction, :Function, :FacetNormal, :split])
export Constant, Expression, VectorFunctionSpace, FunctionSpace, TrialFunction, TestFunction, Function, FacetNormal, split

############################################
#    ufi      https://fenicsproject.org/pub/documents/ufl/ufl-user-manual/ufl-user-manual.pdf
############################################

import Base: div, sqrt  # must be explicitly imported to be extended

push!(fenicsclass, :Identity)
push!(fenicsfunc, [:div, :grad, :nabla_grad, :dot, :tr, :sqrt, :inner, :sym, :lhs, :rhs])
push!(uflfunc, :nabla_div)
export Identity, div, nabla_div, grad, nabla_grad, dot, tr, sqrt, inner, sym, lhs, rhs

fenicsclass = vcat(fenicsclass...)
fenicsfunc = vcat(fenicsfunc...)
mshrclass = vcat(mshrclass...)
mshrfunc = vcat(mshrfunc...)
uflclass = vcat(uflclass...)
uflfunc = vcat(uflfunc...)

