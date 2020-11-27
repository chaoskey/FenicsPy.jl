mshrfunc = []
mshrclass = []

fenicsfunc = []
fenicsclass = []

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

push!(fenicsclass, [:Point, :UnitSquareMesh])
export Point, UnitSquareMesh

############################################
#    dolfin.cpp.fem module
############################################

push!(fenicsfunc, :assemble)
export assemble

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
push!(fenicsfunc, [:VectorFunctionSpace, :FunctionSpace, :TrialFunction, :TestFunction, :Function, :FacetNormal])
export Constant, Expression, VectorFunctionSpace, FunctionSpace, TrialFunction, TestFunction, Function, FacetNormal

############################################
#    ufi      https://fenicsproject.org/pub/documents/ufl/ufl-user-manual/ufl-user-manual.pdf
############################################

import Base: div  # must be explicitly imported to be extended

push!(fenicsclass, [:Identity])
push!(fenicsfunc, [:div, :grad, :nabla_grad, :dot, :inner, :sym, :lhs, :rhs])
export Identity, div, grad, nabla_grad, dot, inner, sym, lhs, rhs

fenicsclass = vcat(fenicsclass...)
fenicsfunc = vcat(fenicsfunc...)
mshrclass = vcat(mshrclass...)
mshrfunc = vcat(mshrfunc...)

