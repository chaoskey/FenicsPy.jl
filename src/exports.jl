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
export TimeSeries

############################################
#    dolfin.cpp.io module
############################################

push!(fenicsclass, [:XDMFFile, :File])
export XDMFFile, File

############################################
#    dolfin.cpp.mesh module
############################################

push!(fenicsclass, :Point)
export Point

############################################
#    dolfin.cpp.fem module
############################################

push!(fenicsfunc, :assemble)
export assemble

#############################################
#    dolfin.fem.module
############################################

push!(fenicsclass, [:DirichletBC, :AutoSubDomain])
push!(fenicsfunc, [:solve])
export DirichletBC, AutoSubDomain, solve

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

