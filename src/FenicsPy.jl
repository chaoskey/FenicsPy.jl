module FenicsPy

using PyCall

function pyfunc(obj, name::Symbol, alias::Symbol=:nothing)
    alias = alias == :nothing ? name : alias
    func = getproperty(obj, name)
    @eval $alias(args...; kwargs...) = $func(args...; kwargs...)
end

function pyclass(obj, name::Symbol, alias::Symbol=:nothing)
    alias = alias == :nothing ? name : alias
    bas = getproperty(obj, name)
    @eval @pydef mutable struct $alias <: $bas
        __init__(self, args...; kwargs...) = begin
            $bas.__init__(self, args...; kwargs...)
        end
    end
end

# python module : fenics,  mshr 
include("exports.jl")

const fenics = PyCall.PyNULL()
const ufl = PyCall.PyNULL()
const mshr = PyCall.PyNULL()

function __init__()

    try
        copy!(fenics, pyimport_conda("fenics", "fenics", "conda-forge"))

	global dx = fenics.dx
	global ds = fenics.ds
	global dS = fenics.dS
	global dP = fenics.dP

	[pyclass(fenics, vcat([cls]...)...) for cls in fenicsclass]
	[pyfunc(fenics, vcat([fun]...)...) for fun in fenicsfunc]

	println("Successfully pyimport fenics")
    catch ex
        println(ex)
    end

    try
        copy!(mshr, pyimport_conda("mshr", "mshr", "conda-forge"))

	[pyclass(mshr, vcat([cls]...)...) for cls in mshrclass]
	[pyfunc(mshr, vcat([fun]...)...) for fun in mshrfunc]

	println("Successfully pyimport mshr")
    catch ex
        println(ex)
    end

    try
        copy!(ufl, pyimport_conda("ufl", "ufl", "conda-forge"))

	[pyclass(ufl, vcat([cls]...)...) for cls in uflclass]
	[pyfunc(ufl, vcat([fun]...)...) for fun in uflfunc]

        println("Successfully pyimport ufl")
    catch ex
        println(ex)
    end
						    
end

end # FenicsPy
