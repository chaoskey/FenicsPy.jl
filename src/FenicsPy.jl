module FenicsPy

using PyCall

# must be explicitly imported to be extended
import Base: getproperty, split, inv, transpose, div, diff, abs, sign, sqrt, exp,  cos, sin, tan, acos, asin, atan,  cosh, sinh, tanh, *, +, -, /, ^, ==, <<
import PyPlot: plot

##################################
#    Base : FeObject
##################################

abstract type
  FeObject
end

##################################
#    input
##################################

FeType = Union{Real, String, Array, Tuple, Dict, FeObject}

to_pyarray(arr::Array) = [isa(a, FeObject) ? a.pyobject : a for a in arr]

to_pytuple(arr::Tuple) = Tuple(isa(a, FeObject) ? a.pyobject : a for a in arr)

to_pydict(arr::Dict) = Dict(i => (isa(arr[i], FeObject) ? arr[i].pyobject : arr[i]) for i in eachindex(arr))

function to_pytype(obj::FeType)
    if isa(obj, FeObject)
        obj.pyobject
    elseif isa(obj, Array)
        to_pyarray(obj)
    elseif isa(obj, Tuple)
        to_pytuple(obj)
    elseif isa(obj, Dict)
        to_pydict(obj)
    else
        obj
    end
end

function args_conv(args...; kwargs...)
    _args = [to_pytype(arg) for arg in args]
    _kwargs = Dict(idx =>  to_pytype(kwargs[idx]) for idx in eachindex(kwargs))
    _args, _kwargs
end

##################################
#    output
##################################

PyType = Union{Real, String, Array, Tuple, Dict, PyObject, Nothing}

function to_feobject(obj::PyObject)
    clsname = match(r"(\w+)'>", string(obj.__class__))[1]
    if haskey(class_dict,clsname)
        clsname = class_dict[clsname]
        getproperty(FenicsPy, clsname)(obj)
    else
        obj
    end
end

to_fearray(arr::Array) = [isa(a, PyObject) ? to_feobject(a) : a for a in arr]

to_fetuple(arr::Tuple) = Tuple(isa(a, PyObject) ? to_feobject(a) : a for a in arr)

to_fedict(arr::Dict) = Dict(i => (isa(arr[i], PyObject) ? to_feobject(arr[i]) : arr[i]) for i in eachindex(arr))

function to_fetype(obj::PyType)
    if isa(obj, PyObject)
        to_feobject(obj)
    elseif isa(obj, Array)
        to_fearray(obj)
    elseif isa(obj, Tuple)
        to_fetuple(obj)
    elseif isa(obj, Dict)
        to_fedict(obj)
    else
        obj
    end
end

##################################
#    @pyfunc
#
# examples :
#   1  @pyfunc fenics TestFunction
#   2  @pyfunc mshr generate_mesh  you_alias
##################################
macro pyfunc(_module::Symbol, name::Symbol, alias::Symbol=:nothing)
    alias = alias === :nothing ? name : alias
    name = string(name)
    esc(quote
        $(alias)(args::FeType...; kwargs...) = begin
            _args, _kwargs = args_conv(args...; kwargs...)
            obj = getproperty($_module, $name)(_args...; _kwargs...)
	    to_fetype(obj)
        end
    end)
end

export pyfunc

##################################
#    @pyclass
#
# examples :
#    1  @pyclass fenics Point
#    2  @pyclass mshr Rectangle CSGGeometry
#    3  @pyclass fenics Function FeObject FeFunction
##################################
const class_dict = Dict()
macro pyclass(_module::Symbol, name::Symbol, _base::Symbol=:FeObject, alias::Symbol=:Nothing)
    alias = alias === :Nothing ? name : alias
    impl = Symbol(alias, "Impl")
    name = string(name)
    class_dict[name] = alias
    esc(quote
        abstract type
            $alias <: $_base
        end

        struct $impl <: $alias
            pyobject::PyObject
        end

        ###############################
        #    Constructors ( copy )
        #
        # examples :
        #     dx = Measure(fenics.dx)
        ###############################
        $(alias)(pyobject::PyObject) = $impl(pyobject)

        ###############################
        #    Constructors ( parameter )
        #
        # examples :
        #     p = Point(1.2, 2.3)
        ###############################
        $(alias)(args::FeType...; kwargs...) = begin
            _args, _kwargs = args_conv(args...; kwargs...)
            name = hasproperty($_module, $name) ? $name : $_base
            _obj = getproperty($_module, name)(_args...; _kwargs...)
            $impl(_obj)
        end

        ###############################
        # examples :
        #     p.x()    # = 1.2
        #     p.pyobject
        ###############################
        function getproperty(_obj::$impl, _sym::Symbol)
            if _sym === :pyobject
	        return getfield(_obj, _sym)
            end
            o = getproperty(_obj.pyobject, _sym)
            if match(r"(\w+)'>", string(o.__class__))[1] == "method"
                return (args::FeType...; kwargs...) ->  begin
                       _args, _kwargs = args_conv(args...; kwargs...)
		       to_fetype(o(_args...; _kwargs...))
                end
            end
            to_fetype(o)
        end

        ###############################
        # https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects
        #
        # examples :
        #    dx = Measure("dx", domain=mesh, subdomain_data=markers)
        #    dx(1)
        ###############################
        function (o::$impl)(args...; kwargs...)
            _args, _kwargs = args_conv(args...; kwargs...)
            to_fetype(o.pyobject(_args...; _kwargs...))
	end
        
    end)
end

export pyclass

##################################
#    python module : fenics,  mshr, ufl
##################################

include("exports.jl")

##################################
#    init
##################################

const fenics = PyCall.PyNULL()
const ufl = PyCall.PyNULL()
const mshr = PyCall.PyNULL()

export fenics, ufl, mshr

function __init__()

    copy!(fenics, pyimport_conda("fenics", "fenics", "conda-forge"))
    copy!(mshr, pyimport_conda("mshr", "mshr", "conda-forge"))
    copy!(ufl, pyimport_conda("ufl", "ufl", "conda-forge"))
    
    global dx = Measure(fenics.dx)
    global ds = Measure(fenics.ds)
    global dS = Measure(fenics.dS)
    global dP = Measure(fenics.dP)
    
    global tetrahedron = Cell(fenics.tetrahedron)
    global hexahedron = Cell(fenics.hexahedron) #matplotlib cannot handle hexahedron elements
    global triangle = Cell(fenics.triangle)
    global quadrilateral = Cell(fenics.quadrilateral)
						    
end

export dx, ds, dS, dP,
       tetrahedron, hexahedron, triangle, quadrilateral

end # FenicsPy
