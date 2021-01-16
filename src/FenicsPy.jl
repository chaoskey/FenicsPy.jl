module FenicsPy

using PyCall

export @pydef, @py_str

# must be explicitly imported to be extended
import Base: split, inv, transpose, div, diff, abs, sign, sqrt, exp,  cos, sin, tan, acos, asin, atan,  cosh, sinh, tanh, log, replace, adjoint, *, +, -, /, ^, ==, <<

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

FeType = Union{Real, String, Array, Tuple, Dict, Function, FeObject}

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

PyType = Union{Real, String, Array, Tuple, Dict, Function, PyObject, Nothing}

function to_feobject(obj::PyObject)
    clsname = match(r"(\w+)'>", string(obj.__class__))[1]
    clsname = get(class_dict, clsname, :FeAny)
    getproperty(FenicsPy, clsname)(obj)
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
#   1  @pyfunc dolfin TestFunction
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
#    1  @pyclass dolfin Point
#    2  @pyclass mshr Rectangle CSGGeometry
#    3  @pyclass dolfin Function FeObject FeFunction
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
        #     dx = Measure(dolfin.dx)
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
        function Base.getproperty(_obj::$impl, _sym::Symbol)
            if _sym === :pyobject
	        return getfield(_obj, _sym)
            end
            o = getproperty(_obj.pyobject, _sym)
            if isa(o, PyObject) && match(r"(\w+)'>", string(o.__class__))[1] == "method"
                return (args::FeType...; kwargs...) ->  begin
                       _args, _kwargs = args_conv(args...; kwargs...)
		       to_fetype(o(_args...; _kwargs...))
                end
            end
            to_fetype(o)
        end

        function Base.setproperty!(_obj::$impl, _sym::Symbol, value)
            if _sym === :pyobject
                return
            end
            setproperty!(_obj.pyobject, _sym, to_pytype(value))
        end

        Base.getindex(o::$impl, key::Union{String, Symbol}) = to_fetype(get(o.pyobject, key))
        Base.getindex(o::$impl, idx::Int) = to_fetype(get(o.pyobject, idx-1))
        Base.setindex!(o::$impl, value, key::Union{String, Symbol}) = set!(o.pyobject, key, to_pytype(value))
        Base.setindex!(o::$impl, value, idx::Int) = set!(o.pyobject, idx-1, to_pytype(value))

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


        ###############################
        # examples :
        #     CellType.Type.quadrilateral
        #     UnitSquareMesh.create
        ###############################
        function Base.getproperty(_obj::Type{$alias}, _sym::Symbol)
            if hasproperty(_obj, _sym)
                return getfield(_obj, _sym)
            end
            name = hasproperty($_module, $name) ? $name : $_base
            pyclass = getproperty($_module, name)
            o = getproperty(pyclass, _sym)
            if isa(o, PyObject) && match(r"(\w+)'>", string(o.__class__))[1] == "method"
                return (args::FeType...; kwargs...) ->  begin
                       _args, _kwargs = args_conv(args...; kwargs...)
                       to_fetype(o(_args...; _kwargs...))
                end
            end
            to_fetype(o)
        end
        
    end)
end

export pyclass

##################################
#    python module : dolfin,  mshr, ufl
##################################

include("common.jl")
include("la.jl")
include("function.jl")
include("fem.jl")
include("mesh.jl")
include("geometry.jl")
include("operator.jl")
include("mshr.jl")
include("io.jl")

##################################
#    init
##################################

const dolfin = PyCall.PyNULL()
const ufl = PyCall.PyNULL()
const mshr = PyCall.PyNULL()

export dolfin, ufl, mshr

function __init__()

    copy!(dolfin, pyimport_conda("dolfin", "fenics-dolfin", "conda-forge"))
    copy!(mshr, pyimport_conda("mshr", "mshr", "conda-forge"))
    copy!(ufl, pyimport_conda("ufl", "fenics-ufl", "conda-forge"))
    
    global dx = Measure(dolfin.dx)
    global ds = Measure(dolfin.ds)
    global dS = Measure(dolfin.dS)
    global dP = Measure(dolfin.dP)
    
    global tetrahedron = Cell(dolfin.tetrahedron)
    global hexahedron = Cell(dolfin.hexahedron) #matplotlib cannot handle hexahedron elements
    global triangle = Cell(dolfin.triangle)
    global quadrilateral = Cell(dolfin.quadrilateral)

    # parameters["linear_algebra_backend"] = backendname  # list_linear_algebra_backends()
    global parameters = GlobalParameters(dolfin.parameters)

    global DOLFIN_EPS = dolfin.DOLFIN_EPS 
end

export dx, ds, dS, dP,
       tetrahedron, hexahedron, triangle, quadrilateral,
       parameters,  
       DOLFIN_EPS

end # FenicsPy
