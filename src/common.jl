@pyclass dolfin FeAny


const DEBUG = 10
const INFO = 20
const WARNING = 30
const ERROR = 40
const CRITICAL = 40
const PROGRESS = 16

export DEBUG, INFO, WARNING, ERROR, CRITICAL, PROGRESS


############################################
#
#    dolfin.cpp.common
#    dolfin.common
#    dolfin.cpp.parameter
#    dolfin.cpp.refinement
#
# https://fenicsproject.org/docs/dolfin/latest/python/_autogenerated/dolfin.cpp.common.html
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/common/index.html
# https://fenicsproject.org/docs/dolfin/latest/python/_autogenerated/dolfin.cpp.parameter.html
# https://fenicsproject.org/docs/dolfin/latest/python/_autogenerated/dolfin.cpp.refinement.html
#
############################################

@pyclass dolfin Variable

# no need to export
@pyclass dolfin Parameters

# no need to export
# not exported in `dolfin` or `ufl`
@pyclass dolfin Parameter
@pyclass dolfin GlobalParameters


@pyfunc dolfin has_mpi
@pyfunc dolfin has_debug
@pyfunc dolfin has_hdf5
@pyfunc dolfin has_hdf5_parallel
@pyfunc dolfin has_parmetis
@pyfunc dolfin has_petsc
@pyfunc dolfin has_scotch
@pyfunc dolfin has_slepc
@pyfunc dolfin has_sundials
@pyfunc dolfin refine

#@pyfunc dolfin plot
plot(u::FeObject, args...; kwargs...) = dolfin.plot(u.pyobject, args...; kwargs...)

export Variable, Parameters, plot, has_mpi, has_debug, has_hdf5, has_hdf5_parallel, has_parmetis,
       has_petsc, has_scotch, has_slepc, has_sundials, refine



############################################
#
#    dolfin.cpp.log
#
# https://fenicsproject.org/docs/dolfin/latest/python/_autogenerated/dolfin.cpp.log.html
#
############################################

@pyclass dolfin Progress

@pyfunc dolfin get_log_level
@pyfunc dolfin set_log_level
@pyfunc dolfin set_log_active
@pyfunc dolfin info
log(level::Integer, str::String) = dolfin.cpp.log.log(dolfin.LogLevel(level), str)


export Progress, get_log_level, set_log_level, set_log_active, info

