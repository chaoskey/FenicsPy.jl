@pyclass dolfin FeAny


############################################
#
#    dolfin.common
#    dolfin.cpp.parameter
#
# https://fenicsproject.org/docs/dolfin/2017.2.0/python/programmers-reference/common/index.html
#
############################################

# no need to export
@pyclass dolfin Parameters

# no need to export
# not exported in `dolfin` or `ufl`
@pyclass dolfin Parameter
@pyclass dolfin GlobalParameters

@pyfunc dolfin plot

export Parameters, plot




