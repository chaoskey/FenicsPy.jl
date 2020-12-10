
############################################
#
#    ufi Operator
#
# https://fenicsproject.org/docs/ufl/1.6.0/ufl.html
# < Operator < Expr 
############################################

# no need to export
# not exported in `feincs` and `ufl` 
@pyclass ufl ListTensor Expression
@pyclass ufl ComponentTensor Expression
@pyclass ufl Indexed Expression

# export
@pyfunc ufl variable

export variable

############################################
#
#    ufi : Tensor algebra operators
#
# https://fenicsproject.org/docs/ufl/1.6.0/ufl.html
# < CompoundTensorOperator < Operator < Expr 
############################################

# no need to export
# not exported in `dolfin` or `ufl`
@pyclass ufl Outer Expression
@pyclass ufl Inner Expression
@pyclass ufl Conj Expression
@pyclass ufl Dot Expression
@pyclass ufl Cross Expression
@pyclass ufl Determinant Expression
@pyclass ufl Inverse Expression
@pyclass ufl Cofactor Expression
@pyclass ufl Transposed Expression
@pyclass ufl Trace Expression
@pyclass ufl Deviatoric Expression
@pyclass ufl Skew Expression
@pyclass ufl Sym Expression
@pyclass ufl Product Expression
@pyclass ufl Sum Expression
@pyclass ufl Division Expression
@pyclass ufl PositiveRestricted Expression
@pyclass ufl NegativeRestricted Expression

# export
@pyfunc ufl outer
@pyfunc ufl inner
@pyfunc ufl dot
@pyfunc ufl cross
@pyfunc ufl perp
@pyfunc dolfin det
inv(u::OpType) = Inverse(ufl.inv(u.pyobject)) # Base.inv
@pyfunc ufl cofac
transpose(u::OpType) = Transposed(ufl.transpose(u.pyobject)) # Base.transpose
@pyfunc ufl tr
@pyfunc ufl diag
@pyfunc ufl diag_vector
@pyfunc ufl dev
@pyfunc ufl skew
@pyfunc ufl sym
@pyfunc ufl avg
@pyfunc ufl jump

export outer, inner, dot, cross, perp, det, cofac, tr, diag, diag_vector, dev, skew, sym, avg, jump

############################################
# Differential operators
#
# https://fenicsproject.org/docs/ufl/1.6.0/ufl.html
# < CompoundDerivative < Derivative < Operator < Expr
############################################

# no need to export
# not exported in `dolfin` or `ufl`
@pyclass ufl VariableDerivative Expression
@pyclass ufl Grad Expression
@pyclass ufl Div Expression
@pyclass ufl NablaGrad Expression
@pyclass ufl NablaDiv Expression
@pyclass ufl Curl Expression

# export
diff(f::OpType,v::Variable) = VariableDerivative(ufl.diff(f.pyobject, v.pyobject))  # Base.diff
@pyfunc dolfin derivative
@pyfunc ufl grad
div(u::OpType) = Div(ufl.div(u.pyobject)) # Base.div
@pyfunc ufl nabla_grad
@pyfunc ufl nabla_div
@pyfunc ufl Dx
@pyfunc ufl Dn
@pyfunc ufl curl
@pyfunc ufl rot

export diff, derivative, grad, nabla_div, nabla_grad, Dx, Dn, curl, rot

############################################
# Nonlinear functions
#
# https://fenicsproject.org/docs/ufl/1.6.0/ufl.html
############################################

# no need to export
# not exported in `dolfin` or `ufl`
@pyclass ufl Ln Expression 
@pyclass ufl Erf Expression

@pyfunc ufl max_value
@pyfunc ufl min_value
abs(u::OpType) = Expression(ufl.abs(u.pyobject) * u.pyobject ) # Base.abs
sign(u::OpType) = Expression(ufl.sign(u.pyobject)) # Base.sign
sqrt(u::OpType) = Expression(ufl.sqrt(u.pyobject)) # Base.sqrt
exp(u::OpType) = Expression(ufl.exp(u.pyobject)) # Base.exp
@pyfunc ufl ln
@pyfunc ufl erf
cos(u::OpType) = Expression(ufl.cos(u.pyobject)) # Base.cos
sin(u::OpType) = Expression(ufl.sin(u.pyobject)) # Base.sin
tan(u::OpType) = Expression(ufl.tan(u.pyobject)) # Base.tan
acos(u::OpType) = Expression(ufl.acos(u.pyobject)) # Base.acos
asin(u::OpType) = Expression(ufl.asin(u.pyobject)) # Base.asin
atan(u::OpType) = Expression(ufl.atan(u.pyobject)) # Base.atan
@pyfunc ufl atan_2
cosh(u::OpType) = Expression(ufl.cosh(u.pyobject)) # Base.cosh
sinh(u::OpType) = Expression(ufl.sinh(u.pyobject)) # Base.sinh
tanh(u::OpType) = Expression(ufl.tanh(u.pyobject)) # Base.tanh
#@pyfunc dolfin bessel_J
#@pyfunc dolfin bessel_Y
#@pyfunc dolfin bessel_I
#@pyfunc dolfin bessel_K

export max_value, min_value, ln, erf, atan_2, bessel_J, bessel_Y, bessel_I, bessel_K

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
/(expr1::Real,expr2::OpType) = Constant(expr1) / expr2

^(expr1::OpType, expr2::Real) = Expression(expr1.pyobject.__pow__(expr2) )
^(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject.__pow__(expr2.pyobject) )


==(expr1::OpType, expr2::OpType) = Expression(expr1.pyobject == expr2.pyobject)
==(expr1::OpType, expr2::Real) = Expression(expr1.pyobject == expr2)
==(expr1::Real, expr2::OpType) = Expression(expr2.pyobject == expr1)

# Base: split
function split(fun::OpType)
     vec = dolfin.split(fun.pyobject)
     expr_vec = [to_fetype(spl) for spl in vec]
     return expr_vec
end

