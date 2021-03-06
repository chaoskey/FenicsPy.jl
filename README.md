# FenicsPy.jl

FenicsPy.jl是FEniCS的Julia封装，有限元法求解PDE。 

缘由：FEniCS.jl不好用， 于是重新写了FEniCS的Julia封装：FenicsPy.jl。

目标：完全和Python版本的FEniCS使用方法一样。

## 安装

1.  安装FEniCS:     https://fenicsproject.org/download/

2.  安装FenicsPy.jl

```julia
import Pkg; Pkg.add("https://github.com/chaoskey/FenicsPy.jl")
```

或

```julia
]
pkg > add https://github.com/chaoskey/FenicsPy.jl
```

3. 测试（目前就是把所有例子执行一遍，时间比较长，慎重执行）

```julia
]
pkg > test FenicsPy
```

更详细的安装说明，不妨参考： [科学计算环境搭建（Win10+WSL2+Ubuntu）](https://chaoskey.gitee.io/notes/docs/julia/0095/)

## 文档

由于使用方法完全和Python版本的FEniCS一样，所以完全可以参考FEniCS官方文档。

FEniCS : https://fenicsproject.org/documentation/

dolfin ：https://fenicsproject.org/olddocs/dolfin

ufl : https://fenicsproject.org/olddocs/ufl

mshr : https://bitbucket.org/fenics-project/mshr/wiki/Home

The FEniCS Book(翻译中) ： https://chaoskey.gitee.io/notes/docs/fem/

## 范例

[有限元法解偏微分方程](https://chaoskey.gitee.io/notes/docs/julia/0094/)

[有限元法求解牛顿流体](https://chaoskey.gitee.io/notes/docs/julia/0096/)

[FenicsPy.jl范例](https://gitee.com/chaoskey/FenicsPy.jl/tree/master/examples)

[FEniCS官方教程](https://fenicsproject.org/tutorial/) 

----------------------------------------

## FEniCS自身问题的解决

`FeincsPy.jl`仅是`FEniCS`的Julia封装，所以`FeincsPy.jl`无法解决`FEniCS`自身问题。

下面的问题基于的`FEniCS`的版本是：

> fenics=2019.1.0=py37_9
> mshr=2019.1.0=py37hf9f41d3_3


- 1） Axes3D currently only supports the aspect argument 'auto'. You passed in 'equal'.

错误复现代码：

```julia
using FenicsPy

mesh = UnitCubeMesh(10, 10, 10)
plot(mesh, title="Unit cube")
```

错误解决方法: 直接修改`plotting.py`对应的代码，根据错误提示的行号，将`ax.set_aspect('equal')` 缩进一个`Tab`或`4个空格`即可。

## FenicsPy 和 FEniCS 的某些差异

我的目标：FenicsPy.jl完全和Python版本的FEniCS使用方法一样。但实际上有少许差异。

### 1）库导入

在`FenicsPy.jl`的库导入，只需要`using FenicsPy`， 就同时支持包括`dolfin`、`ufl`、`mshr`、 `plot`。

### 2）关于`Function`的命名冲突 

`FEniCS`的`Function`, 由于命名冲突的问题，在`FenicsPy.jl`中改名为`FeFunction`。

### 3）关于自定义表达式，即继承于`UserExpression`的实现：

由于待实现函数`eval_cell(self, values, x, cell)`中的参数values要求`共享传参`， 但在`PyCall.jl`的实现中没有做到这点，所以通过写`julia`代码继承`UserExpression`做到values`共享传参`很麻烦。 好在`PyCall.jl`支持调用整段`Python`代码。所以可以写出如下代码：

```julia
# Define magnetic permeability
py"""
from math import pi
from dolfin import UserExpression
class Permeability(UserExpression):
    def __init__(self, markers, **kwargs):
        self.markers = markers
        super().__init__(**kwargs)
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 0:
            values[0] = 4*pi*1e-7 # vacuum
        elif self.markers[cell.index] == 1:
            values[0] = 1e-5      # iron (should really be 6.3e-3)
        else:
            values[0] = 1.26e-6   # copper
    def value_shape(self):
        return ()
"""

μ = Expression(py"Permeability"(markers.pyobject, degree=1))
```

完整的代码，参考[FenicsPy.jl/examples/ft11_magnetostatics.jl](https://gitee.com/chaoskey/FenicsPy.jl/blob/master/examples/ft11_magnetostatics.jl)

### 4）继承Python类的第二种方法

形如：

```julia
# Define Dirichlet boundary
@pydef mutable struct DirichletBoundary <: dolfin.SubDomain
    function inside(self, x, on_boundary)
        return on_boundary
    end
end

bc = SubDomain(DirichletBoundary())
```

完整的代码，参考[FenicsPy.jl/examples/demo_biharmonic.jl](https://gitee.com/chaoskey/FenicsPy.jl/blob/master/examples/demo_biharmonic.jl)


此第二种方法只适用于：成员函数的参数不涉及`共享传参`的情况，这是因为`PyCall.jl` 不支持跨`共享传参`。 

我更乐于统一采用第一种方法（即3）对应的方法）。

### 5）关于边界条件的设置的差异

`FEniCS`的范例代码：

```python
# Sub domain for clamp at left end
def left(x, on_boundary):
    return near(x[0], 0.) and on_boundary

# Sub domain for rotation at right end
def right(x, on_boundary):
    return near(x[0], 1.) and on_boundary

bc = DirichletBC(V, zero, left)

force_boundary = AutoSubDomain(right)
```

但是，我没有找到合适的方法，可以在`FenicsPy.jl`中仿照这个代码。  目前建议采用下面的代码实现相同的功能：


```julia

# Sub domain for clamp at left end
left = "near(x[0], 0.) && on_boundary"

# Sub domain for rotation at right end
right = "near(x[0], 1.) && on_boundary"

bc = DirichletBC(V, zero, left)

force_boundary = CompiledSubDomain(right)
```
完整的代码，参考：[examples/demo_elastodynamics.jl](https://gitee.com/chaoskey/FenicsPy.jl/blob/master/examples/demo_elastodynamics.jl)


### 6）关于绘图

`FenicsPy.jl`关于绘图的接口函数只有`plot(...)`。这个接口函数的第一个参数类型只支持：`Mesh`, `MeshFunction`, `FeFunction`, `Expression`, `DirichletBC`, `FiniteElement`,`MultiMesh`， 所以与`PyPlot.jl`的原本的`plot`不冲突。

事实上， `plot`后台调用的就是`PyPlot.plot`。 目前也只支持`PyPlot` 

如果只使用`plot`，无论绘制上面的特定`FEniCS`数据，还是绘制普通数据，都无需额外`using PyPlot`。 

但是，如果希望用到`PyPlot`的其它API，那么必须额外`using PyPlot`。

比如：

```julia
using PyPlot
plot(time, u_tip)
xlabel("Time")
ylabel("Tip displacement")
ylim(-0.5, 0.5)
```

完整的代码，参考[examples/demo_elastodynamics.jl](https://gitee.com/chaoskey/FenicsPy.jl/blob/master/examples/demo_elastodynamics.jl) 末尾注释掉的代码。



