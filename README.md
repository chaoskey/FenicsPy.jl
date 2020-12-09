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

dolfin ：https://fenicsproject.org/docs/dolfin

ufl : https://fenicsproject.org/docs/ufl

mshr : https://bitbucket.org/fenics-project/mshr/wiki/Home

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

错误现象如图：

![](https://chaoskey.github.io/notes/001.jpg)

解决方法，直接修改`plotting.py`对应的代码，如图：

![](https://chaoskey.github.io/notes/003.jpg)

修改后执行如图：

![](https://chaoskey.github.io/notes/002.jpg)


## FenicsPy 和 FEniCS 的某些差异

我的目标：FenicsPy.jl完全和Python版本的FEniCS使用方法一样。但实际上比如有少许差异。

- 1） 在`FenicsPy.jl`的库导入，只需要`using FenicsPy`， 就同时支持包括`dolfin`、`ufl`、`mshr`、 `plot`。

- 2） `FEniCS`的`Function`, 由于命名冲突的问题，在`FenicsPy.jl`中改名为`FeFunction`。

- 3） 关于自定义表达式，即继承于`UserExpression`的实现：

由于待实现函数`eval_cell(self, values, x, cell)`中的参数values要求`共享传参`， 但在`PyCall.jl`的实现中没有做到这点，所以通过写`julia`代码继承`UserExpression`做到参数values`共享传参很麻烦。 好在`PyCall.jl`支持调用整段`Python`代码。所以可以写出如下代码：

```julia
 Define magnetic permeability
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


