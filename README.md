# FenicsPy.jl

FenicsPy.jl是FEniCS的Julia封装，有限元法求解PDE。 

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




