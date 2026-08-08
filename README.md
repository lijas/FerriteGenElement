# FerriteGenElement.jl

This is a package for generating shape values for various elements and to cmpute the mapping from reference cell to the physical cell, for certain "exotic/classic" elements like Argyris and Morley. The package is in large based on the paper [Kirby, R. C. (2017). A general approach to transforming finite elements](https://doi.org/10.48550/arXiv.1706.09017).

# Introduction
We will now show how we can generate shape function for some elements. We will start of simple and generate the shape values a simple Linear triangle (i.e. Ferrite.Lagrange{RefTriangle,1}). We need to define three things. The main struct for the element, what kind of dofs the element has (and where they are placed), and the function space that the shape function span:

```
struct MyLinearTriangle <: Ferrite.Interpolation{RefTriangle,1}
end

FerriteGenElement.interpolation_dofs(::MyLinearTriangle) = (
    PointEvaluation(Vec((0.0, 0.0))),
    PointEvaluation(Vec((1.0, 0.0))),
    PointEvaluation(Vec((0.0, 1.0))),
)
FerriteGenElement.monomials(::MyLinearTriangle) = (
    Monomial((0,0)), # 1
    Monomial((1,0)), # x
    Monomial((0,1)), # y
)

```

Here, PointEvaluation( Vec((0.0, 0.0)) ) means that the the DoF is the "value" at the coordinate (0.0, 0.0). The linear Lagrange element has three value-dofs in each corner (which follows the definition of Ferrite's RefTriangle). The space spanned by the triangle is defined by the three monimals.

This is all that is needed to generate shape functions for this referenc element. We simply call `compute_reference_shape_values`

```
shape_values = compute_reference_shape_values(MyLinearTriangle())
```

If we want to see the shape functions, we can print them:

```
print_reference_shape_values(shape_values)
```

In this case, we get the familiar shape values:

```
function Ferrite.reference_shape_value(ip::MyLinearTriangle, ξ::Vec{2,T}, i::Int) where T
    x, y = ξ
    i==1 && return 1 - x - y
    i==2 && return x
    i==3 && return y
    throw(ArgumentError("    no shape function $i for interpolation $ip"))
end
```

The way the shape functions are computed is... TODO: Explain math how they are computed.

# The morley element

The same approach can be used generate the shape values for more special types of elements. We will now take a look at the classical [Morley](https://defelement.org/elements/examples/triangle-morley-2.html) element. This element is a non-conforming element with 6 DoF: Three value dofs in each cornder, and 3 DoFs ot the center of edge which represent the normal gradient at that point. This element was developed to study plate bending problems.

```
struct Morley{RefTriangle,order} <: Ferrite.Interpolation{RefTriangle,order}
end

FerriteGenElement.interpolation_dofs(::Morley{RefTriangle,2}) = (
    PointEvaluation(Vec((0.0, 0.0))),
    PointEvaluation(Vec((1.0, 0.0))),
    PointEvaluation(Vec((0.0, 1.0))), 
    NormalGradientEvaluation(Vec(1/2, 1/2), Vec(1/√2, 1/√2)),
    NormalGradientEvaluation(Vec(0.0, 1/2), Vec(1.0, 0.0)),
    NormalGradientEvaluation(Vec(1/2, 0.0), Vec((0.0,-1.0))),
)

FerriteGenElement.monomials(::Morley{RefTriangle,2}) = (
    Monomial((0,0)), # 1
    Monomial((1,0)), # x^1
    Monomial((2,0)), # x^2
    Monomial((0,1)), # y
    Monomial((1,1)), # x*y
    Monomial((0,2)), # y^2
)
```

Other the the standard PointEvaluation in each corner, we also define NormalGradientEvaluation(), which takes the coordinate (center of the edge) and the edge normal. Just as before, we can construct and print the shape functions with:

```
shape_values = compute_reference_shape_values(Morley{RefTriangle,2}())
print_reference_shape_values(shape_values)
```

which prints:

```
function Ferrite.reference_shape_value(ip::Morley{RefTriangle, 2}, ξ::Vec{2,T}, i::Int) where T
    x, y = ξ
    i==1 && return 1 - x - y + 2*x*y
    i==2 && return T(0.5)*x + T(0.5)*x^2 + T(0.5)*y - x*y - T(0.5)*y^2
    i==3 && return T(0.5)*x - T(0.5)*x^2 + T(0.5)*y - x*y + T(0.5)*y^2
    i==4 && return -T(0.7071067811865476)*x + T(0.7071067811865476)*x^2 - T(0.7071067811865476)*y + T(1.4142135623730951)*x*y + T(0.7071067811865476)*y^2
    i==5 && return x - x^2
    i==6 && return -y + y^2
    throw(ArgumentError("    no shape function $i for interpolation $ip"))
end
```

The result can be compared agains the shape values on [DefElement](https://defelement.org/elements/examples/triangle-morley-2.html), but note that Ferrite.RefTriangle follows a different vertex and edge ordering than DefElement.

Unfortinotly, these shape values can (currently) not just be implemented and used in Ferrite with standard element mapping tenchnology. The reason for this is becuase the NormalGradientEvaluation() makes the shape function couple when mapping the shape function from reference to the physical element. Kirby et al. show taht there is a relation between the shape function on the physical element and the reference element, via a matrix M

$$
N_{x} = M N_{\xi} 
$$
 
where $N_{\xi}$ is a vector of the shape function on the reference element and N_{x} is the shape functions on the physical element. For more details on how the matrix $M$ is computed, we referer to paper by Kirby. In this package, however, we can compute it like this:

```
# First create a description of our element (A linear triangle with some random coordinates).
mapping = Mapping(
    Lagrange{RefTriangle,1}(), 
    [Vec((0.01, 1.1)), Vec((1.0, 0.0)), Vec((0.1, 0.1))]
)

# Create the M matrix
M = FerriteGenElement.build_M_matrix(Morley{RefTriangle,2}(), mapping)
```

This generates a 6x6 Matrix which maps the shape values:

```
6×6 adjoint(::Matrix{Float64}) with eltype Float64:
 1.0  0.0  0.0   0.0           0.219512  -0.219512
 0.0  1.0  0.0  -4.44089e-17   0.0        0.219512
 0.0  0.0  1.0   4.44089e-17  -0.219512   0.0
 0.0  0.0  0.0   0.8           0.0        0.0
 0.0  0.0  0.0   0.0           0.883452   0.0
 0.0  0.0  0.0   0.0           0.0        0.883452
```

We can note that the normal gradient dofs on the edges couple with the value dofs in the corners. Note that this only maps the shape functions for this particular definition of the (physical) triangle with the specified vertecies. In practice, we would derive this matrix analtically (or symbolically), which we would then implement in Ferrite such that it works for any triangle in the mesh.

# Stuff to implement

The package is still under development, and there are some things that would be nice to have/improve.

- Some framework for handling Macro elements
- Integral moment DoFs



