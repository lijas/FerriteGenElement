# FerriteGenElement.jl

`FerriteGenElement.jl` is a Julia package for generating shape functions for finite elements and, in particular, for handling the mapping from a reference cell to a physical cell for elements whose degrees of freedom (DoFs) are more general than standard "value" DoFs.

The package is particularly useful for "exotic" or classical finite elements such as the **Morley** and **Argyris** elements, for which the usual finite-element mapping machinery is not sufficient. The implementation is largely based on:

> Kirby, R. C. (2017). *A general approach to transforming finite elements*.  
> [https://doi.org/10.48550/arXiv.1706.09017](https://doi.org/10.48550/arXiv.1706.09017)

The central idea is to describe a finite element in terms of:

1. its interpolation type,
2. its degrees of freedom, represented as linear functionals, and
3. the polynomial (or other finite-dimensional) function space spanned by the shape functions.

From this information, `FerriteGenElement.jl` can automatically construct the shape functions on the reference cell. It can also construct the transformation needed to map these shape functions to a physical cell.

---

# Introduction

We will first demonstrate how shape functions can be generated for a simple element: the linear triangle, corresponding to `Ferrite.Lagrange{RefTriangle,1}`.

To define an element, we need to specify three things:

1. the main interpolation type,
2. the degrees of freedom of the element and where they are located, and
3. the function space spanned by the shape functions.

For example:

```julia
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

Here, `PointEvaluation(Vec((0.0, 0.0)))` means that the corresponding DoF is the value of the function evaluated at the point `(0.0, 0.0)`.

The linear Lagrange triangle therefore has three value DoFs, one at each vertex. The vertices follow the definition and ordering used by `Ferrite.RefTriangle`.

The function space is defined by the three monomials

$$
\mathcal{P}_1 = \{1,x,y\}.
$$

These are all the ingredients needed to generate the shape functions for this reference element.

We simply call:

```julia
shape_values = compute_reference_shape_values(MyLinearTriangle())
```

If we want to inspect the generated shape functions, we can use:

```julia
print_reference_shape_values(shape_values)
```

This produces:

```julia
function Ferrite.reference_shape_value(ip::MyLinearTriangle, ξ::Vec{2,T}, i::Int) where T
    x, y = ξ
    i==1 && return 1 - x - y
    i==2 && return x
    i==3 && return y
    throw(ArgumentError("    no shape function $i for interpolation $ip"))
end
```

As expected, these are the familiar linear triangular shape functions.

## How are the shape functions generated?

The shape functions are constructed directly from the definition of the finite element.

Suppose that the finite-dimensional function space is spanned by basis functions

$$
\{\phi_1,\phi_2,\ldots,\phi_n\},
$$

and that the element has $n$ degrees of freedom represented by linear functionals

$$
\{\ell_1,\ell_2,\ldots,\ell_n\}.
$$

We seek shape functions $N_i$ satisfying the usual Kronecker-delta property

$$
\ell_j(N_i)=\delta_{ij}.
$$

Each shape function can be written as a linear combination of the basis functions,

$$
N_i = \sum_{k=1}^n c_{ki}\phi_k.
$$

Applying the $j$-th DoF gives

$$
\ell_j(N_i)=\sum_{k=1}^n c_{ki}\ell_j(\phi_k).
$$

Define the matrix

$$
A_{jk} = \ell_j(\phi_k).
$$

Then the coefficients of the shape function satisfy

$$
\sum_{k=1}^n A_{jk}c_{ki} = \delta_{ji},
$$

or, in matrix form,

$$
A C = I.
$$

Consequently,

$$
C=A^{-1}.
$$

Thus, the shape functions are obtained by evaluating all DoF functionals on all basis functions, constructing the resulting matrix $A$, and solving the corresponding linear system.

This approach is completely general as long as the specified DoFs form a unisolvent set for the chosen function space. Importantly, the DoFs do not have to be point evaluations: they can be derivatives, normal derivatives, moments, or other linear functionals.

This is what makes it possible to use the same machinery for more complicated finite elements.

---

## The Morley element

The same approach can be used to generate shape functions for more specialised finite elements. We will now consider the classical [Morley element](https://defelement.org/elements/examples/triangle-morley-2.html).

The Morley element is a non-conforming triangular element with six DoFs:

- three value DoFs, one at each vertex, and
- three normal-gradient DoFs, one at the midpoint of each edge.

The element was originally developed for plate-bending problems.

The element can be described as follows:

```julia
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
    Monomial((1,0)), # x
    Monomial((2,0)), # x^2
    Monomial((0,1)), # y
    Monomial((1,1)), # x*y
    Monomial((0,2)), # y^2
)
```

In addition to the standard `PointEvaluation` at each corner, we define `NormalGradientEvaluation`, which takes the coordinate of the evaluation point and the corresponding edge normal.

The Morley function space is

$$
\mathcal{P}_2 = \{1,x,x^2,y,xy,y^2\}.
$$

As before, we can construct and print the shape functions using:

```julia
shape_values = compute_reference_shape_values(Morley{RefTriangle,2}())
print_reference_shape_values(shape_values)
```

which produces:

```julia
function Ferrite.reference_shape_value(ip::Morley{RefTriangle, 2}, ξ::Vec{2,T}, i::Int) where T
    x, y = ξ
    i==1 && return 1 - x - y + 2*x*y
    i==2 && return T(0.5)*x + T(0.5)*x^2 + T(0.5)*y - x*y - T(0.5)*y^2
    i==3 && return T(0.5)*x - T(0.5)*x^2 + T(0.5)*y - x*y + T(0.5)*y^2
    i==4 && return -T(0.7071067811865476)*x + T(0.7071067811865476)*x^2 -
                  T(0.7071067811865476)*y + T(1.4142135623730951)*x*y +
                  T(0.7071067811865476)*y^2
    i==5 && return x - x^2
    i==6 && return -y + y^2
    throw(ArgumentError("    no shape function $i for interpolation $ip"))
end
```

The result can be compared with the shape functions listed on [DefElement](https://defelement.org/elements/examples/triangle-morley-2.html). Note, however, that `Ferrite.RefTriangle` uses a different vertex and edge ordering from DefElement, so the ordering of the resulting shape functions may differ.

---

## Mapping to the physical element

Unfortunately, these shape functions cannot (currently) be implemented and used in Ferrite using the standard finite-element mapping machinery.

The reason is that the `NormalGradientEvaluation` DoFs cause the shape functions to couple when they are mapped from the reference element to a physical element.

For standard Lagrange elements, each physical shape function can be obtained directly from a corresponding reference shape function. For more general finite elements, like the Morley element, this is no longer necessarily the case.

Kirby shows that there is nevertheless a linear relationship between the physical and reference shape functions. Writing the shape functions as vectors,

$$
\mathbf{N}_x =
M\,\mathbf{N}_\xi,
$$

where

- $\mathbf{N}_\xi$ is the vector of shape functions on the reference element,
- $\mathbf{N}_x$ is the vector of shape functions on the physical element, and
- $M$ is a transformation matrix.

The matrix $M$ depends on the geometry of the physical element and on the particular finite element being transformed.

## Computing the transformation matrix
For standard Lagrange elements, this matrix reduces to the familiar one-to-one mapping between reference and physical shape functions. For elements with derivative-based DoFs, such as the Morley element, the matrix is generally non-diagonal and therefore couples the shape functions.

In this package, the transformation matrix can be computed with:

```julia
# First create a description of our element:
# a linear triangle with some arbitrary coordinates.
mapping = Mapping(
    Lagrange{RefTriangle,1}(),
    [Vec((0.01, 1.1)), Vec((1.0, 0.0)), Vec((0.1, 0.1))]
)

# Create the M matrix.
M = FerriteGenElement.build_M_matrix(Morley{RefTriangle,2}(), mapping)
```

This generates a `6 × 6` matrix such as:

```text
6×6 adjoint(::Matrix{Float64}) with eltype Float64:
 1.0  0.0  0.0   0.0         0.219512  -0.219512
 0.0  1.0  0.0  -4.44089e-17 0.0         0.219512
 0.0  0.0  1.0   4.44089e-17 -0.219512    0.0
 0.0  0.0  0.0   0.8         0.0         0.0
 0.0  0.0  0.0   0.0         0.883452    0.0
 0.0  0.0  0.0   0.0         0.0         0.883452
```

The important feature is that the normal-gradient DoFs on the edges couple with the value DoFs at the vertices. This coupling is precisely what prevents the standard scalar shape-function mapping from being used directly.

The matrix shown above is specific to the particular physical triangle defined by the coordinates

```julia
[
    Vec((0.01, 1.1)),
    Vec((1.0, 0.0)),
    Vec((0.1, 0.1))
]
```

In an actual finite-element implementation, the transformation should instead be derived analytically (or symbolically) as a function of the geometry of the physical cell. That expression can then be implemented in Ferrite so that the mapping works for every element in a mesh.

---

# Planned work

The package is still under development. Some features that would be useful to add or improve include:

- A framework for handling macro elements.
- Integral and moment DoFs.
- Analytical or symbolic generation of physical-element transformation matrices.

