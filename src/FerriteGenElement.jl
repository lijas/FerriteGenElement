module FerriteGenElement


using Reexport: @reexport
@reexport using Tensors
@reexport using Ferrite: Ferrite, Interpolation, getnbasefunctions
using LinearAlgebra: Transpose

#using Symbolics

abstract type Dof end

"""
PointEvaluation
"""
struct PointEvaluation <: Dof
    ξ::Vec{2,Float64}
end
apply(::PointEvaluation, f, x) = f(x)
apply2(a::PointEvaluation, f) = f(a.ξ)
ncomponents(::PointEvaluation) = 1

"""
DerivativeEvaluation
"""
struct DerivativeEvaluation <: Dof
    ξ::Vec{2,Float64}
    dir::Vec{2,Float64}
end
apply(::DerivativeEvaluation, f, x, dir) = Tensors.gradient(f, x) ⋅ dir
apply2(a::DerivativeEvaluation, f) = Tensors.gradient(f, a.ξ) ⋅ a.dir
ncomponents(::DerivativeEvaluation) = 1

XDerivativeEvaluation(x) = DerivativeEvaluation(x, Vec((1.0, 0.0)))
YDerivativeEvaluation(x) = DerivativeEvaluation(x, Vec((0.0, 1.0)))

"""
NormalGradientEvaluation
"""
struct NormalGradientEvaluation <: Dof
    ξ::Vec{2,Float64}
    normal::Vec{2,Float64}
end
apply(::NormalGradientEvaluation, f, x, normal) = Tensors.gradient(f, x) ⋅ normal
apply2(a::NormalGradientEvaluation, f) = Tensors.gradient(f, a.ξ) ⋅ a.normal
ncomponents(::NormalGradientEvaluation) = 1

"""
TangentGradientEvaluation
"""
struct TangentGradientEvaluation <: Dof
    ξ::Vec{2,Float64}
    tangent::Vec{2,Float64}
end
apply(::TangentGradientEvaluation, f, x, tangent) = Tensors.gradient(f, x) ⋅ tangent
apply2(a::TangentGradientEvaluation, f) = Tensors.gradient(f, a.ξ) ⋅ a.tangent
ncomponents(::TangentGradientEvaluation) = 1

"""
HessianEvaluation
"""
struct HessianEvaluation <: Dof
    ξ::Vec{2,Float64}
    dir1::Vec{2,Float64}
    dir2::Vec{2,Float64}
end
apply(::HessianEvaluation, f, x, dir1, dir2) = dir1 ⋅ Tensors.hessian(f, x) ⋅ dir2
apply2(a::HessianEvaluation, f) = a.dir1 ⋅ Tensors.hessian(f, a.ξ) ⋅ a.dir2
ncomponents(::HessianEvaluation) = 1

XXHessianEvaluation(x) = HessianEvaluation(x, Vec((1.0, 0.0)), Vec((1.0, 0.0)))
YYHessianEvaluation(x) = HessianEvaluation(x, Vec((0.0, 1.0)), Vec((0.0, 1.0)))
XYHessianEvaluation(x) = HessianEvaluation(x, Vec((1.0, 0.0)), Vec((0.0, 1.0)))

#####
# Mapping
#####
struct Mapping{T}
    ip # Interpolation
    coords::Vector{Vec{2,T}}
end

function compute_jacobian(mapping::Mapping, ξ::Vec{dim, T}) where {dim,T}
    N = getnbasefunctions(mapping.ip)
    J = zero(Tensor{2,dim,T})
    for i in 1:N
        dNdξ = Ferrite.reference_shape_gradient(mapping.ip, ξ, i)
        J += dNdξ ⊗ mapping.coords[i]
    end
    return J
end

#####
# Push forward operations for each DoF
#####
function push_forward(dof::PointEvaluation, ip::Interpolation, shape_nr::Int, ::Mapping)
    return Ferrite.reference_shape_value(ip, dof.ξ, shape_nr)
end
function push_forward(dof::DerivativeEvaluation, ip::Interpolation, shape_nr::Int, mapping::Mapping)
    # Chain rule: ∇_x N = ∇_ξ N ⋅ J⁻¹
    ∇N = Ferrite.reference_shape_gradient(ip, dof.ξ, shape_nr)
    J = compute_jacobian(mapping, dof.ξ)
    ∇x_N = ∇N ⋅ inv(J)
    
    return ∇x_N ⋅ dof.dir
end
function push_forward(dof::NormalGradientEvaluation, ip::Interpolation, shape_nr::Int, mapping::Mapping)
    ∇N = Ferrite.reference_shape_gradient(ip, dof.ξ, shape_nr)
    J = compute_jacobian(mapping, dof.ξ)
    ∇x_N = ∇N ⋅ inv(J)
    
    # Transform the normal
    normal_ref = inv(J)' ⋅ dof.normal
    normal = normal_ref / norm(normal_ref)
    
    return ∇x_N ⋅ normal
end

"""
Monomial
"""
struct Monomial{N}
    exponents::NTuple{N,Int}
end
function (m::Monomial{N})(ξ::Vec{N,T}) where N where T
    out = one(T)
    for i in 1:N
        out *= ξ[i]^m.exponents[i]
    end
    return out
end

include("element_interface.jl")
include("create_B_matrix.jl")
include("compute_shape_functions.jl")

export PointEvaluation, DerivativeEvaluation, HessianEvaluation, NormalGradientEvaluation, TangentGradientEvaluation
export Mapping
export Monomial
export compute_reference_shape_values, print_reference_shape_values, build_B_matrix, build_M_matrix
end # module FerriteGenElement
