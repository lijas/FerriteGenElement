
struct Morley{RefTriangle,order} <: Ferrite.Interpolation{RefTriangle,order}
end

FerriteGenElement.interpolation_dofs(::Morley{RefTriangle,2}) = (
    ValueDof(Vec((0.0, 0.0))),
    ValueDof(Vec((1.0, 0.0))),
    ValueDof(Vec((0.0, 1.0))), 
    NormalGradientDof(Vec(1/2, 1/2), Vec(1/√2, 1/√2)),
    NormalGradientDof(Vec(0.0, 1/2), Vec(1.0, 0.0)),
    NormalGradientDof(Vec(1/2, 0.0), Vec((0.0,-1.0))),
)
FerriteGenElement.monomials(::Morley{RefTriangle,2}) = (
    Monomial((0,0)), # 1
    Monomial((1,0)), # x^1
    Monomial((2,0)), # x^2
    Monomial((0,1)), # y
    Monomial((1,1)), # x*y
    Monomial((0,2)), # y^2
)

shape_values = compute_reference_shape_values(Morley{RefTriangle,2}())
print_reference_shape_values(shape_values)

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

ip = Lagrange{RefTriangle,1}()
coords = [Vec((0.0, 1.0)), Vec((1.0, 0.0)), Vec((0.1, 0.1))]
mapping = Mapping(ip, coords)
B = FerriteGenElement.build_B_matrix(Morley{RefTriangle,2}(), mapping)
M = inv(B)'