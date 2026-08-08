
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

shape_values = compute_reference_shape_values(MyLinearTriangle())
print_reference_shape_values(shape_values)
