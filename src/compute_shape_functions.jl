
"""
    ShapeValues

Stores the computed shape functions for an interpolation.
Shape function can be printing with print_shape_functions(::ShapeValues)
"""
struct ShapeValues{IP}
    ip::IP # <: Ferrite.Interpolation
    C::Transpose{Float64, Matrix{Float64}} # Coefficients for the shape functions
end

function _coeff_string(coeff; atol=1e-12)
    T = typeof(coeff)
    # Check if approximately integer
    n = round(coeff)
    if isapprox(coeff, n; atol=atol)
        return string(Int(n))
    else
        coeff = coeff# round(coeff, sigdigits=5)
        return "T(" * string(coeff) * ")"
    end
end

function _monomial_string(a::Int, b::Int)
    if a == 0 && b == 0
        return "1"
    end

    xpart = a == 0 ? "" : a == 1 ? "x" : "x^$a"
    ypart = b == 0 ? "" : b == 1 ? "y" : "y^$b"

    if xpart == ""
        return ypart
    elseif ypart == ""
        return xpart
    else
        return xpart * "*" * ypart
    end
end

function _term_string(coeff, a::Int, b::Int)
    monomials_string = _monomial_string(a, b)

    if monomials_string == "1"
        return isapprox(coeff, 1.0; atol=1e-12) ? "1" : _coeff_string(coeff)
    end

    if isapprox(coeff, 1.0; atol=1e-12)
        return monomials_string
    else
        return _coeff_string(coeff) * "*" * monomials_string
    end
end

"""
    construct_reference_basis(element::Ferrite.Interpolation)

Computes the shape values for an given interpolation/element.
The function computes a coefficent matrix C, for which the shape values are obtained 
as "shape_functions = C * monos(element)"
"""
function compute_reference_shape_values(element::Ferrite.Interpolation)
    dofs = interpolation_dofs(element)
    monos = monomials(element)
    N = sum(d -> ncomponents(d), dofs)
    P = zeros(Float64, N, N)
    i = 1
    for dof in dofs
        j=1
        for mono in monos
            v = apply2(dof, mono)
            P[i, j] = v
            j += 1
        end
        i += ncomponents(dof)
    end
    
    C = transpose(inv(P))

    return ShapeValues(element, C)
end

function print_reference_shape_values(sv::ShapeValues)

    C = sv.C 
    element = sv.ip
    monos = monomials(element)
    N = sum(d -> ncomponents(d), interpolation_dofs(element))

    shape_functions = String[]
    for row in 1:N
        f = ""
        first_time = true
        for col in 1:N
            a, b = monos[col].exponents
            if !(isapprox(C[row,col], 0.0, atol = 1e-13)) #TODO: better tolerance handling
                coeff = C[row,col]
                term = _term_string(abs(coeff), a, b)
                if first_time
                    first_time = false
                    f *= coeff < 0 ? "-" * term : term
                else
                    f *= coeff < 0 ? " - " : " + "
                    f *= term
                end
            end
        end
        push!(shape_functions, f)
    end

    println("function Ferrite.reference_shape_value(ip::$(typeof(element)), ξ::Vec{2,T}, i::Int) where T")
    println("    x, y = ξ")
    for i in eachindex(shape_functions)
        println("    i==$(i) && return $(shape_functions[i])")
    end
    println("    throw(ArgumentError(\"    no shape function \$i for interpolation \$ip\"))")
    println("end")

end