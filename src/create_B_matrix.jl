
function build_B_matrix(element::Interpolation, mapping::Mapping{T}) where T
    dofs = interpolation_dofs(element)
    N = length(dofs)
    M = N # getnbasefunctions(element)
    B = zeros(T, N, M)
    for i in 1:N, j in 1:M
        node = dofs[i]
        tmp = push_forward(node, element, j, mapping)
        B[i,j] = tmp
    end
    return B
end

build_M_matrix(element::Interpolation, mapping::Mapping) = transpose(inv(build_B_matrix(element, mapping)))