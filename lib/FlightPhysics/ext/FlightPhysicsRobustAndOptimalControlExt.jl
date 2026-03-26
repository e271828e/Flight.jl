module FlightPhysicsRobustAndOptimalControlExt

using RobustAndOptimalControl
using ComponentArrays
using FlightPhysics.Linearization


function RobustAndOptimalControl.named_ss(lss::LinearizedSS)
    x_labels, u_labels, y_labels = map(collect ∘ propertynames, (lss.x0, lss.u0, lss.y0))
    named_ss(ss(lss), x = x_labels, u = u_labels, y = y_labels)
end


function Linearization.subsystem(nss::RobustAndOptimalControl.NamedStateSpace;
                  x = nss.x, u = nss.u, y = nss.y)

    #to do: generalize for scalars

    x_axis = Axis(nss.x)
    u_axis = Axis(nss.u)
    y_axis = Axis(nss.y)

    A_nss = ComponentMatrix(nss.A, x_axis, x_axis)
    B_nss = ComponentMatrix(nss.B, x_axis, u_axis)
    C_nss = ComponentMatrix(nss.C, y_axis, x_axis)
    D_nss = ComponentMatrix(nss.D, y_axis, u_axis)

    A_sub = A_nss[x, x]
    B_sub = B_nss[x, u]
    C_sub = C_nss[y, x]
    D_sub = D_nss[y, u]

    ss_sub = ss(A_sub, B_sub, C_sub, D_sub)

    named_ss(ss_sub; x, u, y)

end




end
