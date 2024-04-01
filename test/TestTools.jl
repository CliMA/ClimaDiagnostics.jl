import SciMLBase

import ClimaCore
import ClimaComms
import ClimaTimeSteppers

function ColumnCenterFiniteDifferenceSpace(
    zelem = 10,
    context = ClimaComms.SingletonCommsContext();
    FT = Float64,
)
    zlim = (FT(0.0), FT(1.0))
    domain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(zlim[1]),
        ClimaCore.Geometry.ZPoint(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    mesh = ClimaCore.Meshes.IntervalMesh(domain, nelems = zelem)
    topology = ClimaCore.Topologies.IntervalTopology(context, mesh)
    return ClimaCore.Spaces.CenterFiniteDifferenceSpace(topology)
end

"""
    create_problem_algo()

An ODE problem for an exponential decay.
"""
function create_problem(; t0 = 0.0, tf = 1.0, dt = 1e-3)
    # Let's solve an exponential decay
    space = ColumnCenterFiniteDifferenceSpace()

    Y = ClimaCore.Fields.FieldVector(; my_var = ones(space))
    p = (; tau = -0.1)

    function exp_tendency!(dY, Y, p, t)
        @. dY.my_var = p.tau * Y.my_var
    end

    prob = SciMLBase.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(T_exp! = exp_tendency!),
        Y,
        (t0, tf),
        p,
    )
    algo = ClimaTimeSteppers.ExplicitAlgorithm(ClimaTimeSteppers.RK4())

    args = prob, algo
    kwargs = Dict(:dt => dt)

    return args, kwargs
end
