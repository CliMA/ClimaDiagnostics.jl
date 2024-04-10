import Dates

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

function SphericalShellSpace(;
    radius = 6371.0,
    height = 10.0,
    nelements = 10,
    zelem = 10,
    npolynomial = 4,
    context = ClimaComms.SingletonCommsContext(),
    FT = Float64,
)
    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(0)),
        ClimaCore.Geometry.ZPoint(FT(height));
        boundary_names = (:bottom, :top),
    )
    vertmesh = ClimaCore.Meshes.IntervalMesh(vertdomain; nelems = zelem)
    vert_center_space = ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = ClimaCore.Domains.SphereDomain(FT(radius))
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, nelements)
    horztopology = ClimaCore.Topologies.Topology2D(context, horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    return ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
end

"""
    create_problem()

An ODE problem for an exponential decay.
"""
function create_problem(space; t0 = 0.0, tf = 1.0, dt = 1e-3)
    # Let's solve an exponential decay

    Y = ClimaCore.Fields.FieldVector(; my_var = ones(space))
    p = (; tau = -0.1, start_date = Dates.DateTime(476, 9, 4))

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
