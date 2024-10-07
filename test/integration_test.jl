using Test
using Profile
using ProfileCanvas

import SciMLBase
import NCDatasets

import ClimaDiagnostics

import ClimaComms
@static if pkgversion(ClimaComms) >= v"0.6"
    ClimaComms.@import_required_backends
end

const context = ClimaComms.context()
ClimaComms.init(context)

include("TestTools.jl")

"""
Set up a full test problem

Increasing `more_compute_diagnostics` adds more copies of a compute diagnostic with no output.
Useful to stress allocations.
"""
function setup_integrator(output_dir; context, more_compute_diagnostics = 0)
    t0 = 0.0
    tf = 10.0
    dt = 1.0
    space = SphericalShellSpace(; context)
    args, kwargs = create_problem(space; t0, tf, dt)

    @info "Writing output to $output_dir"

    dummy_writer = ClimaDiagnostics.Writers.DummyWriter()
    h5_writer = ClimaDiagnostics.Writers.HDF5Writer(output_dir)
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        space,
        output_dir;
        num_points = (10, 5, 3),
    )

    function compute_my_var!(out, u, p, t)
        if isnothing(out)
            return u.my_var
        else
            out .= u.my_var
            return nothing
        end
    end

    simple_var = ClimaDiagnostics.DiagnosticVariable(;
        compute! = compute_my_var!,
        short_name = "YO",
        long_name = "YO YO",
    )

    average_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
        reduction_time_func = (+),
        output_schedule_func = ClimaDiagnostics.Schedules.DivisorSchedule(2),
        pre_output_hook! = ClimaDiagnostics.average_pre_output_hook!,
    )
    inst_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
    )
    inst_every3s_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = nc_writer,
        output_schedule_func = ClimaDiagnostics.Schedules.EveryDtSchedule(
            3.0,
            t_last = t0,
        ),
    )
    inst_diagnostic_h5 = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = h5_writer,
    )
    scheduled_diagnostics = [
        average_diagnostic,
        inst_diagnostic,
        inst_diagnostic_h5,
        inst_every3s_diagnostic,
    ]

    # Add more weight, useful for stressing allocations
    compute_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = dummy_writer,
    )
    scheduled_diagnostics = [
        scheduled_diagnostics...,
        [compute_diagnostic for _ in 1:more_compute_diagnostics]...,
    ]

    return ClimaDiagnostics.IntegratorWithDiagnostics(
        SciMLBase.init(args...; kwargs...),
        scheduled_diagnostics,
    )
end

@testset "A full problem" begin
    mktempdir() do output_dir
        output_dir = ClimaComms.bcast(context, output_dir)

        integrator = setup_integrator(output_dir; context)

        SciMLBase.solve!(integrator)

        if ClimaComms.iamroot(context)
            NCDatasets.NCDataset(joinpath(output_dir, "YO_1it_inst.nc")) do nc
                @test nc["YO"].attrib["short_name"] == "YO"
                @test nc["YO"].attrib["long_name"] == "YO YO, Instantaneous"
                @test size(nc["YO"]) == (11, 10, 5, 10)
            end

            NCDatasets.NCDataset(
                joinpath(output_dir, "YO_2it_average.nc"),
            ) do nc
                @test nc["YO"].attrib["short_name"] == "YO"
                @test nc["YO"].attrib["long_name"] ==
                      "YO YO, average within every 2 iterations"
                @test size(nc["YO"]) == (5, 10, 5, 10)
            end

            NCDatasets.NCDataset(joinpath(output_dir, "YO_3s_inst.nc")) do nc
                @test nc["YO"].attrib["short_name"] == "YO"
                @test nc["YO"].attrib["long_name"] == "YO YO, Instantaneous"
                @test size(nc["YO"]) == (4, 10, 5, 10)
            end
        end
    end
end

@testset "Performance" begin
    mktempdir() do output_dir
        output_dir = ClimaComms.bcast(context, output_dir)

        # Flame
        integrator = setup_integrator(output_dir; context)
        prof = Profile.@profile SciMLBase.solve!(integrator)
        ClimaComms.iamroot(context) && (results = Profile.fetch())
        ClimaComms.iamroot(context) &&
            ProfileCanvas.html_file("flame.html", results)

        # Allocations
        integrator = setup_integrator(output_dir; context)
        prof = Profile.Allocs.@profile SciMLBase.solve!(integrator)
        ClimaComms.iamroot(context) && (results = Profile.Allocs.fetch())
        ClimaComms.iamroot(context) &&
            (allocs = ProfileCanvas.view_allocs(results))
        ClimaComms.iamroot(context) &&
            ProfileCanvas.html_file("allocs.html", allocs)
    end
end
