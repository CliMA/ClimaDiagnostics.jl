using Test
import ClimaDiagnostics
import ClimaDiagnostics.Schedules
import ClimaDiagnostics.ScheduledDiagnostics
import ClimaDiagnostics.Writers

import NCDatasets

import SciMLBase

include("TestTools.jl")

output_dir = mktempdir(pwd())

@testset "Pressure diagnostics" begin
    NUM = 10
    t0 = 0.0
    tf = 1.0
    dt = 0.1
    space = SphericalShellSpace()
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

    function compute_field!(out, u, p, t)
        if isnothing(out)
            return Y.my_var
        else
            out .= Y.my_var
            return nothing
        end
    end

    # A simple diagnostic that just saves the time at every timestep
    simple_var = ClimaDiagnostics.DiagnosticVariable(;
        compute! = compute_field!,
        short_name = "YO",
        long_name = "YO YO",
    )

    netcdf_writer = Writers.NetCDFWriter(
        space,
        output_dir;
        num_points = (NUM, 2NUM, 3NUM),
        sync_schedule = ClimaDiagnostics.Schedules.DivisorSchedule(2),
        z_sampling_method = ClimaDiagnostics.Writers.FakePressureLevelsMethod(),
    )

    time_reduction_diagnostic_every_step = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = netcdf_writer,
        reduction_time_func = max,
        output_short_name = "yo_max",
    )

    inst_diagnostic_every_step = ClimaDiagnostics.ScheduledDiagnostic(
        variable = simple_var,
        output_writer = netcdf_writer,
        output_short_name = "yo_inst",
    )

    pfull_levels = ClimaDiagnostics.era5_pressure_levels()
    pfull_diagnostic_handler = ClimaDiagnostics.PfullCoordsDiagnosticsHandler(
        [time_reduction_diagnostic_every_step, inst_diagnostic_every_step],
        Y,
        p,
        t0,
        compute_field!;
        dt,
        pfull_levels,
    )
    diagnostic_handler = ClimaDiagnostics.DiagnosticsHandler(
        [time_reduction_diagnostic_every_step, inst_diagnostic_every_step],
        Y,
        p,
        t0;
        dt,
    )

    pfull_diag_cb =
        ClimaDiagnostics.DiagnosticsCallback(pfull_diagnostic_handler)
    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    SciMLBase.solve(
        prob,
        algo,
        dt = dt,
        callback = SciMLBase.CallbackSet((), (diag_cb, pfull_diag_cb)),
    )

    # TODO: SUPER IMPORTANT
    # TODO: SUPER IMPORTANT
    # TODO: SUPER IMPORTANT
    # TODO: Test that pressure levels are actually being computed the way that I think it is
    #       I don't know how to do this

    close(netcdf_writer)

    @test isfile(joinpath(output_dir, "_pfull_coords", "yo_inst.nc"))

    NCDatasets.NCDataset(
        joinpath(output_dir, "_pfull_coords", "yo_max.nc"),
    ) do ds
        num_time = length(ds["time"])
        num_pfull_levels = length(ds["pressure_level"])
        num_horizontal_indices = length(ds["horizontal_index"])
        @test size(ds["YO"]) ==
              (num_time, num_pfull_levels, num_horizontal_indices)

        lats = Array(ds["lat"])
        lons = Array(ds["lon"])
        @test all(lat -> -90.0 <= lat <= 90.0, lats)
        @test all(lon -> -180.0 <= lon <= 180.0, lons)

        # Testing interpolate with layers instead...
        # import ClimaCore
        # single_layer = ds["YO"][1,1,:]

        # pfull_field = pfull_diagnostic_handler.pfull_field
        # pfull_field_level = ClimaCore.level(pfull_field, 1)
        # values = ClimaCore.DataLayouts.array2data(single_layer, ClimaCore.Fields.field_values(pfull_field_level))
        # field = ClimaCore.Fields.Field(values, axes(pfull_field_level))

        # longpts = range(-180.0, 180.0, 21)
        # latpts = range(-80.0, 80.0, 21)
        # hcoords = [ClimaCore.Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]

        # interpolated_array = ClimaCore.Remapping.interpolate(field, hcoords, nothing)


        # Test level conversion
        # import Random
        # pfull_field = pfull_diagnostic_handler.pfull_field
        # rand!(parent(pfull_field))

        # pfull_field_level = ClimaCore.level(pfull_field, 1)

        # reshape_to_cols(f) = reshape(parent(f), ClimaCore.Spaces.nlevels(axes(f)), ClimaCore.Spaces.ncolumns(axes(f)))
        # cols = reshape_to_cols(pfull_field)
        # first_level = cols[1,:]

        # values = ClimaCore.DataLayouts.array2data(first_level, ClimaCore.Fields.field_values(pfull_field_level))
        # same_as_pfull_field_level = ClimaCore.Fields.Field(values, axes(pfull_field_level))
    end
end
