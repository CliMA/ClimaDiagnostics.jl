using Test
import Base
import Dates
using Profile

using BenchmarkTools
import ProfileCanvas
import NCDatasets
import ClimaCore.Fields

import ClimaDiagnostics
import ClimaDiagnostics.Writers

include("TestTools.jl")

# The temporary directory where we write the file cannot be in /tmp, it has
# to be on disk
output_dir = mktempdir(".")

@testset "DictWriter" begin
    writer = Writers.DictWriter()

    # Test with some strings and floats instead of actual Fields and ScheduledDiagnostics
    Writers.write_field!(writer, 10.0, "mytest", nothing, nothing, 0.0)
    @test writer.dict["mytest"][0.0] == 10.0
    Writers.write_field!(writer, 20.0, "mytest", nothing, nothing, 2.0)
    @test writer.dict["mytest"][2.0] == 20.0
    Writers.write_field!(writer, 50.0, "mytest2", nothing, nothing, 8.0)
    @test writer.dict["mytest2"][8.0] == 50.0
end

@testset "NetCDFWriter" begin
    space = SphericalShellSpace()
    field = Fields.coordinate_field(space).z

    # Number of interpolation points
    NUM = 100

    writer = Writers.NetCDFWriter(
        space,
        output_dir;
        num_points = (NUM, NUM, NUM),
        sync_schedule = ClimaDiagnostics.Schedules.DivisorSchedule(2),
    )

    writer_no_vert_interpolation = Writers.NetCDFWriter(
        space,
        output_dir;
        num_points = (NUM, NUM, NUM),
        disable_vertical_interpolation = true,
    )

    # Check Base.show
    @test occursin("0 files open", "$writer")

    u = (; field)
    # FIXME: We are hardcoding the start date
    p = (; start_date = Dates.DateTime(1453, 5, 29))
    t = 10.0

    function compute!(out, u, p, t)
        if isnothing(out)
            return copy(u.field)
        else
            out .= u.field
        end
    end

    diagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
        variable = ClimaDiagnostics.DiagnosticVariable(;
            compute!,
            short_name = "ABC",
        ),
        output_short_name = "my_short_name",
        output_long_name = "My Long Name",
        output_writer = writer,
    )
    Writers.interpolate_field!(writer, field, diagnostic, u, p, t)
    Writers.write_field!(writer, field, diagnostic, u, p, t)

    @test writer.unsynced_datasets ==
          Set((writer.open_files[joinpath(output_dir, "my_short_name.nc")],))

    Writers.sync(writer)
    @test writer.unsynced_datasets == Set{NCDatasets.NCDataset}()

    NCDatasets.NCDataset(joinpath(output_dir, "my_short_name.nc")) do nc
        @test nc["ABC"].attrib["short_name"] == "ABC"
        @test nc["ABC"].attrib["long_name"] == "My Long Name"
        @test nc["ABC"].attrib["units"] == ""
        @test nc["ABC"].attrib["start_date"] ==
              string(Dates.DateTime(1453, 5, 29))
        @test size(nc["ABC"]) == (1, NUM, NUM, NUM)
        @test nc["time"][1] == 10.0
        @test nc["date"][1] ==
              string(Dates.DateTime(1453, 5, 29) + Dates.Second(10.0))
    end

    # Disable vertical interpolation
    diagnostic_novert = ClimaDiagnostics.ScheduledDiagnostic(;
        variable = ClimaDiagnostics.DiagnosticVariable(;
            compute!,
            short_name = "ABC_novert",
        ),
        output_short_name = "my_short_name_novert",
        output_long_name = "My Long Name",
        output_writer = writer_no_vert_interpolation,
    )
    Writers.interpolate_field!(
        writer_no_vert_interpolation,
        field,
        diagnostic_novert,
        u,
        p,
        t,
    )
    Writers.write_field!(
        writer_no_vert_interpolation,
        field,
        diagnostic_novert,
        u,
        p,
        t,
    )
    # Write a second time
    Writers.write_field!(
        writer_no_vert_interpolation,
        field,
        diagnostic_novert,
        u,
        p,
        t,
    )

    ###############
    # Performance #
    ###############

    # Profile interpolate
    Profile.@profile Writers.interpolate_field!(
        writer,
        field,
        diagnostic,
        u,
        p,
        t,
    )
    ProfileCanvas.html_file("flame_interpolate_netcdf.html", Profile.fetch())
    Profile.clear()

    # Profile write
    Profile.@profile Writers.write_field!(writer, field, diagnostic, u, p, t)
    ProfileCanvas.html_file("flame_write_netcdf.html", Profile.fetch())

    # Benchmark write
    timing_write_field = @benchmark Writers.write_field!(
        $writer,
        $field,
        $diagnostic,
        $u,
        $p,
        $t,
    )

    # Compare against pure NCDatasets
    function add_nc(nc, outarray, p)
        v = nc["my_short_name"]
        temporal_size, spatial_size... = size(v)
        time_index = temporal_size + 1
        t = 10.0 * time_index
        nc["time"][time_index] = t
        nc["date"][time_index] = string(p.start_date + Dates.Second(round(t)))
        v[time_index, :, :, :] = outarray
    end

    output_path = joinpath(output_dir, "clean_netcdf.nc")
    nc = NCDatasets.NCDataset(output_path, "c")
    NCDatasets.defDim(nc, "time", Inf)
    NCDatasets.defVar(nc, "time", Float64, ("time",))
    NCDatasets.defVar(nc, "date", String, ("time",))
    NCDatasets.defDim(nc, "x", NUM)
    NCDatasets.defDim(nc, "y", NUM)
    NCDatasets.defDim(nc, "z", NUM)
    v = NCDatasets.defVar(
        nc,
        "my_short_name",
        Float64,
        ("time", "x", "y", "z"),
        deflatelevel = writer.compression_level,
    )
    outarray = Array(writer.remappers["ABC"]._interpolated_values)
    v[1, :, :, :] = outarray

    timing_ncdataset = @benchmark $add_nc($nc, $outarray, $p)

    println("Our writer")
    show(stdout, MIME"text/plain"(), timing_write_field)
    println()
    println("NCDatasets")
    show(stdout, MIME"text/plain"(), timing_ncdataset)
    println()
end
