using Test
import Base
import Dates
using Profile

using BenchmarkTools
import ProfileCanvas
import NCDatasets
import ClimaCore
import ClimaCore.Fields
import ClimaCore.Spaces
import ClimaCore.Geometry
import ClimaComms

import ClimaDiagnostics
import ClimaDiagnostics.Writers

include("TestTools.jl")

# The temporary directory where we write the file cannot be in /tmp, it has
# to be on disk
output_dir = mktempdir(pwd())

@testset "DictWriter" begin
    writer = Writers.DictWriter()

    # Test with some strings and floats instead of actual Fields and ScheduledDiagnostics
    Writers.write_field!(writer, 10.0, "mytest", nothing, nothing, 0.0)
    @test writer.dict["mytest"][0.0] == 10.0
    Writers.write_field!(writer, 20.0, "mytest", nothing, nothing, 2.0)
    @test writer.dict["mytest"][2.0] == 20.0
    Writers.write_field!(writer, 50.0, "mytest2", nothing, nothing, 8.0)
    @test writer.dict["mytest2"][8.0] == 50.0

    @test issorted(writer.dict["mytest"])
end

@testset "NetCDFWriter" begin
    space = SphericalShellSpace()
    field = Fields.coordinate_field(space).z

    # Number of interpolation points
    NUM = 50

    writer = Writers.NetCDFWriter(
        space,
        output_dir;
        num_points = (NUM, 2NUM, 3NUM),
        sync_schedule = ClimaDiagnostics.Schedules.DivisorSchedule(2),
        z_sampling_method = ClimaDiagnostics.Writers.FakePressureLevelsMethod(),
    )

    writer_no_vert_interpolation = Writers.NetCDFWriter(
        space,
        output_dir;
        num_points = (NUM, 2NUM, 3NUM),
        z_sampling_method = ClimaDiagnostics.Writers.LevelsMethod(),
    )

    # Check Base.show
    @test occursin("0 files open", "$writer")

    u = (; field)
    # FIXME: We are hardcoding the start date
    p = (; start_date = Dates.DateTime(1453, 5, 29))
    t = 10.0

    function compute!(out, u, p, t)
        if isnothing(out)
            return u.field
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
        @test size(nc["ABC"]) == (1, NUM, 2NUM, 3NUM)
        @test nc["time"][1] == 10.0
        @test nc["date"][1] ==
              string(Dates.DateTime(1453, 5, 29) + Dates.Second(10.0))
        @test nc["time"].attrib["standard_name"] == "time"
        @test nc["time"].attrib["long_name"] == "Time"
        @test nc["time"].attrib["axis"] == "T"
        @test nc["lon"].attrib["standard_name"] == "longitude"
        @test nc["lon"].attrib["long_name"] == "Longitude"
        @test nc["lon"].attrib["axis"] == "X"
        @test nc["lat"].attrib["standard_name"] == "latitude"
        @test nc["lat"].attrib["long_name"] == "Latitude"
        @test nc["lat"].attrib["axis"] == "Y"
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

    # Check boxes
    xyboxspace = BoxSpace()
    xyboxfield = Fields.coordinate_field(xyboxspace).z

    xyboxwriter = Writers.NetCDFWriter(
        xyboxspace,
        output_dir;
        num_points = (NUM, 2NUM, 3NUM),
    )
    xyboxdiagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
        variable = ClimaDiagnostics.DiagnosticVariable(;
            compute!,
            short_name = "ABC",
        ),
        output_short_name = "my_short_name_xybox",
        output_long_name = "My Long Name",
        output_writer = xyboxwriter,
    )
    xyboxu = (; xyboxfield)
    Writers.interpolate_field!(
        xyboxwriter,
        xyboxfield,
        xyboxdiagnostic,
        xyboxu,
        p,
        t,
    )
    Writers.write_field!(xyboxwriter, xyboxfield, xyboxdiagnostic, xyboxu, p, t)

    longlatboxspace = BoxSpace(; lonlat = true)
    longlatboxfield = Fields.coordinate_field(longlatboxspace).z

    longlatboxwriter = Writers.NetCDFWriter(
        longlatboxspace,
        output_dir;
        num_points = (NUM, 2NUM, 3NUM),
    )
    longlatboxdiagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
        variable = ClimaDiagnostics.DiagnosticVariable(;
            compute!,
            short_name = "ABC",
        ),
        output_short_name = "my_short_name_longlatbox",
        output_long_name = "My Long Name",
        output_writer = longlatboxwriter,
    )
    longlatboxu = (; longlatboxfield)
    Writers.interpolate_field!(
        longlatboxwriter,
        longlatboxfield,
        longlatboxdiagnostic,
        longlatboxu,
        p,
        t,
    )
    Writers.write_field!(
        longlatboxwriter,
        longlatboxfield,
        longlatboxdiagnostic,
        longlatboxu,
        p,
        t,
    )
    # Write a second time, to check consistency
    Writers.write_field!(
        longlatboxwriter,
        longlatboxfield,
        longlatboxdiagnostic,
        longlatboxu,
        p,
        t,
    )

    # Check columns
    if pkgversion(ClimaCore) >= v"0.14.23"
        # Center space
        for (i, colspace) in enumerate((
            ColumnCenterFiniteDifferenceSpace(),
            ColumnFaceFiniteDifferenceSpace(),
        ))
            colfield = Fields.coordinate_field(colspace).z

            colwriter =
                Writers.NetCDFWriter(colspace, output_dir; num_points = (NUM,))
            coldiagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
                variable = ClimaDiagnostics.DiagnosticVariable(;
                    compute!,
                    short_name = "ABC",
                ),
                output_short_name = "my_short_name_c$(i)",
                output_long_name = "My Long Name",
                output_writer = colwriter,
            )
            colu = (; colfield)
            Writers.interpolate_field!(
                colwriter,
                colfield,
                coldiagnostic,
                colu,
                p,
                t,
            )
            Writers.write_field!(colwriter, colfield, coldiagnostic, colu, p, t)
            # Write a second time, to check consistency
            Writers.write_field!(colwriter, colfield, coldiagnostic, colu, p, t)
        end
    end

    ###############
    # Point Space #
    ###############
    point_val = 3.14
    point_space =
        Spaces.PointSpace(ClimaComms.context(), Geometry.ZPoint(point_val))
    point_field = Fields.coordinate_field(point_space)
    point_writer = Writers.NetCDFWriter(point_space, output_dir)

    point_u = (; field = point_field)

    point_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
        variable = ClimaDiagnostics.DiagnosticVariable(;
            compute!,
            short_name = "ABC",
        ),
        output_short_name = "my_short_name_point",
        output_long_name = "My Long Name Point",
        output_writer = point_writer,
    )

    Writers.interpolate_field!(
        point_writer,
        point_field,
        point_diagnostic,
        point_u,
        p,
        t,
    )
    Writers.write_field!(
        point_writer,
        point_field,
        point_diagnostic,
        point_u,
        p,
        t,
    )
    # Write a second time
    Writers.write_field!(
        point_writer,
        point_field,
        point_diagnostic,
        point_u,
        p,
        t,
    )
    close(point_writer)

    NCDatasets.NCDataset(joinpath(output_dir, "my_short_name.nc")) do nc
        nc["ABC"] == [point_val, point_val]
    end

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
    NCDatasets.defDim(nc, "y", 2NUM)
    NCDatasets.defDim(nc, "z", 3NUM)
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
