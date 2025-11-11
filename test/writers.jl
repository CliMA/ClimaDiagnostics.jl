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
import ClimaCore.CommonSpaces
import ClimaCore.Meshes
import ClimaComms

import ClimaDiagnostics
import ClimaDiagnostics.Writers

import ClimaUtilities.TimeManager: ITime

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

@testset "Target coordinates" begin
    NUM = 50

    # The order returned from target_coordinates should always be ((lon, lat), vertical) or ((x, y), vertical).
    # As a result, hcoords_from_horizontal_space should return a matrix of LatLongPoint or XYPoint where each axis of the matrix corresponds to lon and lat or x and y respectively.
    # Test spherical shell space
    space = SphericalShellSpace()
    hpts, vpts = Writers.target_coordinates(
        space,
        (NUM, 2NUM, 3NUM),
        ClimaDiagnostics.Writers.FakePressureLevelsMethod(),
    )
    lons, lats = hpts
    @test lons == range(-180.0, 180.0, NUM)
    @test lats == range(-90.0, 90.0, 2NUM)
    # It is a bit difficult to test the results of FakePressureLevelsMethod, so
    # we just check the length of vpts
    @test length(vpts) == 3NUM

    horizontal_space = Spaces.horizontal_space(space)
    hcoords = Writers.hcoords_from_horizontal_space(
        horizontal_space,
        Meshes.domain(Spaces.topology(horizontal_space)),
        hpts,
    )
    @test size(hcoords) == (NUM, 2NUM)
    @test hcoords ==
          [Geometry.LatLongPoint(lat, lon) for lon in lons, lat in lats]

    hpts, vpts = Writers.target_coordinates(
        space,
        (NUM, 2NUM, 3NUM),
        ClimaDiagnostics.Writers.LevelsMethod(),
    )
    lons, lats = hpts
    @test lons == range(-180.0, 180.0, NUM)
    @test lats == range(-90.0, 90.0, 2NUM)
    # LevelsMethod override 3NUM to choose points that correspond to the center of the cells
    @test length(vpts) == 10

    # Test BoxSpace with lonlat = false
    xyboxspace = BoxSpace(; ylim = (-Float64(1), Float64(2)))
    hpts, vpts = Writers.target_coordinates(
        xyboxspace,
        (NUM, 2NUM, 3NUM),
        ClimaDiagnostics.Writers.LevelsMethod(),
    )
    xpts, ypts = hpts
    @test xpts == range(-1, 1, NUM)
    @test ypts == range(-1, 2, 2NUM)
    @test length(vpts) == 10

    horizontal_space = Spaces.horizontal_space(xyboxspace)
    hcoords = Writers.hcoords_from_horizontal_space(
        horizontal_space,
        Meshes.domain(Spaces.topology(horizontal_space)),
        hpts,
    )
    @test size(hcoords) == (NUM, 2NUM)
    @test hcoords == [Geometry.XYPoint(x, y) for x in xpts, y in ypts]

    # Test BoxSpace with lonlat = true
    longlatboxspace =
        BoxSpace(; lonlat = true, ylim = (-Float64(1), Float64(2)))
    hpts, vpts = Writers.target_coordinates(
        longlatboxspace,
        (NUM, 2NUM, 3NUM),
        ClimaDiagnostics.Writers.LevelsMethod(),
    )
    lons, lats = hpts
    @test lons == range(-1, 2, NUM)
    @test lats == range(-1, 1, 2NUM)
    @test length(vpts) == 10

    horizontal_space = Spaces.horizontal_space(longlatboxspace)
    hcoords = Writers.hcoords_from_horizontal_space(
        horizontal_space,
        Meshes.domain(Spaces.topology(horizontal_space)),
        hpts,
    )
    @test size(hcoords) == (NUM, 2NUM)
    @test hcoords ==
          [Geometry.LatLongPoint(lat, lon) for lon in lons, lat in lats]

    # Test column spaces
    colcenterspace = ColumnCenterFiniteDifferenceSpace()
    vpts = Writers.target_coordinates(
        colcenterspace,
        (NUM,),
        ClimaDiagnostics.Writers.LevelsMethod(),
    )
    @test length(vpts) == 10

    colfacespace = ColumnFaceFiniteDifferenceSpace()
    vpts = Writers.target_coordinates(
        colfacespace,
        (NUM,),
        ClimaDiagnostics.Writers.LevelsMethod(),
    )
    @test length(vpts) == 10

    # Test horizontal space
    horizontal_space = ClimaCore.Spaces.level(space, 1)
    hpts = Writers.target_coordinates(horizontal_space, (NUM, 2NUM))
    lons, lats = hpts
    @test lons == range(-180, 180, NUM)
    @test lats == range(-90, 90, 2NUM)

    hcoords = Writers.hcoords_from_horizontal_space(
        horizontal_space,
        Meshes.domain(Spaces.topology(horizontal_space)),
        hpts,
    )
    @test size(hcoords) == (NUM, 2NUM)
    @test hcoords ==
          [Geometry.LatLongPoint(lat, lon) for lon in lons, lat in lats]
end

@testset "NetCDFWriter" begin
    @testset "default_num_points" begin
        @test Writers.default_num_points(
            CommonSpaces.ExtrudedCubedSphereSpace(;
                z_elem = 10,
                z_min = 0,
                z_max = 1,
                radius = 10,
                h_elem = 10,
                n_quad_points = 4,
                staggering = CommonSpaces.CellCenter(),
            ),
        ) == (120, 60, 10)

        @test Writers.default_num_points(
            CommonSpaces.SliceXZSpace(;
                z_elem = 10,
                x_min = 0,
                x_max = 1,
                z_min = 0,
                z_max = 1,
                periodic_x = false,
                n_quad_points = 4,
                x_elem = 4,
                staggering = CommonSpaces.CellCenter(),
            ),
        ) == (12, 10)
        @test Writers.default_num_points(
            CommonSpaces.Box3DSpace(;
                z_elem = 10,
                x_min = 0,
                x_max = 1,
                y_min = 0,
                y_max = 1,
                z_min = 0,
                z_max = 10,
                periodic_x = false,
                periodic_y = false,
                n_quad_points = 4,
                x_elem = 3,
                y_elem = 4,
                staggering = CommonSpaces.CellCenter(),
            ),
        ) == (9, 12, 10)
    end

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
        @test nc["date"][1] == Dates.DateTime(1453, 5, 29) + Dates.Second(10.0)
        @test nc["time"].attrib["standard_name"] == "time"
        @test nc["time"].attrib["long_name"] == "Time"
        @test nc["time"].attrib["axis"] == "T"
        @test nc["lon"].attrib["standard_name"] == "longitude"
        @test nc["lon"].attrib["long_name"] == "Longitude"
        @test nc["lon"].attrib["axis"] == "X"
        @test nc["lat"].attrib["standard_name"] == "latitude"
        @test nc["lat"].attrib["long_name"] == "Latitude"
        @test nc["lat"].attrib["axis"] == "Y"

        # Test dimensions
        hpts, vpts = Writers.target_coordinates(
            space,
            (NUM, 2NUM, 3NUM),
            ClimaDiagnostics.Writers.FakePressureLevelsMethod(),
        )
        lon, lat = hpts
        @test nc["lon"][:] == lon
        @test nc["lat"][:] == lat
        @test nc["z"][:] == vpts

        # Test bounds
        @test nc["time_bnds"][:, 1] == [0.0; 10.0]
        @test nc["date_bnds"][:, 1] == [
            Dates.DateTime(1453, 5, 29)
            Dates.DateTime(1453, 5, 29) + Dates.Second(10.0)
        ]
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

    NCDatasets.NCDataset(joinpath(output_dir, "my_short_name_novert.nc")) do nc
        # Test dimensions
        hpts, vpts = Writers.target_coordinates(
            space,
            (NUM, 2NUM, 3NUM),
            ClimaDiagnostics.Writers.LevelsMethod(),
        )
        lon, lat = hpts
        @test nc["lon"][:] == lon
        @test nc["lat"][:] == lat
        @test nc["z"][:] == vpts
    end

    # Check boxes
    xyboxspace = BoxSpace(; ylim = (-Float64(1), Float64(2)))
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
    Writers.write_field!(xyboxwriter, xyboxfield, xyboxdiagnostic, xyboxu, p, t)

    NCDatasets.NCDataset(joinpath(output_dir, "my_short_name_xybox.nc")) do nc
        # Test dimensions
        hpts, vpts = Writers.target_coordinates(
            xyboxspace,
            (NUM, 2NUM, 3NUM),
            ClimaDiagnostics.Writers.LevelsMethod(),
        )
        lon, lat = hpts
        @test nc["x"][:] == lon
        @test nc["y"][:] == lat
        @test nc["z"][:] == vpts
    end

    longlatboxspace =
        BoxSpace(; lonlat = true, ylim = (-Float64(1), Float64(2)))
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

    NCDatasets.NCDataset(
        joinpath(output_dir, "my_short_name_longlatbox.nc"),
    ) do nc
        # Test dimensions
        hpts, vpts = Writers.target_coordinates(
            longlatboxspace,
            (NUM, 2NUM, 3NUM),
            longlatboxwriter.z_sampling_method,
        )

        lon, lat = hpts
        @test nc["lon"][:] == lon
        @test nc["lat"][:] == lat
        @test nc["z"][:] == vpts
    end

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
            NCDatasets.NCDataset(
                joinpath(output_dir, "my_short_name_c$(i).nc"),
            ) do nc
                # Test dimensions
                vpts = Writers.target_coordinates(
                    colspace,
                    (NUM,),
                    colwriter.z_sampling_method,
                )

                @test nc["z"][:] == vpts
            end
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
    point_writer.preallocated_output_arrays[point_diagnostic] = [point_val]
    # No interpolation needed for point space
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

    NCDatasets.NCDataset(joinpath(output_dir, "my_short_name_point.nc")) do nc
        @test nc["ABC"][:] == [point_val, point_val]
    end

    ###################
    # Horizontal Space#
    ###################

    horizontal_space = ClimaCore.Spaces.level(space, 1)
    horizontal_field = Fields.coordinate_field(horizontal_space).z
    horizontal_writer = Writers.NetCDFWriter(
        horizontal_space,
        output_dir;
        num_points = (NUM, 2NUM),
    )
    horizontal_u = (; field = horizontal_field)

    horizontal_diagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
        variable = ClimaDiagnostics.DiagnosticVariable(;
            compute!,
            short_name = "ABC",
        ),
        output_short_name = "my_short_name_horizontal",
        output_long_name = "My Long Name Point Horizontal",
        output_writer = horizontal_writer,
    )

    Writers.interpolate_field!(
        horizontal_writer,
        horizontal_field,
        horizontal_diagnostic,
        horizontal_u,
        p,
        t,
    )
    Writers.write_field!(
        horizontal_writer,
        horizontal_field,
        horizontal_diagnostic,
        horizontal_u,
        p,
        t,
    )
    # Write a second time
    Writers.write_field!(
        horizontal_writer,
        horizontal_field,
        horizontal_diagnostic,
        horizontal_u,
        p,
        t,
    )
    close(horizontal_writer)
    NCDatasets.NCDataset(
        joinpath(output_dir, "my_short_name_horizontal.nc"),
    ) do nc
        @test size(nc["ABC"]) == (2, NUM, 2NUM)

        lon, lat = Writers.target_coordinates(horizontal_space, (NUM, 2NUM))
        @test nc["lon"][:] == lon
        @test nc["lat"][:] == lat
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

@testset "NetCDFWriter write field time test" begin
    space = SphericalShellSpace(FT = Float32)
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

    u = (; field)
    # FIXME: We are hardcoding the start date
    p = (; start_date = Dates.DateTime(2010, 1))
    # `ITime` should save the times as `Float64`
    # Float32(1000 * 20995200.0) |> Dates.Millisecond is not equal to
    # 20_995_200_000 milliseconds, but it is true for Float64
    t = ITime(
        20_995_200,
        period = Dates.Second(1),
        epoch = Dates.DateTime(2010, 1),
    )

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
        output_short_name = "timetest",
        output_long_name = "My Long Name",
        output_writer = writer,
    )
    Writers.interpolate_field!(writer, field, diagnostic, u, p, t)
    Writers.write_field!(writer, field, diagnostic, u, p, t)

    # This is div(2^53, 1000), since 2^53 + 1 is the first integer that cannot be
    # represented exactly as a Float64 due to rounding error
    # Divide by 1000 because time conversion function will go through
    # milliseconds (e.g. see the conversion from seconds to dates in
    # write_field!)
    t = ITime(
        div(2^53, 1000),
        period = Dates.Second(1),
        epoch = Dates.DateTime(2010, 1),
    )
    Writers.interpolate_field!(writer, field, diagnostic, u, p, t)
    Writers.write_field!(writer, field, diagnostic, u, p, t)

    NCDatasets.NCDataset(joinpath(output_dir, "timetest.nc")) do nc
        times = nc["time"][:]
        @test eltype(times) == Float64
        @test Dates.Second(Dates.Millisecond(round(1000 * times[1]))) ==
              Dates.Second(20_995_200)

        @test Dates.Second(Dates.Millisecond(round(1000 * times[2]))) ==
              Dates.Second(div(2^53, 1000))
    end
end

@testset "NetCDFWriter time reductions behavior" begin
    space = SphericalShellSpace(FT = Float32)
    field = Fields.coordinate_field(space).z

    NUM = 10

    start_date = Dates.DateTime(2010, 1, 1)

    function compute!(out, u, p, t)
        if isnothing(out)
            return u.field
        else
            out .= u.field
        end
    end

    u = (; field)
    p = (; start_date = start_date)

    @testset "Instantaneous diagnostics" begin
        writer = Writers.NetCDFWriter(
            space,
            output_dir;
            num_points = (NUM, 2NUM, 3NUM),
            start_date = start_date,
        )

        diagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
            variable = ClimaDiagnostics.DiagnosticVariable(;
                compute!,
                short_name = "INST",
            ),
            output_short_name = "instant_test",
            output_long_name = "Instantaneous Test",
            output_writer = writer,
        )

        t1 = 10.0
        t2 = 20.0
        t3 = 30.0
        for t in [t1, t2, t3]
            Writers.interpolate_field!(writer, field, diagnostic, u, p, t)
            Writers.write_field!(writer, field, diagnostic, u, p, t)
        end

        NCDatasets.NCDataset(joinpath(output_dir, "instant_test.nc")) do nc
            @test nc["time"] == [t1, t2, t3]

            dates = [
                start_date + Dates.Millisecond(round(1000 * t1)),
                start_date + Dates.Millisecond(round(1000 * t2)),
                start_date + Dates.Millisecond(round(1000 * t3)),
            ]
            @test nc["date"] == dates

            @test nc["time_bnds"][:, 1] == [0.0; t1]
            @test nc["time_bnds"][:, 2] == [t1; t2]
            @test nc["time_bnds"][:, 3] == [t2; t3]

            @test nc["date_bnds"][:, 1] ==
                  [start_date, start_date + Dates.Second(t1)]
            @test nc["date_bnds"][:, 2] ==
                  [start_date + Dates.Second(t1), start_date + Dates.Second(t2)]
            @test nc["date_bnds"][:, 3] ==
                  [start_date + Dates.Second(t2), start_date + Dates.Second(t3)]
        end

        close(writer)
    end

    @testset "Reduced diagnostics (average)" begin
        writer = Writers.NetCDFWriter(
            space,
            output_dir;
            num_points = (NUM, 2NUM, 3NUM),
            start_date = start_date,
        )

        diagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
            variable = ClimaDiagnostics.DiagnosticVariable(;
                compute!,
                short_name = "REDUCED",
            ),
            output_short_name = "reduced_test",
            output_long_name = "Reduced Test",
            output_writer = writer,
            reduction_time_func = (+),
        )

        t1 = 10.0
        t2 = 20.0
        t3 = 30.0
        for t in [t1, t2, t3]
            Writers.interpolate_field!(writer, field, diagnostic, u, p, t)
            Writers.write_field!(writer, field, diagnostic, u, p, t)
        end

        NCDatasets.NCDataset(joinpath(output_dir, "reduced_test.nc")) do nc
            @test nc["time"] == [0.0, t1, t2]

            @test nc["date"][1] == start_date
            @test nc["date"][2] == start_date + Dates.Second(t1)
            @test nc["date"][3] == start_date + Dates.Second(t2)

            @test nc["time_bnds"][:, 1] == [0.0; t1]
            @test nc["time_bnds"][:, 2] == [t1; t2]
            @test nc["time_bnds"][:, 3] == [t2; t3]

            @test nc["date_bnds"][:, 1] ==
                  [start_date, start_date + Dates.Second(t1)]
            @test nc["date_bnds"][:, 2] ==
                  [start_date + Dates.Second(t1), start_date + Dates.Second(t2)]
            @test nc["date_bnds"][:, 3] ==
                  [start_date + Dates.Second(t2), start_date + Dates.Second(t3)]
        end

        close(writer)
    end

    @testset "Calendar-based reductions (monthly)" begin
        writer = Writers.NetCDFWriter(
            space,
            output_dir;
            num_points = (NUM, 2NUM, 3NUM),
            start_date = start_date,
        )

        # Test monthly to verify variable-length periods work correctly
        # (Jan=31 days, Feb=28 days). This covers the general case since
        # the code doesn't distinguish between different calendar periods.
        diagnostic = ClimaDiagnostics.ScheduledDiagnostic(;
            variable = ClimaDiagnostics.DiagnosticVariable(;
                compute!,
                short_name = "monthly",
            ),
            output_writer = writer,
            reduction_time_func = (+),
            output_schedule_func = ClimaDiagnostics.Schedules.EveryCalendarDtSchedule(
                Dates.Month(1);
                start_date,
            ),
        )

        times = [31, 59] .* 86400.0  # End of Jan and Feb in seconds
        for t in times
            Writers.interpolate_field!(writer, field, diagnostic, u, p, t)
            Writers.write_field!(writer, field, diagnostic, u, p, t)
        end

        NCDatasets.NCDataset(joinpath(output_dir, "monthly_1M_+.nc")) do nc
            @test nc["date"][:] ==
                  [Dates.DateTime(2010, 1, 1), Dates.DateTime(2010, 2, 1)]
            @test nc["date_bnds"][1, :] ==
                  [Dates.DateTime(2010, 1, 1), Dates.DateTime(2010, 2, 1)]
            @test nc["date_bnds"][2, :] ==
                  [Dates.DateTime(2010, 2, 1), Dates.DateTime(2010, 3, 1)]
            @test nc["time_bnds"][1, :] == [0.0, first(times)]
            @test nc["time_bnds"][2, :] == times
        end

        close(writer)
    end
end
