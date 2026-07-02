using Test

import ClimaCore.Fields
import ClimaCore.Remapping
import ClimaDiagnostics: Interpolators
import Dates

import ClimaUtilities.TimeManager: ITime

include("TestTools.jl")

@testset "PressureInterpolator" begin
    shell_space = SphericalShellSpace(FT = Float32)
    col_space = ColumnCenterFiniteDifferenceSpace(FT = Float32)

    mutable struct MyTime
        t::Float32
    end

    Base.float(time::MyTime) = float(time.t)
    Base.:(==)(x::MyTime, y::MyTime) = x.t == y.t

    function compute_pressure_field(u, t)
        u.field .*= float(t)
    end

    for space in (shell_space, col_space)
        u = Fields.FieldVector(; field = ones(space))
        p = (;)
        t = MyTime(0.0)

        pfull_intp = Interpolators.PressureInterpolator(
            u.field,
            t;
            pressure_levels = Interpolators.era5_pressure_levels(),
        )

        # Call update! three times
        # The first and third update! tests that the update! change the pressure
        # field and the second update! tests if the caching works.
        t1 = MyTime(1.0)
        compute_pressure_field(u, t1)
        Interpolators.update!(pfull_intp, t1)
        # Compute function is called twice when constructing and updating
        @test all(
            isone,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == t1

        t2 = MyTime(1.0)
        compute_pressure_field(u, t2)
        Interpolators.update!(pfull_intp, t2)
        @test all(
            isone,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] === t1
        @test !(pfull_intp.last_t[] === t2)

        t3 = MyTime(2.0)
        compute_pressure_field(u, t3)
        Interpolators.update!(pfull_intp, t3)
        @test all(
            x -> x == 2.0,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == t3

        t4 = MyTime(0.0)
        compute_pressure_field(u, t4)
        Interpolators.force_update!(pfull_intp, t4)
        @test all(
            iszero,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == t4

        # interpolate_to_pressure_coords and interpolate_to_pressure_coords! gives the
        # same result and it is defined on space that we want
        field = copy(u.field)
        dest1 = Interpolators.interpolate_to_pressure_coords(field, pfull_intp)
        field = copy(u.field)
        dest2 = zeros(axes(dest1))
        Interpolators.interpolate_to_pressure_coords!(dest2, field, pfull_intp)
        @test isequal(dest1, dest2)

        # Test coordinate field is what we expect
        @test isequal(
            Array(
                Fields.field2array(Fields.coordinate_field(axes(dest1)).p)[
                    :,
                    1,
                ],
            ),
            Interpolators.era5_pressure_levels(),
        )
    end
end

@testset "PressureInterpolator show" begin
    space = ColumnCenterFiniteDifferenceSpace(FT = Float32)
    u = Fields.FieldVector(; field = ones(space))
    pfull_intp = Interpolators.PressureInterpolator(u.field, 0.0)

    out = sprint(show, MIME("text/plain"), pfull_intp)
    @test occursin("PressureInterpolator", out)
    @test count(==('\n'), out) <= 10

    out2 = sprint(show, pfull_intp)
    @test occursin("PressureInterpolator", out2)
    @test !occursin('\n', out2)

    out3 =
        sprint(show, MIME("text/plain"), pfull_intp; context = :compact => true)
    @test out2 == out3

    out_summary = sprint(summary, pfull_intp)
    @test occursin("PressureInterpolator", out_summary)
    @test !occursin('\n', out_summary)
end
