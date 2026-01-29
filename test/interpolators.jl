using Test

import ClimaCore.Fields
import ClimaCore.Remapping
import ClimaDiagnostics: Interpolators
import Dates

import ClimaInterpolations

include("TestTools.jl")

@testset "PfullInterpolator" begin
    shell_space = SphericalShellSpace(FT = Float32)
    col_space = ColumnCenterFiniteDifferenceSpace(FT = Float32)

    mutable struct FakePfullCompute
        """Keep track of the number of computes that are done"""
        num_compute::Int
    end

    function (f::FakePfullCompute)(out, u, p, t)
        f.num_compute += 1
        if isnothing(out)
            return t .* u.field
        else
            out .= t .* u.field
            return nothing
        end
    end

    for space in (shell_space, col_space)
        u = Fields.FieldVector(; field = ones(space))
        p = (;)
        t = 0

        fake_compute_pfull! = FakePfullCompute(0)
        pfull_intp = Interpolators.PfullInterpolator(
            fake_compute_pfull!,
            u,
            p,
            t;
            pfull_levels = Interpolators.era5_pressure_levels(),
        )

        # Call update! three times
        # The first and third update! tests that the update! change the pressure
        # field and the second update! tests if the caching works.
        Interpolators.update!(pfull_intp, u, p, 1.0)
        # Compute function is called twice when constructing and updating
        @test fake_compute_pfull!.num_compute == 2
        @test all(
            isone,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == 1.0

        Interpolators.update!(pfull_intp, u, p, 1.0)
        @test fake_compute_pfull!.num_compute == 2 # because of cached result
        @test all(
            isone,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == 1.0

        Interpolators.update!(pfull_intp, u, p, 2.0)
        @test fake_compute_pfull!.num_compute == 3
        @test all(
            x -> x == 2.0,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == 2.0

        Interpolators.force_update!(pfull_intp, u, p, 0.0)
        @test fake_compute_pfull!.num_compute == 4 # because of cached result
        @test all(
            iszero,
            parent(Remapping.pfull_field(pfull_intp.pressure_intp)),
        )
        @test pfull_intp.last_t[] == 0.0

        # interpolate_to_pfull_coords! and interpolate_to_pfull_coords!! gives the
        # same result and it is defined on space that we want
        field = copy(u.field)
        dest1 = Interpolators.interpolate_to_pfull_coords!(field, pfull_intp)
        field = copy(u.field)
        dest2 = zeros(axes(dest1))
        Interpolators.interpolate_to_pfull_coords!!(dest2, field, pfull_intp)
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
