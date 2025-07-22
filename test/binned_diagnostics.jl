using Test
import ClimaDiagnostics
import ClimaDiagnostics.Schedules: EveryDtSchedule
import ClimaDiagnostics.Writers: DictWriter, BinnedWriter
import ClimaCore
import ClimaCore.Fields
import ClimaCore.Geometry
import ClimaCore.Spaces
import NCDatasets
using StaticArrays

@testset "Binned Diagnostics" begin

    @testset "BinnedAccumulator" begin
        bins = [0.0, 1.0, 2.0, 5.0]
        acc = ClimaDiagnostics.BinnedAccumulator(bins)

        # Correct type and parameters
        @test acc isa ClimaDiagnostics.BinnedAccumulator{Float64, 3, 4}

        # Bin edges and counts are SVectors of the same type
        @test acc.bin_edges isa SVector{4, Float64}
        @test acc.bin_counts isa SVector{3, Float64}
        @test acc.bin_edges == SVector{4}(bins)
        @test all(acc.bin_counts .== 0.0)
        @test acc.total_count == 0.0

        # Copy constructor
        acc_copy = copy(acc)
        @test acc_copy.bin_edges == acc.bin_edges
        @test acc_copy.bin_counts == acc.bin_counts
        @test acc_copy.total_count == acc.total_count

        # Test show
        io = IOBuffer()
        show(io, acc)
        str = String(take!(io))
        @test occursin("BinnedAccumulator", str)
        @test occursin("3 bins", str)
        @test occursin("0 values", str)
    end

    @testset "bin_index" begin
        bins = SVector(0.0, 1.0, 2.0, 5.0)

        @testset "Values within bins" begin
            # Test values strictly within bin ranges
            @test ClimaDiagnostics.bin_index(bins, 0.5) == 1
            @test ClimaDiagnostics.bin_index(bins, 1.5) == 2
            @test ClimaDiagnostics.bin_index(bins, 3.0) == 3
        end

        @testset "Values on bin edges" begin
            # Test values that fall exactly on bin edges
            # searchsortedlast makes bins left-inclusive: [left, right)
            @test ClimaDiagnostics.bin_index(bins, 0.0) == 1
            @test ClimaDiagnostics.bin_index(bins, 1.0) == 2
            @test ClimaDiagnostics.bin_index(bins, 2.0) == 3
        end

        @testset "Values outside range" begin
            # Test values outside the overall bin range
            @test isnothing(ClimaDiagnostics.bin_index(bins, -1.0))
            @test isnothing(ClimaDiagnostics.bin_index(bins, 5.0)) # right edge is exclusive
            @test isnothing(ClimaDiagnostics.bin_index(bins, 6.0))
        end
    end

    @testset "binned_reduction" begin
        bins = [0.0, 1.0, 2.0, 5.0]
        acc = ClimaDiagnostics.BinnedAccumulator(bins)

        # Test adding values to different bins
        acc1 = ClimaDiagnostics.binned_reduction(acc, 0.5)  # bin 1
        @test acc1.bin_counts[1] == 1
        @test acc1.bin_counts[2] == 0
        @test acc1.bin_counts[3] == 0
        @test acc1.total_count == 1

        acc2 = ClimaDiagnostics.binned_reduction(acc1, 1.5)  # bin 2
        @test acc2.bin_counts[1] == 1
        @test acc2.bin_counts[2] == 1
        @test acc2.bin_counts[3] == 0
        @test acc2.total_count == 2

        acc3 = ClimaDiagnostics.binned_reduction(acc2, 3.0)  # bin 3
        @test acc3.bin_counts[1] == 1
        @test acc3.bin_counts[2] == 1
        @test acc3.bin_counts[3] == 1
        @test acc3.total_count == 3

        # Test adding multiple values to same bin
        acc4 = ClimaDiagnostics.binned_reduction(acc3, 0.7)  # bin 1 again
        @test acc4.bin_counts[1] == 2
        @test acc4.bin_counts[2] == 1
        @test acc4.bin_counts[3] == 1
        @test acc4.total_count == 4

        # Test values outside range (should not change accumulator)
        acc_outside = ClimaDiagnostics.binned_reduction(acc4, -1.0)
        @test acc_outside.bin_counts == acc4.bin_counts
        @test acc_outside.total_count == acc4.total_count

        acc_outside2 = ClimaDiagnostics.binned_reduction(acc_outside, 10.0)
        @test acc_outside2.bin_counts == acc4.bin_counts
        @test acc_outside2.total_count == acc4.total_count
    end

    @testset "create_binned_diagnostic" begin
        # Dummy data and parameters
        var = ClimaDiagnostics.DiagnosticVariable(
            short_name = "test_var",
            long_name = "Test Variable",
            units = "none",
            compute = (state) -> 1.0,
        )
        bins = [0.0, 1.0, 2.0, 5.0]
        writer = ClimaDiagnostics.Writers.DummyWriter()

        # Create the diagnostic
        binned_diag = ClimaDiagnostics.create_binned_diagnostic(
            var,
            bins,
            ClimaDiagnostics.EveryDtSchedule(1),
            ClimaDiagnostics.EveryDtSchedule(1),
            writer,
        )

        # Basic validation
        @test binned_diag isa ClimaDiagnostics.ScheduledDiagnostic
        @test binned_diag.variable.short_name == "test_var_binned"
        @test binned_diag.output_writer isa
              ClimaDiagnostics.Writers.BinnedWriter

        # Test the reducer
        @test binned_diag.reduction_time_func isa ClimaDiagnostics.BinnedReducer

        # Test the identity
        identity_accumulator = ClimaDiagnostics.identity_of_reduction(
            binned_diag.reduction_time_func,
        )
        @test identity_accumulator.bin_edges == bins
        @test all(identity_accumulator.bin_counts .== 0)
    end

    @testset "pre_output_hooks" begin
        # Create a test field of accumulators
        bins = [0.0, 1.0, 2.0, 5.0]
        initial_acc = ClimaDiagnostics.BinnedAccumulator(bins)

        # Using a Vector instead of a Field to simplify testing
        test_field = [deepcopy(initial_acc) for _ in 1:2]

        # Simulate a reduction to populate the accumulator
        test_field[1] = ClimaDiagnostics.binned_reduction(test_field[1], 0.5)
        test_field[2] = ClimaDiagnostics.binned_reduction(test_field[2], 1.5)

        @test test_field[1].total_count == 1.0
        @test test_field[1].bin_counts == [1.0, 0.0, 0.0]

        # Apply the hook
        ClimaDiagnostics.binned_pre_output_hook!(test_field, 1)

        # Check if the hook reset the accumulators
        @test test_field[1].total_count == 0.0
        @test all(test_field[1].bin_counts .== 0.0)
        @test test_field[2].total_count == 0.0
        @test all(test_field[2].bin_counts .== 0.0)
    end
end
