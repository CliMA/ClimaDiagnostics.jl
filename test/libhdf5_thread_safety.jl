#=
Exercise the exact code pattern that would SIGSEGV without the
`LIBHDF5_LOCK` fix in v0.3.4: concurrent NCDatasets writes through
`NetCDFWriter` interleaved with HDF5.jl writes, driven from multiple
Julia threads.

HDF5_jll is built without `--enable-threadsafe` (`Threadsafety: no`
in `libhdf5.settings`). Before v0.3.4, `NCDatasets.NETCDF_LOCK` and
`HDF5.API.liblock` were independent, so two Julia threads could enter
libhdf5 simultaneously and corrupt its internal skip-list /
property-list state, yielding:

    NetCDF: HDF error (NetCDF error code: -101)
    double free or corruption (fasttop) + SIGSEGV in H5SL_insert

This test writes through `NetCDFWriter` from `nthreads()` worker
threads while a side task spins an HDF5.jl write loop, for
`NCYCLES` cycles. If the lock is wired correctly, both paths
serialise on `HDF5.API.liblock` and the test exits cleanly. If the
lock regresses (e.g. someone points `LIBHDF5_LOCK` at a private
`ReentrantLock()` again), this test SIGSEGVs or raises
`NetCDFError(-101)` — either outcome is fatal to the process and
visible as a failed CI step.

The MT race test is skipped when `Threads.nthreads() < 2` — single-
threaded Julia cannot reproduce the race. Dedicated Buildkite steps
`cpu_mt_tests` and `gpu_tests` (both launched with `--threads=8`)
exercise it.
=#
using Test
import Dates
import NCDatasets
import HDF5
using Base.Threads: @threads, Atomic, atomic_add!

import ClimaCore
import ClimaCore.Fields
import ClimaCore.Spaces
import ClimaCore.InputOutput
import ClimaComms

import ClimaDiagnostics
import ClimaDiagnostics.Writers

include("TestTools.jl")

const NCYCLES = 40
const OUTDIR = mktempdir(pwd(); prefix = "libhdf5_mt_")

@info "libhdf5 thread-safety MT test" nthreads = Threads.nthreads() OUTDIR

# Invariant: LIBHDF5_LOCK must equal HDF5.API.liblock. This is the
# "shared barrier" the fix depends on.
@testset "LIBHDF5_LOCK is shared with HDF5.jl" begin
    @test isdefined(Writers, :LIBHDF5_LOCK)
    @test Writers.LIBHDF5_LOCK isa ReentrantLock
    @test Writers.LIBHDF5_LOCK === HDF5.API.liblock
end

if Threads.nthreads() < 2
    @info "Skipping MT race test — Threads.nthreads() < 2"
else
    @testset "Concurrent NetCDFWriter + HDF5.jl writes do not corrupt libhdf5" begin
        # Tiny BoxSpace — just enough to construct a real NetCDFWriter.
        space = BoxSpace(; nelements = (3, 3), zelem = 4, npolynomial = 3)
        FT = ClimaCore.Spaces.undertype(space)

        writer = Writers.NetCDFWriter(space, OUTDIR;
            num_points = (4, 4, 4),
            start_date = Dates.DateTime(2020, 1, 1),
        )

        # Build NFIELDS real ScheduledDiagnostic + Field pairs.
        NFIELDS = max(2, Threads.nthreads())
        fields = [Fields.Field(FT, space) for _ in 1:NFIELDS]
        for (i, f) in enumerate(fields)
            parent(f) .= FT(i)
        end

        # Simple compute! that copies the pre-filled field.
        diagnostics = map(1:NFIELDS) do i
            pinned_field = fields[i]
            compute_i!(out, u, p, t) = (out .= pinned_field; out)
            ClimaDiagnostics.ScheduledDiagnostic(;
                variable = ClimaDiagnostics.DiagnosticVariable(;
                    compute! = compute_i!,
                    short_name = "var$(lpad(i, 2, '0'))",
                ),
                output_short_name = "var$(lpad(i, 2, '0'))_short",
                output_long_name = "Field $i",
                output_writer = writer,
            )
        end

        # Prime the remappers / prealloc buffers (first call allocates).
        for i in 1:NFIELDS
            Writers.interpolate_field!(
                writer, fields[i], diagnostics[i], nothing, nothing, 0.0,
            )
        end

        # Concurrent ClimaCore.InputOutput.HDF5Writer checkpoint loop on
        # a side task — this is the exact writer the production
        # `nan_checking_callback` uses, so we race NetCDFWriter against
        # the real ClimaCore HDF5 path, not just raw `HDF5.h5open`.
        context = ClimaComms.context(fields[1])
        h5_stop = Atomic{Bool}(false)
        h5_writes = Atomic{Int}(0)
        h5_task = Threads.@spawn begin
            k = 0
            while !h5_stop[]
                k += 1
                path = joinpath(OUTDIR, "ckpt_$(lpad(k, 4, '0')).hdf5")
                InputOutput.HDF5Writer(path, context) do hdfwriter
                    InputOutput.write!(hdfwriter, fields[1], "state")
                end
                atomic_add!(h5_writes, 1)
                yield()
            end
        end

        nc_writes = Atomic{Int}(0)
        try
            for cycle in 1:NCYCLES
                @threads for i in 1:NFIELDS
                    Writers.write_field!(
                        writer, fields[i], diagnostics[i],
                        nothing, nothing, Float64(cycle),
                    )
                    atomic_add!(nc_writes, 1)
                end
                Writers.sync(writer)
            end
            @test nc_writes[] == NCYCLES * NFIELDS
            @info "Survived $NCYCLES cycles × $NFIELDS threads. \
                   NC writes=$(nc_writes[]), HDF5.jl writes=$(h5_writes[])"
        finally
            h5_stop[] = true
            wait(h5_task)
            Base.close(writer)
        end
    end
end
