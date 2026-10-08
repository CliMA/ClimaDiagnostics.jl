"""
The `Writers` module defines generic writers for diagnostics with opinionated defaults.

Currently, it implements:
- `DummyWriter`, does nothing, useful for testing ClimaDiagnostics
- `DictWriter`, for in-memory, dictionary-based writing
- `HDF5Writer`, to save raw `ClimaCore` `Fields` to HDF5 files
- `NetCDFWriter`, to save remapped `ClimaCore` `Fields` to NetCDF files

Writers can implement:
- `interpolate_field!`
- `write_field!`
- `Base.close`
- `sync`
- `z_sampling_method`
"""
module Writers

import ClimaCore

import ..AbstractWriter, ..ScheduledDiagnostic
import ..ScheduledDiagnostics: output_short_name, output_long_name

function write_field!(writer::AbstractWriter, field, diagnostic, u, p, t)
    # Nothing to be done here :)
    return nothing
end

function Base.close(writer::AbstractWriter)
    # Nothing to be done here :)
    return nothing
end

function interpolate_field!(writer::AbstractWriter, field, diagnostic, u, p, t)
    # Nothing to be done here :)
    return nothing
end

function sync(writer::AbstractWriter)
    # Nothing to be done here :)
    return nothing
end

"""
    z_sampling_method(writer::AbstractWriter)

Return the `AbstractZSamplingMethod` used by `writer` to sample the vertical
direction, or `nothing` if `writer` does not interpolate along the vertical
(the default).
"""
function z_sampling_method(writer::AbstractWriter)
    return nothing
end

include("dummy_writer.jl")
include("dict_writer.jl")
include("hdf5_writer.jl")
include("netcdf_writer.jl")

end
