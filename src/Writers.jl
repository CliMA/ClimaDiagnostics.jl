"""
The `Writers` module defines generic writers for diagnostics with opinionated defaults.

Currently, it implements:
- `DictWriter`, for in-memory, dictionary-based writing
- `HDF5Writer`, to save raw `ClimaCore` `Fields` to HDF5 files
- `NetCDFWriter`, to save remapped `ClimaCore` `Fields` to NetCDF files
"""
module Writers

import ..AbstractWriter, ..ScheduledDiagnostic
import ..ScheduledDiagnostics: output_short_name

include("dict_writer.jl")
include("hdf5_writer.jl")
include("netcdf_writer.jl")

end
