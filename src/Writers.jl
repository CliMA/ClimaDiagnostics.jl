# Writers.jl
#
# This file defines generic writers for diagnostics with opinionated defaults.

module Writers

import ..AbstractWriter, ..ScheduledDiagnostic
import ..ScheduledDiagnostics: output_short_name

include("dict_writer.jl")
include("hdf5_writer.jl")
include("netcdf_writer.jl")

end
