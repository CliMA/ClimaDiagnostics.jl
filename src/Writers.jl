# Writers.jl
#
# This file defines generic writers for diagnostics with opinionated defaults.

module Writers

include("hdf5_writer.jl")
include("netcdf_writer.jl")

end
