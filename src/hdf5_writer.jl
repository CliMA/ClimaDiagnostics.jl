import ClimaComms
import ClimaCore.InputOutput

"""
    HDF5Writer()

Save a `ScheduledDiagnostic` to a HDF5 file inside the `output_dir` of the simulation.

TODO: This is a very barebone HDF5Writer!

We need to implement the following features/options:
- Toggle for write new files/append
- Checks for existing files
- Check for new subfolders that have to be created
- More meaningful naming conventions (keeping in mind that we can have multiple variables
  with different reductions)
- All variables in one file/each variable in its own file
- All timesteps in one file/each timestep in its own file
- Writing the correct attributes
- Overriding simulation.output_dir (e.g., if the path starts with /)
- ...more features/options

"""
struct HDF5Writer <: AbstractWriter
    output_dir::String
end

"""
    close(writer::HDF5Writer)

Close all the files open in `writer`. (Currently no-op.)
"""
Base.close(writer::HDF5Writer) = nothing

function write_field!(writer::HDF5Writer, field, diagnostic, u, p, t)
    var = diagnostic.variable
    time = t

    output_path = joinpath(
        writer.output_dir,
        "$(diagnostic.output_short_name)_$(time).h5",
    )

    comms_ctx = ClimaComms.context(u.c)
    hdfwriter = InputOutput.HDF5Writer(output_path, comms_ctx)
    InputOutput.write!(hdfwriter, field, "$(diagnostic.output_short_name)")
    attributes = Dict(
        "time" => time,
        "long_name" => diagnostic.output_long_name,
        "variable_units" => var.units,
        "standard_variable_name" => var.standard_name,
    )

    for (k, v) in attributes
        InputOutput.HDF5.write_attribute(hdfwriter.file, k, v)
    end

    Base.close(hdfwriter)
    return nothing
end
