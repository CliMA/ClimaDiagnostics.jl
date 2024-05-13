import ClimaComms
import ClimaCore.InputOutput

"""
    HDF5Writer(output_dir)

Save a `ScheduledDiagnostic` to a HDF5 file inside the `output_dir`.

TODO: This is a very barebone HDF5Writer!
"""
struct HDF5Writer <: AbstractWriter
    output_dir::String
end

"""
    close(writer::HDF5Writer)

Close all the files open in `writer`. (Currently no-op.)
"""
function Base.close(writer::HDF5Writer)
    return nothing
end

"""
    HDF5Writer(output_dir)

Save a `ScheduledDiagnostic` to a HDF5 file inside the `output_dir`.

The name of the file is determined by the `output_short_name` of the output
`ScheduledDiagnostic`. New files are created for each timestep.

`Field`s can be read back using the `InputOutput` module in `ClimaCore`.
"""
function write_field!(writer::HDF5Writer, field, diagnostic, u, p, t)
    var = diagnostic.variable
    time = t

    output_path = joinpath(
        writer.output_dir,
        "$(diagnostic.output_short_name)_$(time).h5",
    )

    comms_ctx = ClimaComms.context(field)
    hdfwriter = InputOutput.HDF5Writer(output_path, comms_ctx)
    InputOutput.write!(hdfwriter, field, "$(diagnostic.output_short_name)")
    attributes = Dict(
        "time" => time,
        "short_name" => diagnostic.output_short_name,
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
