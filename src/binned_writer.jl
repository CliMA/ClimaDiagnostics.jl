import StaticArrays.SVector
import ..DiagnosticVariables: DiagnosticVariable
import ..ScheduledDiagnostics: ScheduledDiagnostic

"""
    BinnedWriter(base_writer)

A writer for binned diagnostics. This is a wrapper around a base writer
(e.g., `NetCDFWriter`) that handles the special structure of binned data.
It writes bin counts as a multi-dimensional variable and bin edges as a
coordinate.
"""
struct BinnedWriter{W <: AbstractWriter} <: AbstractWriter
    base_writer::W
end

function interpolate_field!(writer::BinnedWriter, field, diagnostic, u, p, t)
    # The BinnedWriter does not need to interpolate the fields. The interpolation will be
    # handled by the base_writer after the binned field has been decomposed into scalars.
    return nothing
end

function write_field!(
    writer::BinnedWriter,
    field, # field of BinnedAccumulators
    diagnostic,
    u,
    p,
    t,
)
    space = axes(field)
    ncols = ClimaCore.Spaces.ncolumns(space)
    nlevs = ClimaCore.Spaces.nlevels(space)

    if ncols == 1 && nlevs == 1
        # Single point/column - access bin_edges directly from the field element
        bin_edges = getproperty(field, :bin_edges)[]
    elseif ncols > 1 && nlevs == 1
        # 2D spectral element space - extract first point
        if space isa Spaces.SpectralElementSpace2D
            single_point_field = ClimaCore.column(field, 1, 1, 1)
        elseif space isa Spaces.SpectralElementSpace1D
            single_point_field = ClimaCore.column(field, 1, 1)
        else
            error("Unexpected 2D space type: $(typeof(space))")
        end
        bin_edges = getproperty(single_point_field, :bin_edges)[]
    elseif nlevs > 1
        # Extruded space - extract first column, then first level
        single_column_field = ClimaCore.column(field, 1, 1, 1)
        # Extract the first level from the column
        single_level_field = ClimaCore.Fields.level(single_column_field, 1)
        bin_edges = getproperty(single_level_field, :bin_edges)[]
    else
        error("Unexpected space dimensions: ncols=$ncols, nlevs=$nlevs")
    end

    n_bins = length(bin_edges) - 1

    # We will write one variable for each bin count. This is simpler than creating
    # a new dimension in the output file, which would require more invasive
    # changes to the writer infrastructure.
    bin_counts_field = getproperty(field, :bin_counts)

    for i in 1:n_bins
        # For each element in the `bin_counts_field` (which are SVectors),
        # take the i-th component. This returns a new field of scalars.
        bin_i_counts = getindex.(bin_counts_field, i)

        # Create a descriptive name and attributes for this bin's variable
        var = diagnostic.variable
        bin_name = "$(var.short_name)_bin_$(i)"
        bin_long_name = "$(diagnostic.output_long_name) (bin $(i))"
        bin_variable = DiagnosticVariable(;
            short_name = bin_name,
            long_name = var.long_name,
            units = "count",
            comments = "bin_range = [$(bin_edges[i]), $(bin_edges[i+1]))",
            compute! = var.compute!,
        )
        bin_diagnostic = ScheduledDiagnostic(;
            variable = bin_variable,
            output_schedule_func = diagnostic.output_schedule_func,
            output_writer = writer.base_writer,
            reduction_time_func = diagnostic.reduction_time_func,
            compute_schedule_func = diagnostic.compute_schedule_func,
            pre_output_hook! = diagnostic.pre_output_hook!,
            output_short_name = bin_name,
            output_long_name = bin_long_name,
        )

        # Interpolate first, so the base_writer can populate its internal state. Then,
        # use the base writer to write this scalar field.
        interpolate_field!(
            writer.base_writer,
            bin_i_counts,
            bin_diagnostic,
            u,
            p,
            t,
        )
        write_field!(writer.base_writer, bin_i_counts, bin_diagnostic, u, p, t)
    end
end
