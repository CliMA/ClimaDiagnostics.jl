"""The `DictWriter` is a writer that does not write to disk, but to memory (in a
dictionary).

This is particularly useful for testing and debugging.

This is not type stable (the underlying dictionary does not know in advance what types might
be used)
"""
struct DictWriter{T} <: AbstractWriter
    dict::T
end

function DictWriter()
    return DictWriter(Dict())
end

function write_field!(writer::DictWriter, field, diagnostic, u, p, t)
    key_name =
        diagnostic isa ScheduledDiagnostic ? output_short_name(diagnostic) :
        diagnostic
    diagnostic_dict = get!(writer.dict, key_name, Dict())
    diagnostic_dict[t] = copy(field)
    return nothing
end

function Base.getindex(writer::DictWriter, key)
    return Base.getindex(writer.dict, key)
end
