"""The `DummyWriter` does nothing. It can be used to temporarily disable output of a
scheduled diagnostic. Mostly useful for testing and debugging ClimaDiagnostics.
"""
struct DummyWriter <: AbstractWriter end

function write_field!(writer::DummyWriter, field, diagnostic, u, p, t)
    # Nothing to be done here :)
    # return nothing
end

function Base.close(writer::DummyWriter)
    # Nothing to be done here :)
    return nothing
end

function interpolate_field!(writer::DummyWriter, field, diagnostic, u, p, t)
    # Nothing to be done here :)
    return nothing
end
