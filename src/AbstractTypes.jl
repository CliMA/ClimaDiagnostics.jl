# Abstract types needed across submodules

"""
    AbstractWriter

An object that knows how to save some output.

`AbstractWriter`s have to provide one function, `write_field!`

The function has to have singature
`write_field!(writer::Writer, field, diagnostic::ScheduledDiagnostic, u, p, t)`

It is up to the writer to implment this
"""
abstract type AbstractWriter end
