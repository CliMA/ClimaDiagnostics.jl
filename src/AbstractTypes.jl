# Abstract types needed across submodules

"""
    AbstractWriter

An object that knows how to save some output.

`AbstractWriter`s have to provide one function, `write_field!`

The function has to have signature
`write_field!(writer::Writer, field, diagnostic, u, p, t, window_start)`,
where `window_start` is the start of the output window that ends at `t`.

It is up to the `Writer` to implement this.
"""
abstract type AbstractWriter end
