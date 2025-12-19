import ClimaCore: Fields

abstract type AbstractHandlerMixin end

"""
    NoMixIn <: AbstractHandlerMixin

No extra functionality is added to `DiagnosticsHandler`.
"""
struct NoMixin <: AbstractHandlerMixin end

"""
    PfullMixin <: AbstractHandlerMixin

This is a companion struct to `PfullCoordsStyle`. Extra functionality, in
particular `compute_fields`, is added to complement `PfullCoordsStyle`.

Note that `pfull_compute!`, `pfull_field`, and `perm_matrix` are stored in
`PfullCoordsStyle`, but these fields appear in this struct for convenience.
"""
struct PfullMixin{
    F <: Function,
    PRESSURE <: Fields.Field,
    FIELDS, # TODO: Add type annotation here!
    PERM_MATRIX <: AbstractMatrix,
} <: AbstractHandlerMixin
    """Container storing the fields from the compute functions."""
    compute_fields::FIELDS

    """Function to compute the pressure field"""
    pfull_compute!::F

    """ClimaCore field of pressure created by pfull_compute!"""
    pfull_field::PRESSURE

    """A permutation matrix, created by sortperm, for
    sorting the pressures for each column"""
    perm_matrix::PERM_MATRIX
end
