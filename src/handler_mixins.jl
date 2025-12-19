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
`PfullCoordsStyle`, but these fields appear in the struct for convenience.
"""
struct PfullMixin{
    F <: Function,
    PRESSURE, # TODO: Add type annotation here!
    FIELDS, # TODO: Add type annotation here!
    PERM_MATRIX <: AbstractMatrix,
} <: AbstractHandlerMixin
    """Container storing the fields from the compute functions."""
    compute_fields::FIELDS

    """Function to compute the pressure field"""
    pfull_compute!::F # TODO: Add check that this is the same across all coords style

    """ClimaCore field of pressure created by pfull_compute!"""
    pfull_field::PRESSURE # TODO: Add check that this is the same across all coords style

    """A permutation matrix, created by sortperm, for
    sorting the pressures for each column"""
    perm_matrix::PERM_MATRIX # TODO: This is the same across all pressure_coords_style by virtue of pfull_field
end
