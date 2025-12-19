abstract type AbstractHandlerMixin end

"""
    NoMixIn <: AbstractHandlerMixin

No extra functionality is added to `DiagnosticsHandler`.
"""
struct NoMixin <: AbstractHandlerMixin end

"""


This is a companion struct to `PfullCoordsStyle`. Extra functionality, in
particular `compute_fields`, is added to complement `PfullCoordsStyle`.
"""
struct PfullMixin{
    F <: Function,
    PRESSURE, # TODO: Add type annotation here!
    FIELDS, # TODO: Add type annotation here!
    PERM_MATRIX <: AbstractMatrix,
} <: AbstractHandlerMixin
    """Container storing the fields from the compute functions."""
    compute_fields::FIELDS

    # TODO: Can remove all of these fields and access from any one of the coordinates
    # style (problem is establishing a singleton then...)
    # TODO: I am not sure about these fields, since they are stored in the PfullCoordsStyle
    # It is helpful to have them in one place though to access instead of accessing them
    # in the writer (convience basically)
    """Function to compute the pressure field"""
    pfull_compute!::F # TODO: Add check that this is the same across all coords style

    """ClimaCore field of pressure created by pfull_compute!"""
    pfull_field::PRESSURE # TODO: Add check that this is the same across all coords style

    """A permutation matrix, created by sortperm, for
    sorting the pressures for each column"""
    perm_matrix::PERM_MATRIX # TODO: This is the same across all pressure_coords_style by virtue of pfull_field
end

# TODO: Simplify this, I think compute_fields and the writer is sufficient and
# this allows me to move the error handling to here instead of in the
# diagnostics handler (if additional functionality is wanted, then the mixin
# should be responsible for checking it)
# function PfullMixin(scheduled_diagnostics, compute_fields)

#     return PfullMixin{
#         typeof(compute_fields),
#         typeof(pfull_compute!),
#         typeof(pfull_field),
#         typeof(perm_matrix),
#     }(
#         compute_fields,
#         pfull_compute!,
#         pfull_field,
#         perm_matrix,
#     )
# end
