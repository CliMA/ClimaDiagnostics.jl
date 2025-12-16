# TODO: Will be used for the output writers
# Structually and logically, the output writers is being told to output the
# data in lon lat or x-y (no conversion style) or pressure style
# (pfullcoordsstyle)

"""
    abstract type AbstractCoordsStyle end

An object that tells the `NetCDFWriter` whether to convert the coordinates or
not to another coordinates system.

This determines how the diagnostic output from the `NetCDFWriter` is represented
spatially. For instance, no conversion is needed for outputting data on a
longitude-latitude-z grid, but conversion is needed for a
longitude-latitude-pressure levels grid.

There is only one `AbstractCoordsStyle` per `NetCDFWriter`.
"""
abstract type AbstractCoordsStyle end

"""
    NoConversionStyle

A coordinate style indicating that no coordinate conversion should be performed.

When using `NoConversionStyle`, diagnostic output is in the model's native
coordinates (e.g., longitude-latitude or x-y horizontal coordinates with the
model vertical levels).
"""
struct NoConversionStyle <: AbstractCoordsStyle end

"""
    PfullCoordsStyle{FT <: AbstractFloat}

A coordinate style for interpolating diagnostic output to pressure levels.

When using `PfullCoordsStyle`, diagnostic data on model vertical levels is
interpolated onto the specified pressure levels.
"""
struct PfullCoordsStyle{F, FT <: AbstractFloat, PFULL_FIELD <: ClimaCore.Fields.Field, DI, PRESSURE_COORDS <: AbstractMatrix,} <: AbstractCoordsStyle
    """A vector of pressure levels to indicate how pressure levels should be
       outputted"""
    pressure_levels::Vector{FT}

    """Compute function for pressure"""
    pfull_compute!::F # This should be the same across all PfullCoordsStyle

    """A ClimaCore.Field representing pressure. This is used at
    the end of the simulation to interpolate offline along the horizontal
    direction. However, this pay the price of allocating a pressure field""" # TODO: This could be made more efficient by singleton design pattern?
    pressure_field::PFULL_FIELD # This should be the same across all PfullCoordsStyle

    """A dictionary mapping diagnostics to arrays on CPU."""
    preallocated_output_arrays::DI

    """Two dimensional array of pressures"""
    pressure_coords::PRESSURE_COORDS
end

# TODO: Might be worth it to stuff things like units, attributes for the netcdf
# writer here too
# TODO: Not sure if I should include pfull_compute! here, since if I can then
# it is easy to compute!
"""
    PfullCoordsStyle(pressure_levels)

Construct a `PfullCoordsStyle` from an iterable of pressure levels.

The pressure levels must be in sorted in ascending or descending order.
This ignore `z_sampling_method`.
"""
function PfullCoordsStyle(Y, p, t, pfull_compute!; pfull_levels = era5_pressure_levels())
    pfull_levels = collect(pfull_levels)
    issorted(pfull_levels, rev = true) && reverse!(pfull_levels)
    issorted(pfull_levels) ||
    error("Pressure levels ($pfull_levels) are not sorted")
    pfull_field = pfull_compute!(nothing, Y, p, t)

    FT = eltype(pfull_field)
    preallocated_output_arrays = Dict{ScheduledDiagnostic, Matrix{FT}}()

    typeofarray = ClimaComms.array_type(pfull_field)
    pressure_coordinates = typeofarray{FT}(
        repeat(pfull_levels, 1, Spaces.ncolumns(axes(pfull_field))),
    )

    return PfullCoordsStyle(
        pfull_levels,
        pfull_compute!,
        pfull_field,
        preallocated_output_arrays,
        pressure_coordinates,
    )
end

#! format: off
"""
    era5_pressure_levels()

Return the pressure levels used by ERA5 whose units are Pa.
"""
era5_pressure_levels() = return 100.0 .* [
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
    925, 950, 975, 1000,
]
#! format: on
