module Interpolators

import ClimaCore: Fields, Operators, Remapping
import ClimaComms

"""
    PressureInterpolator

A vertical interpolator for interpolating from data on model vertical levels
using a pressure field defined at model full-levels onto the specified pressure
levels.

This is a wrapper over `ClimaCore.Remapping.PressureInterpolator`.
"""
struct PressureInterpolator{
    CENTER_FIELD,
    FACE_FIELD,
    REMAPPER <: Remapping.PressureInterpolator,
    TIME,
}
    """Scratch field for compute!. This is reused across all diagnostics whose
    compute! functions return a field on a center space."""
    center_scratch_field::CENTER_FIELD

    """Scratch field for compute!. This is reused across all diagnostics whose
    compute! functions return a field on a face space."""
    face_scratch_field::FACE_FIELD

    """Pressure interpolator object for interpolating fields to pressure
    coordinates from ClimaCore"""
    pressure_intp::REMAPPER

    """Last time the pressure interpolator was computed. Used for caching to
    avoid redundant updates."""
    last_t::TIME
end

"""
    PressureInterpolator(
        pfull_field,
        t;
        pressure_levels = era5_pressure_levels(),
        pressure_intp_kwargs = (;),
    )

Construct a `PressureInterpolator` from
- `pfull_field`, the pressure field defined on a center space,
- `t`: the current simulation time.

The keyword argument `pressure_levels` is a vector of pressure levels in
ascending or descending order. The default is the ERA5 pressure levels.
The keyword argument `pressure_intp_kwargs` is a `NamedTuple` of keyword
arguments which is passed to `ClimaCore.Remapper.PressureInterpolator`.
"""
function PressureInterpolator(
    pfull_field,
    t;
    pressure_levels = era5_pressure_levels(),
    pressure_intp_kwargs = (;),
)
    intp_c2f = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    face_pfull_field = @. intp_c2f(pfull_field)
    center_scratch_field = deepcopy(pfull_field)
    face_scratch_field = deepcopy(face_pfull_field)
    pfull_intp = Remapping.PressureInterpolator(
        pfull_field,
        pressure_levels;
        pressure_intp_kwargs...,
    )
    ref_t = Ref(t)
    return PressureInterpolator(
        center_scratch_field,
        face_scratch_field,
        pfull_intp,
        ref_t,
    )
end

"""
    force_update!(pfull_intp::PressureInterpolator, t)

Force update `pfull_intp` for interpolation.
"""
function force_update!(pfull_intp::PressureInterpolator, t)
    (; pressure_intp, last_t) = pfull_intp
    Remapping.update!(pressure_intp)
    last_t[] = t
    return nothing
end

"""
    update!(pfull_intp::PressureInterpolator, t)

Update `pfull_intp` for interpolation.

If called with the same `t`, then no computation is done.
"""
function update!(pfull_intp::PressureInterpolator, t)
    (; last_t) = pfull_intp
    last_t[] == t && return nothing
    force_update!(pfull_intp, t)
    return nothing
end

"""
    interpolate_to_pressure_coords!(
        dest::Fields.Field,
        field::Fields.Field,
        pfull_intp::PressureInterpolator,
    )

Interpolate `field` onto `dest` in-place according to `pfull_intp`.

It is the user's responsibility to call `update!` or `force_update!` before
calling this function.
"""
function interpolate_to_pressure_coords!(
    dest::Fields.Field,
    field::Fields.Field,
    pfull_intp::PressureInterpolator,
)
    (; pressure_intp) = pfull_intp
    Remapping.interpolate_pressure!(dest, field, pressure_intp)
    return nothing
end

"""
    interpolate_to_pressure_coords(
        field::Fields.Field,
        pfull_intp::PressureInterpolator,
    )

Interpolate `field` to pressure coordinates using `pfull_intp` and return
the interpolated field on a new space with the same horizontal grid as `field`,
but the vertical coordinate is pressure.

It is the user's responsibility to call `update!` or `force_update!` before
calling this function.
"""
function interpolate_to_pressure_coords(
    field::Fields.Field,
    pfull_intp::PressureInterpolator,
)
    (; pressure_intp) = pfull_intp
    dest = Remapping.interpolate_pressure(field, pressure_intp)
    return dest
end

#! format: off
"""
    era5_pressure_levels()

Return the pressure levels used by
[ERA5](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5) in Pa.

"""
era5_pressure_levels() = return 100.0 .* [
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
    925, 950, 975, 1000,
]
#! format: on

end
