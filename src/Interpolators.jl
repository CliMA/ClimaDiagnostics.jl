module Interpolators

import ClimaCore: Fields, Operators, Remapping
import ClimaComms

"""
    PfullInterpolator

A vertical interpolator for interpolating from data on model vertical levels
using a pressure field defined at model full-levels onto the specified pressure
levels.
"""
struct PfullInterpolator{
    F,
    FIELD,
    FACE_FIELD,
    REMAPPER <: Remapping.PressureInterpolator,
    TIME,
}
    """In-place compute function for pressure field"""
    compute_pfull!::F

    """Scratch field for compute!. This is reused across all diagnostics whose
    compute! functions return a field on a center space."""
    center_scratch_field::FIELD

    """Scratch field for compute!. This is reused across all diagnostics whose
    compute! functions return a field on a face space."""
    face_scratch_field::FACE_FIELD

    """Pressure interpolator object for interpolating fields to pressure
    coordinates from ClimaCore"""
    pressure_intp::REMAPPER

    """Last time the pressure field is computed. This is used to avoid repeated
    computation of the pressure field."""
    last_t::TIME
end

"""
    PfullInterpolator(
        compute_pfull!,
        Y,
        p,
        t;
        pfull_levels = era5_pressure_levels(),
    )

Construct a `PfullInterpolator` from
- `compute_pfull!`, a function or struct that compute a pressure field,
- `Y`, the state vector,
- `p`: the cache and parameters,
- `t`: current simulation time.

The keyword argument `pfull_levels` is a vector of pressure levels in ascending
order. The default is the ERA5 pressure levels.

The signature of `compute_pfull!` should be of the form:
```julia
function compute_pfull!(out, state, cache, time)
    if isnothing(out)
        return state.pfull
    else
        out .= state.pfull
    end
end
```

Currently, the argument `compute_pfull!` does not support functions or structs
that return a `Base.Broadcast.Broadcasted`.
"""
function PfullInterpolator(
    compute_pfull!,
    Y,
    p,
    t;
    pfull_levels = era5_pressure_levels(),
)
    intp_c2f = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    center_pfull_field = compute_pfull!(nothing, Y, p, t)
    face_pfull_field = @. intp_c2f(center_pfull_field)
    center_scratch_field = deepcopy(center_pfull_field)
    face_scratch_field = deepcopy(face_pfull_field)
    pfull_intp =
        Remapping.PressureInterpolator(center_pfull_field, pfull_levels)
    ref_t = Ref(t)
    return PfullInterpolator(
        compute_pfull!,
        center_scratch_field,
        face_scratch_field,
        pfull_intp,
        ref_t,
    )
end

"""
    force_update!(pfull_intp::PfullInterpolator, Y, p, t)

Force update `pfull_intp` for interpolation.
"""
function force_update!(pfull_intp::PfullInterpolator, Y, p, t)
    (; compute_pfull!, pressure_intp, last_t) = pfull_intp
    pfull_field = Remapping.pfull_field(pressure_intp)
    compute_pfull!(pfull_field, Y, p, t)
    Remapping.update!(pressure_intp)
    last_t[] = t
    return nothing
end

"""
    update!(pfull_intp::PfullInterpolator, Y, p, t)

Update `pfull_intp` for interpolation.

If called with the same `t`, then no computation is done.
"""
function update!(pfull_intp::PfullInterpolator, Y, p, t)
    (; last_t) = pfull_intp
    last_t[] == t && return nothing
    force_update!(pfull_intp, Y, p, t)
    return nothing
end

"""
    interpolate_to_pfull_coords!!(
        dest::Fields.Field,
        field::Fields.Field,
        pfull_intp::PfullInterpolator,
    )

Interpolate `field` onto `dest` in-place according to `pfull_intp`.

This mutates both `dest` and `field`.

It is the user's responsibility to call `update!` or `force_update!` before
calling this function.
"""
function interpolate_to_pfull_coords!!(
    dest::Fields.Field,
    field::Fields.Field,
    pfull_intp::PfullInterpolator,
)
    (; pressure_intp) = pfull_intp
    Remapping.interpolate_pressure!!(dest, field, pressure_intp)
    return nothing
end

"""
    interpolate_to_pfull_coords!(
        field::Fields.Field,
        pfull_intp::PfullInterpolator,
    )

Interpolate `field` to pressure coordinates using `pfull_intp` and return
the interpolated field on a new space with the same horizontal grid as `field`,
but the vertical direction is pressure.

It is the user's responsibility to call `update!` or `force_update!` before
calling this function.
"""
function interpolate_to_pfull_coords!(
    field::Fields.Field,
    pfull_intp::PfullInterpolator,
)
    (; pressure_intp) = pfull_intp
    dest = Remapping.interpolate_pressure!(field, pressure_intp)
    return dest
end

#! format: off
"""
    era5_pressure_levels()

Return the pressure levels used by ERA5 in Pa.

See [here](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5).
"""
era5_pressure_levels() = return 100.0 .* [
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
    925, 950, 975, 1000,
]
#! format: on

end
