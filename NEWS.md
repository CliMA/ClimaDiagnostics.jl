# NEWS

v0.2.3
-------

- Detect and ignore duplicated diagnostics.


v0.2.2
-------

- Fix support for caches that are not NamedTuples in `NetCDFWriter`.


v0.2.1
-------

- Fix support for purely horziontal spaces in `NetCDFWriter`.

v0.2.0
-------

- ![][badge-💥breaking] The `NetCDFWriter` now outputs points on the vertical levels by default.
- ![][badge-💥breaking] `disable_vertical_interpolation` is removed in favor of `z_sampling_method`.

[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
