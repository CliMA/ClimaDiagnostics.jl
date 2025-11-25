using ClimaInterpolations
import ClimaInterpolations.Interpolation1D: interpolate1d!, Linear, Flat

xsource = [1.0, 2.0, 3.0, 4.0]
fsource = [10.0, 3.0, 1.0, 4.0]

xtarget = [2.0, 3.0]
ftarget = [1000.0, 123123.0]
interpolate1d!(ftarget, xsource, xtarget, fsource, Linear(), Flat())

xsource = [[1.0, 2.0] [3.0, 4.0]]
fsource = [[10.0, 3.0] [1.0, 4.0]]

xtarget = [[2.0, 3.0] [3.0, 4.0]]
ftarget = [[0.0, 0.0] [0.0, 0.0]]
interpolate1d!(ftarget, xsource, xtarget, fsource, Linear(), Flat())

xsource = permutedims(
    [[1.0, 10.0, 100.0] [2.0, 11.0, 101.0] [3.0, 12.0, 102.0]],
    (2, 1),
)
fsource = permutedims(
    [[100.0, 1.0, 11.0] [50.0, 2.0, 12.0] [10.0, 3.0, 13.0]],
    (2, 1),
)

# xtarget = reshape([1.5, 11.5, 100.5], 1, 3)
xtarget = reshape([2.5, 12.5, 101.0], 1, 3)

ftarget = reshape([-1.0, -1.0, -1.0], 1, 3)
ClimaInterpolations.Interpolation1D.interpolate1d!(
    ftarget,
    xsource,
    xtarget,
    fsource,
    ClimaInterpolations.Interpolation1D.Linear(),
    ClimaInterpolations.Interpolation1D.Flat(),
)
