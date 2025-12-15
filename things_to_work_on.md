1. Add the functionality to close
2. Depending on where it is saved and how ClimaAnalysis work, ClimaAnalysis
   SimDir might need to be updated
3. Use conservative regridding instead of nearest neighbor
4. Uniform DiagnosticsHandler and the pressure version (or not I guess?)
5. Add error handling for spaces

Look into why the diagnostics doesn't work with both of them
    - I think it is because of the NetCDFWriter and how I can't share both of
      them




1. Get a level of a climacore field as scratch space (ClimaCore.Fields.level)
2. ClimaCore.DataLayouts.array2data(interp_arr[2,:], ClimaCore.Fields.field_values(surface_field))




Another thing that can be done to reduce the code 
