module BlueCloudPlankton

using DIVAnd
using PyPlot
using PyCall
using Dates
using Statistics
using DelimitedFiles
using NCDatasets

# Import stuff from cartopy
ccrs = pyimport("cartopy.crs")
gridliner = pyimport("cartopy.mpl.gridliner")
cfeature = pyimport("cartopy.feature")
mticker = pyimport("matplotlib.ticker")
myproj = ccrs.PlateCarree()

mpl = pyimport("matplotlib");
cartopyticker = pyimport("cartopy.mpl.ticker")
lon_formatter = cartopyticker.LongitudeFormatter()
lat_formatter = cartopyticker.LatitudeFormatter()
mpl.rc("axes", linewidth=2)
mpl.rc("font", weight="light", size=14)

"""
    read_data(datafile)

Add coastline and ticklabels on the figure as a global map

## Examples
```julia-repl
julia> lon, lat, dates, abundance, scientificNames = read_data(datafile)
```
"""
function read_data(datafile::String)
    getcolumn(name) = data[:,findfirst(columnnames[:] .== name)]

    data, columnnames = readdlm(datafile, ',', header=true);

    # Extract coordinates
    lat = getcolumn("decimalLatitude");
    lon = getcolumn("decimalLongitude");
    abundance = Float64.(getcolumn("abundance"));
    scientificNames = getcolumn("scientificName")

    # Parse dates
    df = DateFormat("y-m-d H:M:S.s");
    dates = DateTime.(getcolumn("eventDate"), df);

    return Float64.(lon), Float64.(lat), dates, Float64.(abundance), String.(scientificNames)

end

# Plotting functions
"""
    decorate_global_map(ax)

Add coastline and ticklabels on the figure as a global map

## Examples
```julia-repl
julia> decorate_global_map(ax)
```
"""
function decorate_global_map(ax, coast)
    PyPlot.grid(linewidth=0.2)
    ax.add_feature(coast, color=".6",
            edgecolor="k", zorder=5)
    ax.set_xticks(-180.:45.5:180.)
    ax.set_yticks(-90.:30:90.)
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

end

"""
    decorate_domain_map(ax)

Add coastline and ticklabels on the figure for the selected domain

## Examples
```julia-repl
julia> decorate_domain_map(ax)
```
"""
function decorate_domain_map(ax, coast)
    PyPlot.grid(linewidth=0.2)
    ax.add_feature(coast, color=".6",
            edgecolor="k", zorder=5)
    ax.set_xticks(-90.:20.:40.)
    ax.set_yticks(30.:10.:80.)
    ax.set_xlim(-90., 40.)
    ax.set_ylim(30., 80.)
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

end

"""
    plot_results(longrid, latgrid, field, speciesname, figname)

Create a pcolor plot of the interpolated field

## Examples
```julia-repl
julia> plot_results(longrid, latgrid, field, speciesname, figname)
```
"""
function plot_results(longrid, latgrid, field, speciesname="", figname="")
    fig = PyPlot.figure(figsize=(10, 10))
    ax = PyPlot.subplot(111, projection=myproj)
    pcm = PyPlot.pcolormesh(longrid, latgrid, field)
    cb = PyPlot.colorbar(pcm, orientation="vertical", shrink=0.35)
    decorate_domain_map(ax)
    title(speciesname)
    if length(figname) > 0
        PyPlot.savefig(figname, dpi=300, bbox_inches="tight")
    end
    PyPlot.close()
end

"""
    create_nc_results(filename, lon, lat, field; valex=-999.9)

Create a netCDF with the results inside.

## Examples
```julia-repl
julia> create_nc_results(filename, lons, lats, times, spm)
```
"""
function create_nc_results(filename::String, lons, lats, times, field,
    speciesname::String=""; valex=-999.9)
    Dataset(filename, "c") do ds

        # Dimensions
        ds.dim["lon"] = length(lons)
        ds.dim["lat"] = length(lats)
        ds.dim["time"] = Inf # unlimited dimension

        # Declare variables
        ncfield = defVar(ds,"LOGSPM", Float64, ("lon", "lat", "time"))
        ncfield.attrib["missing_value"] = Float64(valex)
        ncfield.attrib["_FillValue"] = Float64(valex)
        ncfield.attrib["long_name"] = "interpolated abundance"
        ncfield.attrib["units"] = "ind/m3"

        nctime = defVar(ds,"time", Float32, ("time",))
        nctime.attrib["missing_value"] = Float32(valex)
        nctime.attrib["units"] = "seconds since 1981-01-01 00:00:00"
        nctime.attrib["time"] = "time"

        nclon = defVar(ds,"lon", Float32, ("lon",))
        nclon.attrib["missing_value"] = Float32(valex)
        nclon.attrib["_FillValue"] = Float32(valex)
        nclon.attrib["units"] = "degrees East"
        nclon.attrib["lon"] = "longitude"

        nclat = defVar(ds,"lat", Float32, ("lat",))
        nclat.attrib["missing_value"] = Float32(valex)
        nclat.attrib["_FillValue"] = Float32(valex)
        nclat.attrib["units"] = "degrees North"
        nclat.attrib["lat"] = "latitude"

        # Global attributes
        ds.attrib["institution"] = "GHER - University of Liege"
        ds.attrib["title"] = "Interpolated field of $(speciesname)"
        ds.attrib["comment"] = "Original data prepared by VLIZ"
        ds.attrib["author"] = "C. Troupin (ctroupin@uliege), A. Barth (a.barth@uliege.be)"
        ds.attrib["tool"] = "create_nc_spm_tile.jl"
        ds.attrib["institution_url"] = "http://labos.ulg.ac.be/gher/"
        ds.attrib["institution_logo_url"] = "http://gher-diva.phys.ulg.ac.be/Images/gher-logo.png"

        # Define variables
        ncfield[:] = field
        nctime[:] = times
        nclon[:] = lons
        nclat[:] = lats;

    end
end;


export decorate_domain_map

end # module
