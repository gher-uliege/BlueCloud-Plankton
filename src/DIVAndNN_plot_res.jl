# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.4
#   kernelspec:
#     display_name: Julia 1.5.1
#     language: julia
#     name: julia-1.5
# ---

# %% [markdown]
# # BlueCloud Zooplankton Demonstrator
#
# This notebook plots the results from the previous notebook DIVAndNN_analysis.ipynb.

# %%
using DIVAnd
using DIVAndNN
push!(LOAD_PATH,@__DIR__)
using BlueCloudPlankton
using Dates
using JSON
using PyPlot
using Glob
using Statistics
using NCDatasets

# %%

include("grid.jl")
datafile = joinpath(datadir, "data-cpr.csv")

# %% [markdown]
# Bathymetry for plotting

# %%
bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,gridlon,gridlat);

# %% [markdown]
# Load observations and analysis

# %%
expdir = joinpath(resdir,"results-ncovars3-epsilon2ap10-len300000.0-niter500-nlayers3-ndimensions2")
lon, lat, dates, value, scientificNames = BlueCloudPlankton.read_data(datafile)
scientificname_accepted = unique(scientificNames)
@show unique(scientificNames)

# %% [markdown]
# Helper functions for plotting

# %%
function decorate()
    colorbar(orientation="horizontal")
    contourf(bx,by,b' .> 0, levels=[0,.5], cmap = "gray")
    gca().set_aspect(1/cosd(mean(gridlat)))
end

# %%
function plotfield(filename)
    fig = figure(figsize = (10,5))
    sname = split(basename(filename),"_")[2]

    ds = Dataset(filename)
    value_analysis = nomissing(ds[sname * "_L1"][:,:],NaN)
    gridlon = ds["lon"][:]
    gridlat = ds["lat"][:]
    close(ds)

    sel = scientificNames .== sname
    XY = DIVAnd.ndgrid(gridlon,gridlat)
    value_binned = DIVAndNN.binobs((lon[sel],lat[sel]),value[sel],XY);

    cl = quantile(value_binned[isfinite.(value_binned)],[0.01, 0.99])
    if cl[1] == 0
        cl = (cl[2]/100,cl[2])
        @warn "setting explicitly lower scale to $(cl[1]) for $sname"
    end

    norm = PyPlot.matplotlib.colors.LogNorm(vmin=cl[1], vmax=cl[2])

    clf()
    fig.suptitle(sname,style="italic")
    subplot(1,2,1)
    pcolormesh(gridlon,gridlat,value_binned', norm = norm)
    title("Binned observations")
    decorate()


    subplot(1,2,2)
    pcolormesh(gridlon,gridlat,value_analysis', norm = norm)
    title("Analysis")
    decorate()
    savefig(joinpath(figdir,"$sname.png"))
end

# %% [markdown]
# Plot the result for the first species

# %%
filenames = glob("*nc",expdir);
plotfield(filenames[1])

# %% [markdown]
# Plot the all species

# %%
plotfield.(filenames)
