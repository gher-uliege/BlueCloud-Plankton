# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# # BlueCloud Zooplankton Demonstrator
#
# This notebook plots the results from the previous notebook DIVAndNN_analysis.ipynb.

# %%
srcdir = @__DIR__
if !isfile(joinpath(srcdir,"grid.jl"))
    srcdir = get(ENV,"SRCDIR","/workspace/VREFolders/Zoo-Phytoplankton_EOV/Zooplankton_EOV/bluecloud-plankton-master/src/")
end
using VideoIO
using Images
using DIVAnd
using DIVAndNN
push!(LOAD_PATH,srcdir)
push!(LOAD_PATH,@__DIR__)
using BlueCloudPlankton
using Dates
using JSON
using PyPlot
using Glob
using Statistics
using NCDatasets

# %%
include(joinpath(srcdir,"grid.jl"))
datafile = joinpath(datadir, "data-cpr.csv")

# %% [markdown]
# Bathymetry for plotting

# %%
bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;
bx,by,b = DIVAnd.extract_bath(bathname,bathisglobal,gridlon,gridlat);

# %% [markdown]
# Load observations

# %%

lon, lat, dates, value, scientificNames = BlueCloudPlankton.read_data(datafile)
scientificname_accepted = unique(scientificNames)

# %% [markdown]
# Load analysis

# %%

function createtime(dirn)
    open(dirn) do f
        return mtime(f)
    end
end

expdir = last(sort(joinpath.(resdir,readdir(resdir)),by = createtime))

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
    value_analysis = nomissing(ds[sname * "_L1"][:],NaN)
    gridlon = ds["lon"][:]
    gridlat = ds["lat"][:]
    gridtime =
        if haskey(ds,"time")
            ds["time"][:]
        else
            # there no time dimension in product
            [nothing]
        end
    close(ds)

    value_binned = zeros(length(gridlon),length(gridlat),length(gridtime))

    for n = 1:length(gridtime)
        sel = scientificNames .== sname

        if gridtime[n] !== nothing
            sel = sel .& (gridtime[n] .<= dates .<= (gridtime[n] + Dates.Year(1)))
        end

        XY = DIVAnd.ndgrid(gridlon,gridlat)
        value_binned[:,:,n] = DIVAndNN.binobs((lon[sel],lat[sel]),value[sel],XY);
    end

    cl = quantile(value_binned[isfinite.(value_binned)],[0.01, 0.99])
    if cl[1] == 0
        cl = (cl[2]/100,cl[2])
        @warn "setting explicitly lower scale to $(cl[1]) for $sname"
    end
    norm = PyPlot.matplotlib.colors.LogNorm(vmin=cl[1], vmax=cl[2])

    for n = 1:length(gridtime)
        clf()
        suptitle =
            if gridtime[n] !== nothing
                "$sname $(Dates.year(gridtime[n]))"
            else
                sname
            end

        fig.suptitle(suptitle,style="italic")
        subplot(1,2,1)
        pcolormesh(gridlon,gridlat,value_binned[:,:,n]', norm = norm)
        title("Binned observations")
        decorate()


        subplot(1,2,2)
        pcolormesh(gridlon,gridlat,value_analysis[:,:,n]', norm = norm)
        title("Analysis")
        decorate()
        savefig(joinpath(figdir,"$sname-$n.png"))
    end
end

# %% [markdown]
# Plot the result for the first species

# %%
PyPlot.ioff()
filenames = glob("*interp.nc",expdir);
plotfield(filenames[1]);

# %% [markdown]
# Plot the all species

# %%
plotfield.(filenames);

# %% [markdown]
# Make an animation of the all species distribution

# %%
animation_format = "mp4"

for filename in filenames
    sname = split(basename(filename),"_")[2]
    ds = Dataset(filename)

    if haskey(ds,"time")
        gridtime = ds["time"][:]

        @info "Encode animation of $sname"

        imgnames = [joinpath(figdir,"$sname-$n.png") for n = 1:length(gridtime)]
        imgstack = [RGB.(Images.load(imgname)) for imgname in imgnames];
        VideoIO.save(joinpath(figdir,"$sname.$animation_format"),imgstack,framerate = 2);
    end

    close(ds)
end

@info "Figures have been saved in $(figdir). Consider to copy the files to a permanent storage."
