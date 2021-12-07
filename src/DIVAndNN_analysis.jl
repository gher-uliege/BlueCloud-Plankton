# -*- coding: utf-8 -*-
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
# <h1 style="text-align: center;">BlueCloud Zooplankton Demonstrator</h1>
#
# <p style="text-align: center;">Alexander Barth, Charles Troupin</p>
# <p style="text-align: center;">GHER, University of Liège, Belgium</p>
#
# The aim of this notebook is to created a gridded dataset of the
# [Continuous Plankton Recorder](https://www.cprsurvey.org/services/the-continuous-plankton-recorder/) within the Virtual Research Environement developped in BlueCloud.
#
#
# ### Method
#
# DIVA (Data Interpolating Variational Analysis) aim to derive a gridded climatology from in situ observations. The  derived field should be
#  __close to the observations__ (it should not necessarily pass through all observations because observations have errors), __close to a first guess estimate__, "__smooth__" (i.e. small first and second ordrer derivatives). This is formalized using the following cost functions.
#
# $\newcommand{\vec}{\mathbf}\newcommand{\mat}{\mathbf}$
# $$
# J(\vec x) =
# \left(\mat H \vec x - \vec y \right)^T \mat R^{-1}  \left(\mat H \vec x - \vec y \right)
# +
# \left(\vec x - \vec x^b \right)^T \mat B^{-1}  \left(\vec x - \vec x^b \right)
# $$
#
# DIVAnd typically use obervations from a single parameter to create gridded datasets.
# In this notebook we explore the use of multivariate analysis using neural networks and DIVAnd.
# Assuming that there is a list of covariables $\vec z_1, \vec z_2, ...$ related to the parameter of interest $\vec x$, we assume that gridded field can be writtes as:
#
# $$
# \vec x = \vec x' + f(\vec z_1,\vec z_2,..., \mat W_1, \vec b_1, \mat W_2, \vec b_2,...)
# $$
#
# where $f$ is a non-linear function of the known covariables and unknown parameters $\mat W_1, \vec b_1, \mat W_2, \vec b_2,...$)
# The structure of the function $f$ is given here by a neural network (multilayer perception).
# The field $\vec x'$ is also unknown. It is subjected to the smoothness contraints from DIVAnd.
#
# We want to minimize:
#
# $$\begin{eqnarray}
# J(\vec x',\mat W_1, \vec b_1, \mat W_2, \vec b_2,...) &=&
# \left(\mat H \vec x - \vec y \right)^T \mat R^{-1}  \left(\mat H \vec x - \vec y \right) \\\\
# && +
# \left(\vec x' - \vec x^b \right)^T \mat B^{-1}  \left(\vec x' - \vec x^b \right)
# \end{eqnarray}$$
#
#
# For every location $j$, the value of vector $\vec v^1$ are the co-variables at the location $j$.
# This vector is linearly transformed by a matrix $\mat W_k$ and a vector $\vec b_k$ and then a non-linear activation function is applied to each element element of the resulting vector (except for the last step).
#
# $$
# \mat v^{(k+1)}_j = g^{(k+1)}(\mat v^{(k)}_j \mat W_k + \vec b_k)
# $$
#
# Here, the weight do not dependent on space but the longitude and latitude are one of the covariables.
#
# The following "Co-variables" are used (with the data source in parenthesis):
# * Sea water temperature (SeaDataCloud)
# * Salinity (SeaDataCloud)
# * Distance from coast (NASA Goddard Space Flight Center)
# * Bathymetry (GEBCO)
# * Nitrate, Silicate and Phosphate (World Ocean Atlas 2018)
#
#
# The following figure gives a high-level overview of the different datasets involed in this notebook.
#
# ![overview](docs/overview.svg)
#
#
#
# _This notebook is release under the terms of the GPL version 2 (or later, at your option)_
#
# The first step is to install all dependencies (if necessary)

# %%
using Pkg

try
    # check if DIVAndNN is already installed
    using DIVAndNN
catch
    # install all dependencies
    pkg"add https://github.com/gher-ulg/DIVAndNN.jl"
    pkg"add JSON PyCall PyPlot DIVAnd Glob DataStructures NCDatasets"
end

# %% [markdown]
# Load the modules

# %%
srcdir = @__DIR__
if !isfile(joinpath(srcdir,"grid.jl"))
    srcdir = get(ENV,"SRCDIR","/workspace/VREFolders/Zoo-Phytoplankton_EOV/DIVAndNN/bluecloud-plankton-master/src/")
end
using DIVAnd
using DIVAndNN
using Random
using NCDatasets
using DelimitedFiles
using Statistics
push!(LOAD_PATH,srcdir)
push!(LOAD_PATH,@__DIR__)
using BlueCloudPlankton
using DataStructures
using Printf
using Dates
using JSON
using PyPlot

# %% [markdown]
# Set seed for random number generator and include grid information

# %%
Random.seed!(1234)
include(joinpath(srcdir,"grid.jl"))

# %% [markdown]
# Data frome the continuous plankton recorder
# https://www.cprsurvey.org/services/the-continuous-plankton-recorder/

# %%
datafile = joinpath(datadir, "data-cpr.csv")
maybedownload("https://dox.ulg.ac.be/index.php/s/ME3U5wPdPK8GRFu/download",
              datafile)


# %% [markdown]
# Check the depth range of the observations

# %%
data, columnnames = readdlm(datafile, ',', header=true);
getcolumn(name) = data[:,findfirst(columnnames[:] .== name)]

# %%
# check the depth range of the data
minimumDepthInMeters = getcolumn("minimumDepthInMeters")
maximumDepthInMeters = getcolumn("maximumDepthInMeters")
@info "range of minimum depth $(extrema(minimumDepthInMeters))"
@info "range of maximum depth $(extrema(maximumDepthInMeters))"

# %% [markdown]
# The prepare following fields:
# * [Temperature and salinity from SeaDataCloud](https://doi.org/10.13155/77512)
# * Nitrate, silicate and phosphate from [World Ocean Atlas 2018](https://www.ncei.noaa.gov/products/world-ocean-atlas)

# %%
data_TS = [
    ("http://www.ifremer.fr/erddap/griddap/SDC_GLO_CLIM_TS_V2_1","Salinity","salinity"),
    ("http://www.ifremer.fr/erddap/griddap/SDC_GLO_CLIM_TS_V2_1","Temperature","temperature"),
    ("https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/nitrate/all/1.00/woa18_all_n00_01.nc","n_an","nitrate"),
    ("https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/silicate/all/1.00/woa18_all_i00_01.nc","i_an","silicate"),
    ("https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/phosphate/all/1.00/woa18_all_p00_01.nc","p_an","phosphate"),
]

maybedownload("https://dox.ulg.ac.be/index.php/s/7zwCEszAPIFeBAm/download",joinpath(datadir,"salinity.nc"))
maybedownload("https://dox.ulg.ac.be/index.php/s/OQMYYGFCEtS3xc9/download",joinpath(datadir,"temperature.nc"))

DIVAndNN.prep_tempsalt(gridlon,gridlat,data_TS,datadir)

# %% [markdown]
# Take all years

# %%
years = 0:3000
ndimensions = 2

#years = 1990:2017
#ndimensions = 3


# %% [markdown]
# Get bathymetry from GEBCO and prepare land-sea mask

# %%
bathname = joinpath(datadir,"gebco_30sec_4.nc");
bathisglobal = true;
maybedownload("https://dox.ulg.ac.be/index.php/s/RSwm4HPHImdZoQP/download",
              joinpath(datadir,"gebco_30sec_4.nc"))

maskname = joinpath(datadir,"mask.nc");

if !isfile(maskname)
    DIVAndNN.prep_mask(bathname,bathisglobal,gridlon,gridlat,years,maskname)
end

ds = Dataset(maskname,"r")
mask = nomissing(ds["mask"][:,:]) .== 1
close(ds)

DIVAndNN.prep_bath(bathname,bathisglobal,gridlon,gridlat,datadir)

if ndimensions == 3
    mask = repeat(mask,inner=(1,1,length(years)))
end

if ndimensions == 3
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat,years);
else
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat);
end;


# %% [markdown]
# Distance to coast

# %%
interp_fname = joinpath(datadir,"dist2coast_subset.nc")
fname_dist2coast = "https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/dist2coast_1deg"
if !isfile(interp_fname)
    DIVAndNN.prep_dist2coast(fname_dist2coast,gridlon,gridlat,interp_fname)
end

# %% [markdown]
# load covariables for the neural network

# %%
covars_coord = false
covars_const = true

covars_fname = [
    ("bathymetry.nc",       "batymetry",  identity),
    ("dist2coast_subset.nc","distance",   identity),
    #("salinity.nc","salinity",log),
    ("salinity.nc",         "salinity",   identity),
    ("temperature.nc",      "temperature",identity),
    ("nitrate.nc",          "nitrate",    identity),
    ("phosphate.nc",        "phosphate",  identity),
    ("silicate.nc",         "silicate",   identity),
]

covars_fname = map(entry -> (joinpath(datadir,entry[1]),entry[2:end]...),covars_fname)

field = DIVAndNN.loadcovar((gridlon,gridlat),covars_fname;
                           covars_coord = covars_coord,
                           covars_const = covars_const)

if ndimensions == 3
    sz = size(field)
    field = repeat(reshape(field,(sz[1],sz[2],1,sz[3])),1,1,length(years),1)
end
DIVAndNN.normalize!(mask,field)


# %% [markdown]
# Read data files

# %%
lon, lat, dates, abundance, scientificNames = BlueCloudPlankton.read_data(datafile);

# %% [markdown]
# all unique scientific names

# %%
scientificname_accepted = unique(scientificNames)
println.(unique(scientificNames));


# %% [markdown]
# Randomly choose cross-validation point

# %%
cvfname = replace(datafile,".csv" => "-cv.nc")
validation_fraction = 0.2

if isfile(cvfname)
    @info "load $cvfname"

    for_cv =
        Dataset(cvfname,"r") do ds
            Bool.(ds["for_cv"][:])
        end
else
    @info "split data"
    for_cv = rand(length(lon)) .< validation_fraction

    Dataset(cvfname,"c") do ds
        defVar(ds,"for_cv",Int8.(for_cv),("observation",),attrib = OrderedDict(
            "long_name" => "0: observation used for analysis; 1 observation used for validation")
               )
    end
end;

# %% [markdown]
# Number of data used for cross-validation and for the analysis

# %%
@info "Number of cross-valiation points: $(sum(for_cv))"
@info "Number of data points for the analysis: $(sum(.!for_cv))"

# %% [markdown]
# time correlation (for 3d analyses)

# %%
lent = 0. # years
if ndimensions == 3
    lent = 3.
end

# %%
plotevery = 100
niter = 500
trainfrac = 1.

epsilon2ap = 10
epsilon2_background = 10
epsilon2_cpme = epsilon2_background

NLayers = [size(field)[end],4,1]

learning_rate = 0.001
L2reg = 0.0001
dropoutprob = 0.6

len = 300e3;


# %%
outdir = joinpath(resdir,"results-ncovars$(length(covars_fname))-epsilon2ap$(epsilon2ap)-len$(len)-niter$(niter)-nlayers$(length(NLayers))-ndimensions$(ndimensions)")
mkpath(outdir)

nameindex = 1
for nameindex in 1:length(scientificname_accepted)

    sname = String(scientificname_accepted[nameindex])
    global loss_iter
    global val_iter
    @info sname

    paramname = joinpath(outdir,"DIVAndNN_$(sname)_interp.json")

    if isfile(paramname)
        #                continue
    end

    s = (sname .== scientificNames) .& .!for_cv
    lon_a,lat_a,obstime_a,value_a = lon[s], lat[s], dates[s], abundance[s]
    @info "number of observation used in analysis: $(sum(s))"

    s = (sname .== scientificNames) .& for_cv
    lon_cv,lat_cv,obstime_cv,value_cv = lon[s], lat[s], dates[s], abundance[s]
    @info "number of observation used for validation: $(sum(s))"

    time_a = Float64.(Dates.year.(obstime_a))
    time_cv = Float64.(Dates.year.(obstime_cv))

    @debug begin
        @show value_a[1:min(end,10)]
        @show extrema(value_a)
        @show length(value_a)
    end

    Random.seed!(1234)

    value_analysis = zeros(size(mask))

    if ndimensions == 3
        xobs_a = (lon_a,lat_a,time_a)
        xobs_cv = (lon_cv,lat_cv,time_cv)
        lenxy = (len,len,lent)
        analysis_grid = (gridlon,gridlat,years)
        analysis_grid2 = (gridlon,gridlat,DateTime.(years,1,1))
    else
        xobs_a = (lon_a,lat_a)
        xobs_cv = (lon_cv,lat_cv)
        lenxy = (len,len)
        analysis_grid = (gridlon,gridlat)
        analysis_grid2 = (gridlon,gridlat)
    end


    value_analysis,fw0 = DIVAndNN.analysisprob(
        mask,pmn,xyi,xobs_a,
        value_a,
        lenxy,epsilon2ap,
        field,
        NLayers,
        costfun = DIVAndNN.regression,
        niter = niter,
        dropoutprob = dropoutprob,
        L2reg = L2reg,
        learning_rate = learning_rate,
        rmaverage = true,
        trainfrac = trainfrac,
        epsilon2_background = epsilon2_background,
    )

    vp = DIVAndNN.validate_regression(analysis_grid,value_analysis,xobs_cv,value_cv)
    outname = joinpath(outdir,"DIVAndNN_$(sname)_interp.nc")

    cpme = DIVAnd_cpme(mask,pmn,xyi,xobs_a,value_a,lenxy,epsilon2_cpme)
    DIVAnd.save(outname,analysis_grid2,value_analysis,sname; relerr = cpme)

    open(paramname,"w") do f
        write(f,JSON.json(
            Dict(
                "validation" => vp,
                "L2reg" =>            L2reg,
                "dropoutprob" =>      dropoutprob,
                "epsilon2ap" =>       epsilon2ap,
                "epsilon2_background" => epsilon2_background,
                "len" =>              len,
                "niter" =>            niter,
                "learning_rate" =>    learning_rate,
                "NLayers" =>    NLayers,
                "name" =>    sname,
                "covars" => first.(covars_fname),
            )
        ))
    end

end

score = DIVAndNN.summary(outdir)

paramname2 = joinpath(outdir,"DIVAndNN.json")

open(paramname2,"w") do f
    write(f,JSON.json(
        Dict(
            "validation" => score,
            "L2reg" =>            L2reg,
            "dropoutprob" =>      dropoutprob,
            "epsilon2ap" =>       epsilon2ap,
            "epsilon2_background" => epsilon2_background,
            "len" =>              len,
            "niter" =>            niter,
            "learning_rate" =>    learning_rate,
            "NLayers" =>    NLayers,
            "covars" => first.(covars_fname),
        )
    ))
end;


@info "Results have been saved in $(outdir). Consider to copy the files to a permanent storage."
