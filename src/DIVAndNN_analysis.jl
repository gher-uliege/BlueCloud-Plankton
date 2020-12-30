using DIVAnd
using DIVAndNN
using Random
using NCDatasets
using DelimitedFiles
using Statistics
push!(LOAD_PATH,@__DIR__)
using BlueCloudPlankton
using DataStructures
using Printf
using Dates
using JSON
using PyPlot

Random.seed!(1234)

include("grid.jl")

datafile = joinpath(datadir, "data.csv")

data, columnnames = readdlm(datafile, ',', header=true);
getcolumn(name) = data[:,findfirst(columnnames[:] .== name)]

# check the depth range of the data
minimumDepthInMeters = getcolumn("minimumDepthInMeters")
maximumDepthInMeters = getcolumn("maximumDepthInMeters")
@show extrema(minimumDepthInMeters)
@show extrema(maximumDepthInMeters)

#maybedownload("https://www.ncei.noaa.gov/thredds-ocean/fileServer/ncei/woa/temperature/decav/0.25/woa18_decav_t00_04.nc",
#              joinpath(datadir,"woa18_decav_t00_04.nc"))



# place a copy of the 24 files (12 month, temperature and salinity) from the following addresses
# in datadir
# https://files.seadatanet.org/climatologies/Global_Ocean/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_T_V2_1_1900_2017_025/
# https://files.seadatanet.org/climatologies/Global_Ocean/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_S_V2_1_1900_2017_025/

months = ["January","February","March","April","May","June","July","August","September","October","November","December"]


datalist = [
   (varname = "Temperature",
    urls = map(m -> replace("https://files.seadatanet.org/climatologies/Global_Ocean/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_T_V2_1_1900_2017_025/SDC_GLO_CLIM_T_V2_1_1900_2017_025_April.nc","April" => m),months)),

   (varname = "Salinity",
    urls = map(m -> replace("https://files.seadatanet.org/climatologies/Global_Ocean/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_TS_V2/SDC_GLO_CLIM_S_V2_1_1900_2017_025/SDC_GLO_CLIM_S_V2_1_1900_2017_025_April.nc","April" => m),months)),
]


for l = 1:length(datalist)
    varname = datalist[l].varname
    urls = datalist[l].urls
    interp_fname = joinpath(datadir,"$(lowercase(varname)).nc")
    if isfile(interp_fname)
        @info "$interp_fname already there"
        continue
    end

    fullnames = joinpath.(datadir,basename.(urls))

    # sadly, you need to manual download
    #maybedownload.(urls,fullnames)

    ds = NCDataset(fullnames,aggdim = "time")

    Tlon = ds["lon"][:]
    Tlat = ds["lat"][:]
    Tz = ds["depth"][:]

    # extract surface value
    k = 1
    @show Tz[k]

    Tmean = mean(nomissing(ds[varname][:,:,k,:],NaN),dims=3)[:,:,1];
    Tmean = DIVAnd.ufill(Tmean,isfinite.(Tmean))

    # check range
    @show extrema(Tmean)

    DIVAndNN.saveinterp((Tlon,Tlat),Tmean,(gridlon,gridlat),varname,interp_fname)
    close(ds)
end

# take all years
years = 0:3000
ndimensions = 2

# land-sea mask and domain parameters

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

# Interpolate the bathymetry
DIVAndNN.prep_bath(bathname,bathisglobal,gridlon,gridlat,datadir)


if ndimensions == 3
    mask = repeat(mask,inner=(1,1,length(years)))
end

if ndimensions == 3
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat,years);
else
    mask2,pmn,xyi = DIVAnd.domain(bathname,bathisglobal,gridlon,gridlat);
end


covars_coord = false
covars_const = true

# load covariables
covars_fname = [
    ("bathymetry.nc","batymetry",identity),
    #            ("dist2coast_subset.nc","distance",identity),
    #("Chlorophyll/chloro_reinterp.nc","chla",identity),
    #("oxygen_reinterp2.nc","oxygen",identity),
    #("salinity.nc","salinity",log),
    ("salinity.nc","salinity",identity),
    ("temperature.nc","temperature",identity),
    #("nitrogen.nc",      "nitrogen",identity),
    #("phosphate.nc",     "phosphate",identity),
    #("silicate.nc",      "silicate",identity),
]
#covars_fname = []
covars_fname = map(entry -> (joinpath(datadir,entry[1]),entry[2:end]...),covars_fname)

field = DIVAndNN.loadcovar((gridlon,gridlat),covars_fname;
                           covars_coord = covars_coord,
                           covars_const = covars_const)

if ndimensions == 3
    sz = size(field)
    field = repeat(reshape(field,(sz[1],sz[2],1,sz[3])),1,1,length(years),1)
end
DIVAndNN.normalize!(mask,field)

lon, lat, dates, abundance, scientificNames = BlueCloudPlankton.read_data(datafile)

scientificname_accepted = unique(scientificNames)

@show unique(scientificNames)

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
end

@show sum(for_cv)
@show sum(.!for_cv)

lent = 0. # years
if ndimensions == 3
    lent = 5.
end

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

len = 75e3
len = 200e3
len = 300e3

#for len = [50e3, 75e3, 100e3, 125e3]
#    for epsilon2ap = [1, 5, 10, 50, 100, 500]
#for len = [100e3]
#    for epsilon2ap = [0.1]

        outdir = joinpath(resdir,"results-ncovars$(length(covars_fname))-epsilon2ap$(epsilon2ap)-len$(len)-niter$(niter)-nlayers$(length(NLayers))-ndimensions$(ndimensions)")
        mkpath(outdir)

        nameindex = 1
        for nameindex in 1:length(scientificname_accepted)
        #for nameindex in 1:1

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

            @show value_a[1:min(end,10)]
            @show extrema(value_a)
            @show length(value_a)

            Random.seed!(1234)

            value_analysis = zeros(size(mask))

            if ndimensions == 3
                xobs_a = (lon_a,lat_a,time_a)
                xobs_cv = (lon_cv,lat_cv,time_cv)
                lenxy = (len,len,lent)
                analysis_grid = (gridlon,gridlat,years)
            else
                xobs_a = (lon_a,lat_a)
                xobs_cv = (lon_cv,lat_cv)
                lenxy = (len,len)
                analysis_grid = (gridlon,gridlat)
            end


            loss_iter = []
            val_iter = []

            function plotres(i,lossi,value_analysis,y,gradloss,out,iobssel,obspos)
                #@show extrema(value_analysis[isfinite.(value_analysis)])
                vp = DIVAndNN.validate_regression(analysis_grid,value_analysis,xobs_cv,value_cv)
                push!(loss_iter,lossi)
                push!(val_iter,vp)
	            @printf("| %10d | %30.5f | %30.5f |\n",i,lossi,vp)
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
	            plotres = plotres,
	            plotevery = plotevery,
                rmaverage = true,
                trainfrac = trainfrac,
                epsilon2_background = epsilon2_background,
            )

            vp = DIVAndNN.validate_regression(analysis_grid,value_analysis,xobs_cv,value_cv)
            @show vp

            outname = joinpath(outdir,"DIVAndNN_$(sname)_interp.nc")

            cpme = DIVAnd_cpme(mask,pmn,xyi,xobs_a,value_a,lenxy,epsilon2_cpme)
            DIVAnd.save(outname,(gridlon,gridlat),value_analysis,sname; relerr = cpme)

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
                        "loss_iter" => loss_iter,
                        "val_iter" => val_iter,
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
        end

#    end
#end


