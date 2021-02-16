
# grid
gridlon = -90.:0.5:40.
gridlat = 30.:0.5:80.

srcdir = dirname(pathof(BlueCloudPlankton))
basedir = joinpath(srcdir,"..");

# avoid network file system on BlueCloud
basedir = expanduser("~/BlueCloud-data")

datadir = joinpath(basedir,"data/")
resdir = joinpath(basedir,"results/")
figdir = joinpath(basedir,"figures/")

if !isdir(datadir)
    mkpath(datadir)
end

if !isdir(resdir)
    mkpath(resdir)
end

if !isdir(figdir)
    mkpath(figdir)
end


function maybedownload(url,fname)
    if !isfile(fname)
        @info "downloading $url"
        download(url,fname)
    else
        @info("$url is already downloaded")
    end
end
