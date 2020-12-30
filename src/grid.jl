
# grid
gridlon = -90.:0.5:40.
gridlat = 30.:0.5:80.

datadir = "../data/"
resdir = "../results/"
figdir = "../figures"

if !isdir(datadir)
    mkdir(datadir)
end

if !isdir(resdir)
    mkdir(resdir)
end

if !isdir(figdir)
    mkdir(figdir)
end


function maybedownload(url,fname)
    if !isfile(fname)
        @info "downloading $url"
        download(url,fname)
    else
        @info("$url is already downloaded")
    end
end
