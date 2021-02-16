using Test
using Dates
using Statistics
push!(LOAD_PATH, joinpath(pwd(), "../src/"))
using BlueCloudPlankton
using Glob

srcdir = dirname(pathof(BlueCloudPlankton))

include(joinpath(srcdir,"grid.jl"))

@info "datadir: $datadir"
if !isdir(datadir)
    mkpath(datadir)
end

@testset "Data reading" begin
    datafiletest = joinpath(datadir, "data_small.csv")
    if !isfile(datafiletest)
        download("https://dox.ulg.ac.be/index.php/s/s561Q6289sRD1xt/download", datafiletest)
    else
        @info("Test data file already downloaded")
    end
    lon, lat, dates, abundance, scientificNames = BlueCloudPlankton.read_data(datafiletest);

    @test length(lon) == 99
    @test lon[1] == -63.183
    @test lat[1] == 44.008
    @test dates[2] == DateTime(2017, 3, 13, 2, 50)
    @test scientificNames[4] == "Metridia lucens"
end

@testset "DIVAndNN" begin
    include(joinpath(srcdir,"DIVAndNN_analysis.jl"))
    include(joinpath(srcdir,"DIVAndNN_plot_res.jl"))
    # check presence of NetCDF files
    @test length(glob("*/*nc",resdir)) > 0
end
