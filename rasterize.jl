"""
Script to generate grids at different resolutions and latitudes
to check the possible resolution of an ICESat-2 and GEDI DTM.

Note we use ICESat-2 ATL08 here instead of ATL03, as it doesn't
change the results (ATL08 is 100 m or 20 m aggregate of ATL03, the same or higher as our highest resolution check)
but saves us several TBs of processing.
"""

using SpaceLiDAR
using GeoArrays
using StaticArrays
using ProgressMeter
using DataFrames
using GeoDataFrames
using TiledIteration
using CSV
using StatsPlots
using GeoParquet

foldernamei = "/mnt/00276e90-b667-456f-8892-e8e66eb238c0/icesat2/ATL08/v05/"
foldernamei2 = "/mnt/ec66e171-5639-4c62-9d2c-08e81c462669/icesat2/ATL08/v05/"
foldernameg = "/mnt/57697bfb-d51e-43c4-a929-db614d893853/gedi/v2/"
foldernameg2 = "/mnt/d3770dc8-0b9d-4a21-84c7-d0a3dc5946c9/gedi/v2/"
offset = 103.0
bbox = (min_x=103.0, min_y=0.0, max_x=103.5, max_y=90.0)

# Useful for downloading the data beforehand
granules = SpaceLiDAR.find(:ICESat2, "ATL08", bbox)
files = [file for file in readdir(foldernamei) if splitext(file)[end] == ".h5"]
fgranules = filter(!g -> g.id in files, granules)
files = [file for file in readdir(foldernamei2) if splitext(file)[end] == ".h5"]
fgranules = filter(!g -> g.id in files, fgranules)
SpaceLiDAR.write_granule_urls!("atl08_103.txt", fgranules)

granules = SpaceLiDAR.find(:GEDI, "GEDI02_A", bbox)
files = [file for file in readdir(foldernameg) if splitext(file)[end] == ".h5"]
fgranules = filter(!g -> g.id in files, granules)
files = [file for file in readdir(foldernameg2) if splitext(file)[end] == ".h5"]
fgranules = filter(!g -> g.id in files, fgranules)

SpaceLiDAR.write_granule_urls!("gedi_103.txt", fgranules)


height = 100
width = 200
hstep = 10_000  # [m] so 100m resolution
wstep = 20_000  # [m] so 100m resolution
_, lat_step = SpaceLiDAR.shift(0, 0, 0, hstep)  # in degrees
lats = 0.0:1:89.0
for lat ∈ lats

    # Checkpointing
    gacfn = joinpath(@__DIR__, "data", "c_$lat.tif")
    gagfn = joinpath(@__DIR__, "data", "g_$lat.tif")
    gaifn = joinpath(@__DIR__, "data", "i_$lat.tif")
    # isfile(gacfn) && isfile(gagfn) && isfile(gaifn) && continue

    lon_step, _ = SpaceLiDAR.shift(offset, lat, 90, wstep)
    lbbox = (min_x=offset, min_y=lat, max_x=lon_step, max_y=lat + lat_step)

    # Used for plotting the map in plot_resolution.py
    # println("mpatches.Rectangle(xy=[$offset, $lat], width=$lon_step, height=$lat_step")

    # Find local data for bbox
    granulesi = find(:ICESat2, "ATL08", lbbox)
    granulesg = find(:GEDI, "GEDI02_A", lbbox)
    granulesg = isempty(granulesg) ? SpaceLiDAR.Granule[] : granulesg
    granulesi = isempty(granulesi) ? SpaceLiDAR.Granule[] : granulesi
    local_granulesi = SpaceLiDAR.instantiate(granulesi, foldernamei)
    local_granulesi = vcat(local_granulesi, SpaceLiDAR.instantiate(granulesi, foldernamei2))
    local_granulesg = SpaceLiDAR.instantiate(granulesg, foldernameg)
    local_granulesg = vcat(local_granulesg, SpaceLiDAR.instantiate(granulesg, foldernameg2))

    # Initialize empty rasters which we'll use for rasterization
    gac = GeoArray(zeros(Int32, width, height))
    GeoArrays.bbox!(gac, lbbox)
    epsg!(gac, 4326)
    gai = GeoArray(zeros(Int32, width, height))
    GeoArrays.bbox!(gai, lbbox)
    epsg!(gai, 4326)
    gag = GeoArray(zeros(Int32, width, height))
    GeoArrays.bbox!(gag, lbbox)
    epsg!(gag, 4326)

    @showprogress "Intersecting ICESat-2" 1 for granule in local_granulesi
        isfile(granule.url) || continue
        df = DataFrame(granule)
        select!(df, [:longitude, :latitude])
        SpaceLiDAR.in_bbox!(df, lbbox)
        for row in eachrow(df)
            try
                i, j = indices(gai, (row.longitude, row.latitude))
                gai[i, j, 1] += 1
                gac[i, j, 1] += 1
            catch y
                if y isa BoundsError
                    nothing
                else
                    throw(y)
                end
            end
        end
        GC.gc()
    end

    @showprogress "Intersecting GEDI" 1 for granule in local_granulesg
        isfile(granule.url) || continue

        # only for the author, who compressed some downloads into another format
        if filesize(granule.url) == 0
            granule.url = replace(granule.url, ".h5" => ".pq")
            df = GeoParquet.read(granule.url)
            select!(df, :x => :longitude, :y => :latitude)
        else
            df = DataFrame(granule)
            select!(df, :longitude, :latitude)
        end
        SpaceLiDAR.in_bbox!(df, lbbox)
        for row in eachrow(df)
            try
                i, j = indices(gag, (row.longitude, row.latitude))
                gag[i, j, 1] += 1
                gac[i, j, 1] += 1
            catch y
                if y isa BoundsError
                    nothing
                else
                    throw(y)
                end
            end
        end
        GC.gc()
    end

    GeoArrays.write(gacfn, gac)
    GeoArrays.write(gagfn, gag)
    GeoArrays.write(gaifn, gai)
end

# These loops are split because the first part can fail
# and we don't want to have to redo everything
dfs = Vector{NamedTuple}()
for lat ∈ lats
    gacfn = joinpath(@__DIR__, "data", "c_$lat.tif")
    gagfn = joinpath(@__DIR__, "data", "g_$lat.tif")
    gaifn = joinpath(@__DIR__, "data", "i_$lat.tif")
    gac = GeoArrays.read(gacfn)
    gag = GeoArrays.read(gagfn)
    gai = GeoArrays.read(gaifn)

    # Loop over 100m resolution raster with sized windows
    # each representing a lower resolution, up to 5km.
    # If any of the cells within a window has a count bigger than zero
    # it is taken as filled at that window resolution.
    for r ∈ (1, 2, 5, 7, 10, 15, 20, 30, 40, 50)
        ti = TileIterator(axes(gac), (r, r, 1))
        i = 0
        g = 0
        c = 0
        for tileaxs in ti
            i += any(gai.A[tileaxs[1:2]...] .> 0)
            g += any(gag.A[tileaxs[1:2]...] .> 0)
            c += any(gac.A[tileaxs[1:2]...] .> 0)
        end
        pc = c / length(ti) * 100
        pg = g / length(ti) * 100
        pi = i / length(ti) * 100
        push!(dfs, (; lat, r=r * 100, icesat2=pi, gedi=pg, combined=pc))
    end
end
df = DataFrame(dfs)
CSV.write("occurence.csv", df)

# Useful for plotting
sd = stack(df, [:icesat2, :gedi, :combined]; variable_name=:mission, value_name=:percentage)
CSV.write("occurencestacked.csv", sd)
