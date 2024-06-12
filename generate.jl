"""
Script to download all relevant ICESat-2 and GEDI data and saves it to geoparquet (.pq) files.
"""

using SpaceLiDAR
using ProgressMeter
using GeoDataFrames
using DataFrames
using GeoParquet
using WellKnownGeometry
using CategoricalArrays

sw = (
    min_x=8.34818106762969,
    min_y=47.1363854217594,
    max_x=9.02445376048887,
    max_y=47.7137808935097,
)
nz = (
    min_x=174.5061981,
    min_y=-37.3015846,
    max_x=175.3114720,
    max_y=-36.8368691,
)
nl = (min_x=3.2, min_y=50.75, max_x=7.22, max_y=53.7)

icesat2_atl03_foldername = "/mnt/ec66e171-5639-4c62-9d2c-08e81c462669/icesat2/ATL03/v05/"
icesat_atl08_foldername1 = "/mnt/ec66e171-5639-4c62-9d2c-08e81c462669/icesat2/ATL08/v05/"
icesat_atl08_foldername2 = "/mnt/00276e90-b667-456f-8892-e8e66eb238c0/icesat2/ATL08/v05/"

gedi_foldername = "/mnt/d3770dc8-0b9d-4a21-84c7-d0a3dc5946c9/gedi/v2/"
areas = (; nz, nl, sw)


function topq(granule::SpaceLiDAR.ICESat2_Granule, fn, area)
    name, ext = splitext(fn)
    for track in SpaceLiDAR.icesat2_tracks
        fn = "$(name)_$track$ext"
        isfile(fn) && continue
        try
            fn8n = replace(granule.id, "ATL03" => "ATL08")
            fn8 = joinpath(icesat_atl08_foldername2, fn8n)
            isfile(fn8) || (fn8 = joinpath(icesat_atl08_foldername1, fn8n))
            isfile(fn8) || (@error "Can't find ATL08 file $fn8"; continue)
            granule8 = granule_from_file(fn8)
            nt = SpaceLiDAR.classify(granule, granule8, tracks=(track,))
            length(nt) == 0 && continue
            df = DataFrame(nt[1])
            length(df.longitude) == 0 && continue
            df.track = Vector(df.track)
            df.strong_beam = Vector(df.strong_beam)
            df.detector_id = Vector(df.detector_id)
            df.track = Vector(df.track)
            df.classification = String.(df.classification)
            df.angle = SpaceLiDAR.track_angle(df.longitude, df.latitude)
            # df.anglet = SpaceLiDAR.track_angle.(Ref(granule), df.latitude)
            SpaceLiDAR.in_bbox!(df, area)
            filter!(:classification => ==("ground"), df)
            filter!(:height => !isnan, df)
            length(df.longitude) == 0 && continue
            select!(df, Not(:classification))
            select!(df, Not(:detector_id))
            df.geom = getwkb.(SpaceLiDAR.Point.(df.longitude, df.latitude, df.height))
            GeoParquet.write(fn, df)
        catch ex
            if ex isa InterruptException |
                      error("Canceled.")
            else
                @error "$fn can't be parsed: $ex"
            end
        end
    end
end

function topq(granule::SpaceLiDAR.GEDI_Granule, fn, area)
    name, ext = splitext(fn)
    for track in SpaceLiDAR.gedi_tracks
        fn = "$(name)_$track$ext"
        isfile(fn) && continue
        try
            nt = SpaceLiDAR.points(granule, canopy=false, tracks=(track,))
            df = DataFrame(nt[1])
            length(df.longitude) == 0 && continue
            df.track = Vector(df.track)
            df.strong_beam = Vector(df.strong_beam)
            df.track = Vector(df.track)
            df.classification = Vector(df.classification)
            df.angle = SpaceLiDAR.track_angle(df.longitude, df.latitude)
            # df.anglet = SpaceLiDAR.track_angle.(Ref(granule), df.latitude)
            SpaceLiDAR.in_bbox!(df, area)
            filter!(:classification => ==("ground"), df)
            filter!(:height => !isnan, df)
            length(df.longitude) == 0 && continue
            select!(df, Not(:classification))
            df.geom = getwkb.(SpaceLiDAR.Point.(df.longitude, df.latitude, df.height))
            GeoParquet.write(fn, df)
        catch ex
            if ex isa InterruptException
                error("Canceled.")
            else
                @error "$fn can't be parsed: $ex"
            end
        end
    end
end

for (name, area) in pairs(areas)
    granules = find(:ICESat2, "ATL03", area)
    fgranules = SpaceLiDAR.instantiate(granules, icesat2_atl03_foldername)

    @showprogress "ICESat-2 granules for $name" for granule in reverse(fgranules)
        id, _ = splitext(granule.id)
        fn = joinpath(@__DIR__, "data/$name/$(id).pq")
        topq(granule, fn, area)
    end
end

for (name, area) in pairs(areas)
    granules = SpaceLiDAR.find(:GEDI, "GEDI02_A", area)
    fgranules = SpaceLiDAR.instantiate(granules, gedi_foldername)
    @showprogress "GEDI granules for $name" for granule in fgranules
        id, _ = splitext(granule.id)
        fn = joinpath(@__DIR__, "data/$(name)/$(id).pq")
        topq(granule, fn, area)
    end
end
