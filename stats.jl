"""
Script to generate statistics and plots.
"""

using SpaceLiDAR
using DataFrames
using GeoParquet
using Statistics
using ProgressMeter
# using CategoricalArrays
using CSV
using Parquet2


areas = ("sw", "nz", "nl")
datadir = "/mnt/ec66e171-5639-4c62-9d2c-08e81c462669/validation/"

outliers = readlines(joinpath(@__DIR__, "outliers.txt"))

const landcoverd = Dict{Int16,String}(
    10 => "Tree cover",
    20 => "Shrubland",
    30 => "Grassland",
    40 => "Cropland",
    50 => "Built-up",
    60 => "Bare / sparse vegetation",
    70 => "Snow and ice",
    80 => "Permanent water bodies",
    90 => "Herbaceous wetland",
    95 => "Mangroves",
    100 => "Moss and lichen",
)

dfs = Vector{NamedTuple}()
for area in areas
    granules = readdir(joinpath(datadir, area), join=true)
    @showprogress "Checking $area" 1 for granule in granules
        name = first(splitext(basename(granule)))
        mission = startswith(name, "ATL03") ? "ICESat-2" : "GEDI"
        id, beam = rsplit(name, '_', limit=2)
        id in outliers && continue
        df = GeoParquet.read(granule)
        diff = df.ref_center
        strong = df.strong_beam
        day = df.sun_angle .> 0
        landcovers = get.(Ref(landcoverd), df.landcover, Ref(""))
        mission = fill(mission, length(landcovers))
        areasf = fill(area, length(landcovers))
        # landcovers = CategoricalArray{String,1,Int16}(get.(Ref(landcoverd), df.landcover, Ref("")))
        # mission = CategoricalArray{String,1,Int8}(fill(mission, length(landcovers)))
        # areasf = CategoricalArray{String,1,Int8}(fill(area, length(landcovers)))
        push!(dfs, (; diff, landcovers, mission, areasf, strong, day))
    end
end
df = reduce(vcat, DataFrame.(dfs))
subset!(df, :diff => x -> .!isnan.(x))  # handful of NaNs...

Parquet2.writefile(joinpath(@__DIR__, "stats.pq"), df, compression_codec=:zstd)

bias(x) = mean(x)
mae(x) = mean(abs.(x))
rmse(x) = sqrt(mean(x .^ 2))

i, g = ("ICESat-2", "GEDI")
# Table difference
sdf = Vector{NamedTuple}()
for landcover in values(landcoverd)
    diffi = @view df.diff[(df.landcovers.==landcover).&(df.mission.==i)]
    diffg = @view df.diff[(df.landcovers.==landcover).&(df.mission.==g)]
    push!(
        sdf,
        (;
            landcover,
            biasi=-bias(diffi),
            biasg=-bias(diffg),
            maei=mae(diffi),
            maeg=mae(diffg),
            rmsei=rmse(diffi),
            rmseg=rmse(diffg),
            lengthi=length(diffi),
            lengthg=length(diffg)
        ),
    )
end
diffi = @view df.diff[df.mission.==i]
diffg = @view df.diff[df.mission.==g]
push!(
    sdf,
    (;
        landcover="Overall",
        biasi=-bias(diffi),
        biasg=-bias(diffg),
        maei=mae(diffi),
        maeg=mae(diffg),
        rmsei=rmse(diffi),
        rmseg=rmse(diffg),
        lengthi=length(diffi),
        lengthg=length(diffg)
    ),
)

ldf = DataFrame(sdf)
CSV.write(joinpath(@__DIR__, "overall_stats.csv"), ldf)

# Overall combined
-bias(df.diff)  # -0.02
mae(df.diff)  # 0.49
rmse(df.diff)  # 1.84
length(df.diff)  # 128013864

# Table for Appendix with difference for each area
for area in areas
    sdf = Vector{NamedTuple}()
    for landcover in values(landcoverd)
        diffi = @view df.diff[(df.landcovers.==landcover).&(df.mission.==i).&(df.areasf.==area)]
        diffg = @view df.diff[(df.landcovers.==landcover).&(df.mission.==g).&(df.areasf.==area)]
        push!(
            sdf,
            (;
                landcover,
                biasi=-bias(diffi),
                biasg=-bias(diffg),
                maei=mae(diffi),
                maeg=mae(diffg),
                rmsei=rmse(diffi),
                rmseg=rmse(diffg),
                lengthi=length(diffi),
                lengthg=length(diffg)
            ),
        )
    end
    diffi = @view df.diff[(df.mission.==i).&(df.areasf.==area)]
    diffg = @view df.diff[(df.mission.==g).&(df.areasf.==area)]
    push!(
        sdf,
        (;
            landcover="Overall",
            biasi=-bias(diffi),
            biasg=-bias(diffg),
            maei=mae(diffi),
            maeg=mae(diffg),
            rmsei=rmse(diffi),
            rmseg=rmse(diffg),
            lengthi=length(diffi),
            lengthg=length(diffg)
        ),
    )

    ldf = DataFrame(sdf)
    CSV.write(joinpath(@__DIR__, "overall_stats_$(area).csv"), ldf)
end


area = "sw"
dfs = Vector{NamedTuple}()
granules = readdir(joinpath(datadir, area), join=true)

@showprogress "Checking slopes" 1 for granule in granules
    name = first(splitext(basename(granule)))
    mission = startswith(name, "ATL03") ? "ICESat-2" : "GEDI"
    id, beam = rsplit(name, '_', limit=2)
    id in outliers && continue
    df = GeoParquet.read(granule)
    diff = df.ref_center
    strong = df.strong_beam
    day = df.sun_angle .> 0

    landcovers = get.(Ref(landcoverd), df.landcover, Ref(""))
    mission = fill(mission, length(landcovers))
    areasf = fill(area, length(landcovers))
    # landcovers = CategoricalArray{String,1,Int16}(get.(Ref(landcoverd), df.landcover, Ref("")))
    # mission = CategoricalArray{String,1,Int8}(fill(mission, length(landcovers)))
    # areasf = CategoricalArray{String,1,Int8}(fill(area, length(landcovers)))
    push!(dfs,
        (;
            diff,
            landcovers,
            mission,
            areasf,
            strong,
            day,
            df.slope,
            df.ref_f25,
            df.ref_f5,
            df.ref_f75,
            df.ref_f10,
            df.ref_b25,
            df.ref_b5,
            df.ref_b75,
            df.ref_b10,
            df.ref_l25,
            df.ref_l5,
            df.ref_l75,
            df.ref_l10,
            df.ref_r25,
            df.ref_r5,
            df.ref_r75,
            df.ref_r10
        ))
end
df = reduce(vcat, DataFrame.(dfs))
subset!(df, :diff => x -> .!isnan.(x))  # handful of NaNs...
Parquet2.writefile(joinpath(@__DIR__, "stats_slope.pq"), df, compression_codec=:zstd)
