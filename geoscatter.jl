# Plot shoreline map of WNAT with scatter of max wind velocity and
# scatter of maximum water level errors at tide gauges
using CSV
using Shapefile, GeoInterface
using Plots
using NetCDF
using Statistics
# GR backend
gr()

filename = "../Downloads/GeoData/gshhg-shp-2.3.7/GSHHS_shp/i/GSHHS_i_L1.shp"
shp = open(filename, "r") do io
    read(io, Shapefile.Handle)
end

# plot the shorelines in the AOI
XR = -60; XL = -100; YD = 20; YU = 50;
plt = plot(shp.shapes[1:500], line=:solid, linewidth=0.1, color=:gray,
    xlims=(XL,XR), ylims=(YD,YU), dpi=600,
    title="Maximum WL difference for 2004")

# open maximum wind file and plot the maximum winds
filename = "ATT/maxwvel.63.nc"
y = ncread(filename, "y")
x = ncread(filename, "x")
W = ncread(filename, "wind_max")
# remove zero velocity regions
y= y[W .> 0]; x = x[W .> 0]; W = W[W .> 0]
# scale so that colors match the error scale in next scatter
        # saturate at 80% of maximum
W = W ./ (0.8*(maximum(W) - minimum(W)))
W = W .- minimum(W) .- 0.5
scatter!(x, y, marker_z=W, legend=false, m=2,
         markerstrokewidth=0.0, color=:speed)

# open the MaxWL file and plot the errors
filename = "ATT/DailyWLPlots_Mv13/MaxWLs.csv"
df = CSV.read(filename)
MLerrors = df.MaxWL_Mod - df.MaxWL_Obs
scatter!(df.Lon, df.Lat, marker_z=MLerrors, legend=false, m=3,
    markerstrokewidth=0.5, color=:balance, colorbar=true,
    clims=(-0.5,0.5), colorbar_title="[m]")
    
display(plt)
savefig("ATT/DailyWLPlots_Mv13/MaxWLerror")
