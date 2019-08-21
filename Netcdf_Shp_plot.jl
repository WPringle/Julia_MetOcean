using NetCDF
using Plots
using Shapefile, GeoInterface
using QHull

gr()
#The high-level interface is quite similar to the Matlab NetCDF interface,
# reading files is done by:
filename = "../Downloads/Globus/UV10_PSFC_2004_Jun-Oct.nc"
z = ncinfo(filename)
NX = 229 #name,"/","south_north")
NY = 311 #ncgetatt(filename,"/","west_east")
y = dropdims(ncread(filename, "XLAT",[1, 1, 1],[NX, NY, 1]);dims=3)
x = dropdims(ncread(filename, "XLONG",[1, 1, 1],[NX, NY, 1]);dims=3)
P = dropdims(ncread(filename, "PSFC",[1, 1, 1],[NX, NY, 1]);dims=3)
ch = chull([x[:] y[:]])

# read the shoreline shapefile
filename = "../Downloads/GeoData/gshhg-shp-2.3.7/GSHHS_shp/i/GSHHS_i_L1.shp"
shp = open(filename, "r") do io
    read(io, Shapefile.Handle)
end

# plot all the mesh points
plt = plot(ch.points[ch.vertices,1],ch.points[ch.vertices,2],color=:black)

plot!(shp.shapes[1:100],line=:solid)

# plot the shoreline
#for ii = 1:length(shp.shapes)
#    plot!(shp.shapes[ii],linestyle=:solid,color=:black)
#    #plot!(GeoInterface.coordinates(shp.shapes[ii])[1],linestyle=:solid,color=:black)
#end

display(plt)
