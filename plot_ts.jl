# Script that plots the time series of the model versus observations
# - observations are automatically extracted from tides and currents
# - only time series where both model and observations are sufficiently complete
#   are plotted
# - Plots daily maximum water levels from the timeseries in addition
#   to full timeseries and lowpass filtered timeseries
using NetCDF
using Plots
using Dates
using HTTP
using CSV
using Statistics
using TimeSeries

using FFTW
include("fourfilt.jl")

# using GR backend, behaves like matlab
gr()

nodatalength = 0.0 # decimal percentage of period (if non-zero then the lowpass filter does not work...)
# low pass filter length
          #days
lpwindow = 28*3600.0*24.0 # seconds
dpiset = 300

# filenames
filename = "ATT/Mv13.61.nc"
sn = ncread(filename, "station_name")
# read parts of the netcdf file, convert to datetimes
time = ncread(filename, "time")
zeta = ncread(filename, "zeta")
zeta[zeta .< -999] .= NaN
lon = ncread(filename, "x")
lat = ncread(filename, "y")
dt = DateTime(2004,6,15) + Dates.Second.(time)
# Make sure only from July 1st
zeta = zeta[:,dt .>= DateTime(2004,7,1)]
deleteat!(dt,dt .< DateTime(2004,7,1))
tmstep = (Dates.value(dt[2]-dt[1]))/1000 #seconds

# This is for reading observations through tidesandcurrents API
base = "https://tidesandcurrents.noaa.gov/api/datagetter?"
suf = "product=hourly_height&datum=MSL&units=metric&time_zone=gmt&application=ports_screen&format=csv"
sd  = join(["begin_date=",Dates.format(dt[1],"yyyymmdd"),"&"])
ed  = join(["end_date=",Dates.format(dt[end],"yyyymmdd"),"&"])

# legend names
linenames = [:Observed, :Model]

# initializing the scatter arrays
xp = Float64[]; yp = Float64[]; mv = Float64[]; ov = Float64[]

# loop over all stations
for ii = 1:size(sn,2)

    # only get stations in our AOI
    #if (lon[ii] > -75 || lon[ii] < -83 ||
    #lat[ii] > 36.5 || lat[ii] < 24); continue; end
    if (lon[ii] > -67 || lon[ii] < -100 ||
        lat[ii] > 46 || lat[ii] < 24); continue; end

    # if modelled doesn't have enough values
    if (sum(isnan.(zeta[ii,:])) > length(dt)*nodatalength); continue; end

    # get and plot the station data
    sta = join(["station=",rstrip(join(sn[:,ii])),"&"])
    res = HTTP.get(join([base,sd,ed,sta,suf]))
    mycsv = CSV.read(res.body)
    # if got error or missing time
    if (ismissing(mycsv[1,1]) || mycsv[1,1][1:5] == "Error"); continue; end
    dtn = DateTime.(mycsv[:,1],"yyyy-mm-dd HH:MM")
    # if observed doesn't have enough non-missing values
    if (sum(ismissing.(mycsv[:,2])) > length(dtn)*nodatalength); continue; end

    # Display the station IDs and indices
    display(ii)
    display(join(sn[:,ii]))

    # Process observed results
    zo = mycsv[:,2]
    zo[ismissing.(mycsv[:,2])] .= NaN
    # remove low-pass filtered timeseries
    tostep = (Dates.value(dtn[2]-dtn[1]))/1000 #seconds
    zofilt = fourfilt(zo,tostep,Inf,lpwindow)
    zoaf = zo .- zofilt
    # convert observed data into the timeseries array
    to = TimeArray(dtn,zoaf)
    # calculate daily maximum
    to = collapse(to,day,first,maximum)
    rename!(to, linenames[1])

    # Process modeled results
    zm = zeta[ii,:]
    # remove low-pass filtered timeseries
    zmfilt = fourfilt(zm,tmstep,Inf,lpwindow)
    zmaf = zm .- zmfilt
    # convert model data into the timeseries array
    tm = TimeArray(dt,zmaf)
    # calculate daily maximum
    tm = collapse(tm,day,first,maximum)
    rename!(tm, linenames[2])

    # plot the observations and model
    plt1 = plot(to,
        title = join(["Daily Maximum Water Levels at Station ID: ",
            rstrip(join(sn[:, ii]))]
        ),
        legend = :topright,
        ylabel = "Water Level (LMSL) [m]",
        xlabel = "DateTime [UTC]"
    )
    plot!(tm)

    plt2 = plot(dtn,zo)
    plot!(dt,zm)
    plot!(dtn,zofilt)
    plot!(dt,zmfilt)
    plt = plot(plt1,plt2,dpi=dpiset,layout = (2,1))

    display(plt) # display in plot pane

    # calculate max WL error
    vtm = values(tm)
    vto = values(to)
    vtm = maximum(vtm[.!isnan.(vtm)])
    vto = maximum(vto[.!isnan.(vto)])

    # add the location plot to show error at stations
    push!(xp,lon[ii])
    push!(yp,lat[ii])
    push!(mv,vtm)
    push!(ov,vto)

    savefig(join(["ATT/DailyWLPlots_Mv13/" rstrip(join(sn[:, ii]))]))

end

# write out the maximum water levels over this time period
using DataFrames
df = DataFrame(Lon = xp, Lat = yp, MaxWL_Obs = ov, MaxWL_Mod = mv)
CSV.write("ATT/DailyWLPlots_Mv13/MaxWLs.csv",df)

# I'm having trouble using scatter here because of clash between Plots and GR package
#plt = scatter(xp,yp,200*ones(length(xp)),error) #,legend=:none)
#display(plt) # display in plot pane
#savefig(join(["ATT/DailyWLPlots_Mv13/MaxWLError"))
