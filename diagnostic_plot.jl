# Illustrations for SolarIrradiance documentation
# Copyright (c) 2021 Pierre Haessig

using SolarIrradiance
using PyPlot


# Load GHI data 

# 2012, day 136: sun & clouds
year = 2012
n = 136
GHI_day = [0.0, 0.0, 0.0, 0.0, 0.0, 39.18, 209.0, 374.0, 240.0, 476.0, 379.0, 619.0, 872.01, 664.0, 680.0, 434.0, 311.0, 313.0, 169.0, 20.88, 0.0, 0.0, 0.0, 0.0]

# 2012, day 137: sun
n = 137
GHI_day = [0.0, 0.0, 0.0, 0.0, 0.0, 64.0, 218.0, 396.0, 566.0, 709.0, 827.01, 892.01, 851.01, 834.01, 756.01, 666.0, 519.0, 346.0, 168.0, 21.88, 0.0, 0.0, 0.0, 0.0]

## Plot


## Transpose GHI for a given slope & azimuth and plot ##

# Time of the measurement
tc = 12 # 0:23
dt = 1/60 # 1 minute ~ instantaneous (0.0 is orbidden)

# Location (Rennes, France):
lat = 48.117 #°
lon = -1.678 #°

# Panel orientation:
slope = 40 #°
azimuth = 0:10:180 #90 #° → west-facing
albedo = 0.5; # in [0,1]

Gplane = zeros(length(tc))

# broadcast to common shape
n = n .+ zero(tc) .+ zero(azimuth)
tc = tc .+ zero(azimuth)
GHI = GHI_day[13] .+ zero(azimuth)

ω = hour_angle.(n, tc, lon)

Gplane = global_radiation_tilt.(GHI, n, tc, dt, lat, lon, slope, azimuth, albedo)

ωud = day_bounds.(n, lat)
ωud_plane = day_bounds.(n, slope, azimuth, lat)
ωud_plane2 = day_bounds2.(n, slope, azimuth, lat)
# split the tuples
ωu = first.(ωud); ωd = last.(ωud);
ωu_plane = first.(ωud_plane); ωd_plane = last.(ωud_plane);
ωu_plane2 = first.(ωud_plane2); ωd_plane2 = last.(ωud_plane2);

cosθ = cosθavg.(n, tc, dt, slope, azimuth, lat, lon)
cosθZ = cosθZavg.(n, tc, dt, lat, lon)

# Plot

x = azimuth
xlabel = "azimuth (°)"

fig, (ax1,ax2,ax3) = subplots(3,1, sharex=true, figsize=(5,6))
ax1.plot(x, GHI, "o-", label="horizontal (data)")
ax1.plot(x, Gplane, "D-", label="tilted (estimate)")
ax1.grid(true)
ax1.set(
    title="Irradiance transposition\non a plane tilted by $(slope)°",
    ylabel="Irradiance (W/m²)",
)

ax2.plot(x, rad2deg.(ω), "+-", color="tab:green", label="ω")
ax2.fill_between(x, rad2deg.(ωd), rad2deg.(ωu), alpha=0.2, ec="k", label="over ground")
ax2.fill_between(x, rad2deg.(ωd_plane), rad2deg.(ωu_plane), alpha=0.2, ec="k", label="over plane")
ax2.plot(x, rad2deg.(ωd_plane2), "+:", color="k", label="ω")
ax2.plot(x, rad2deg.(ωu_plane2), "+:", color="k", label="ω")
ax2.grid(true)

ax2.set(
    title="Sun position",
    ylabel="angle (°)",
)
ax2.legend()

ax3.plot(x, cosθZ, "+-", label="<cosθZ>")
ax3.plot(x, cosθ, "+-", label="<cosθ>")
ax3.grid(true)
dt_min = round(60*dt; digits=1)
ax3.set(
    title="Averaged cosine incidence angles (dt=$dt_min min)",
    xlabel=xlabel
)
ax3.legend()

ax1.legend()
fig.tight_layout()