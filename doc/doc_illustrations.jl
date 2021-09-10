# Illustrations for SolarIrradiance documentation
# Copyright (c) 2021 Pierre Haessig

using SolarIrradiance
using PyPlot
#pygui(true) # for interactive plot windows

## Load GHI data ##

# 2012, day 136: sun & clouds
year = 2012
n = 136
GHI_day = [0.0, 0.0, 0.0, 0.0, 0.0, 39.18, 209.0, 374.0, 240.0, 476.0, 379.0, 619.0, 872.01, 664.0, 680.0, 434.0, 311.0, 313.0, 169.0, 20.88, 0.0, 0.0, 0.0, 0.0]

# 2012, day 137: sun
n = 137
GHI_day = [0.0, 0.0, 0.0, 0.0, 0.0, 64.0, 218.0, 396.0, 566.0, 709.0, 827.01, 892.01, 851.01, 834.01, 756.01, 666.0, 519.0, 346.0, 168.0, 21.88, 0.0, 0.0, 0.0, 0.0]


## Custom colormaps for sets of license

# https://gka.github.io/palettes/#/9|s|0030dd,ef272d|ffffe0,ff005e,93003a|1|1
colorsBR = ["#0030dd", "#622cc7", "#8628b1", "#9f259b", "#b32486", "#c52370", "#d4235b", "#e22545", "#ef272d"]
cmapBR = matplotlib.colors.LinearSegmentedColormap.from_list("BlueRed", colorsBR)

# https://gka.github.io/palettes/#/9|s|3b2fa7,007300,887600,eb212a|ffffe0,ff005e,93003a|0|0
colorsBGR = ["#3b2fa7", "#254968", "#0f622a", "#117300", "#447400", "#777600", "#a1610b", "#c6411a", "#eb212a"]
cmapBGR = matplotlib.colors.LinearSegmentedColormap.from_list("BlueGreenRed", colorsBGR)

# https://gka.github.io/palettes/#/9|d|0030dd,4d4d4d|4d4d4d,df0822|0|0
colorsBKR = ["#0030dd", "#1337b9", "#263e95", "#3a4671", "#4d4d4d", "#713c42", "#962a37", "#ba192d", "#df0822"]
# https://gka.github.io/palettes/#/9|d|0024cf,505050|505050,df0822|0|1
colorsBKR = ["#0024cf", "#3d31ae", "#4d3d8e", "#52476f", "#505050", "#794c45", "#9c443a", "#be342e", "#df0822"]
cmapBKR = matplotlib.colors.LinearSegmentedColormap.from_list("BlueGrayRed", colorsBKR)

# https://gka.github.io/palettes/#/9|d|00aa00,b60053|b60053,0000aa|0|1
colorsGRB = ["#00aa00", "#659129", "#89743b", "#a25148", "#b60053", "#9e0069", "#82007e", "#5e0094", "#0000aa"]
cmapGRB = matplotlib.colors.LinearSegmentedColormap.from_list("GreenRedBlue", colorsGRB)


## Transpose GHI for a given slope & azimuth ##

# Time of the measurement
tc_range = 0:23
dt = 1/60 # 1 minute ~ instantaneous (0.0 is orbidden)

# Location (Rennes, France):
lat = 48.117 #°
lon = -1.678 #°

# Panel orientation:
slope = 40 #°
azimuth = 90 #° → west-facing
albedo = 0.5; # in [0,1]

Gplane = zeros(length(tc_range))

for (k, tck) in enumerate(tc_range)
    GHIk = GHI_day[k]
    Gplane[k] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
end

## Plot the single transpose output ##

fig, ax = subplots(figsize=(5,3.5))
ax.plot(tc_range, GHI_day, "o-", label="horizontal (data)")
ax.plot(tc_range, Gplane, "D-", label="tilted (estimate)")
ax.grid(true)
ax.set(
    title="Irradiance transposition\non a west-facing plane tilted by $(slope)°",
    xlabel="time (UTC hours)",
    ylabel="Irradiance (W/m²)",
)
ax.legend()
fig.tight_layout()
fname = "GHI_$(year)-$(n)_trans_s$(slope)_a$(azimuth).png"
fig.savefig(fname, dpi=300)

## Plot a range of azimuths ##

azi_range = -90:15:90

Gplane_azi = zeros(length(tc_range), length(azi_range))

for (i, azimuth) in enumerate(azi_range)
    for (k, tck) in enumerate(tc_range)
        GHIk = GHI_day[k]
        Gplane_azi[k,i] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
    end
end

# Plot azimuth range
fig, ax = subplots(figsize=(5,3.5))

ax.plot(tc_range, GHI_day, "-", label="GHI data",
        color=(0,0.7,0, 0.5), ms=3, lw=3, alpha=0.5)

for (i, azimuth) in enumerate(azi_range)
    lw=1
    ls="-"
    if azimuth == -90
        label = "−90° (E)"
    elseif azimuth == 0
        label = "0° (S)"
    elseif azimuth == +90
        label = "+90° (W)"
    else
        label = ""
        lw=0.5
        ls="-"
    end
    color = cmapBKR((i-1)/(length(azi_range)-1))
    ax.plot(tc_range, Gplane_azi[:,i], ls, ms=3,
            color = color, label=label, lw=lw)
end

ax.grid(true)
ax.set(
    title="Irradiance transposition\non a plane tilted by $(slope)° with varying azimuths",
    xlabel="time (UTC hours)",
    ylabel="Irradiance (W/m²)",
)
ax.set_ylim(0, 1400)
ax.legend(ncol=2)
fig.tight_layout()
fname = "GHI_$(year)-$(n)_trans_s$(slope)_azim_range.png"
fig.savefig(fname, dpi=300)


## Plot a range of slopes ##

azimuth = 0 #°
slo_range = 0:15:90

Gplane_slo = zeros(length(tc_range), length(slo_range))

for (i, slope) in enumerate(slo_range)
    for (k, tck) in enumerate(tc_range)
        GHIk = GHI_day[k]
        Gplane_slo[k,i] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
    end
end

# Plot slope range
fig, ax = subplots(figsize=(5,3.5))

ax.plot(tc_range, GHI_day, "-", label="GHI data",
        color=(0,0.7,0, 0.5), ms=3, lw=3, alpha=0.5)

for (i, slope) in enumerate(slo_range)
    lw=1
    ls="-"
    if slope == 0
        label = "0°"
    elseif slope == +45
        label = "45°"
    elseif slope == +90
        label = "90°"
    else
        label = ""
        lw=0.5
        ls="-"
    end
    color = cmapGRB((i-1)/(length(slo_range)-1))
    ax.plot(tc_range, Gplane_slo[:,i], ls, ms=3,
            color = color, label=label, lw=lw)
end

ax.grid(true)
ax.set(
    title="Irradiance transposition\non a south-facing plane with varying slopes",
    xlabel="time (UTC hours)",
    ylabel="Irradiance (W/m²)",
)
ax.set_ylim(0, 1400)
ax.legend(ncol=2, loc="upper center")
fig.tight_layout()
fname = "GHI_$(year)-$(n)_trans_a$(azimuth)_slope_range.png"
fig.savefig(fname, dpi=300)

## Animation of a range of azimuths ##

slope = 40 #°
azi_range = -90:3:90 # → 2s @ 30 fps

lineG_data = zeros(2, length(tc_range))
lineG_data[1,:] = tc_range

function update_Gplane_azi(i, lineG)
    azimuth = azi_range[i]
    # Compute Gplane
    for (k, tck) in enumerate(tc_range)
        GHIk = GHI_day[k]
        lineG_data[2,k] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
    end
    ax.set_title("Irradiance transposition on a plane\n tilted by $(slope)°, azimuth $(azimuth)°")
    # Update plot
    lineG.set_data(lineG_data)
    return (lineG, )
end

# Figure creation with null Gplane
fig, ax = subplots(figsize=(5,3.5))
ax.plot(tc_range, GHI_day, "o-", label="horizontal (data)")
lineG = ax.plot(tc_range, lineG_data[2,:], "D-", label="tilted (estimate)")[1]
ax.grid(true)
ax.set(
    title="Irradiance transposition on a plane\n tilted by $(slope)°, azimuth XX°",
    xlabel="time (UTC hours)",
    ylabel="Irradiance (W/m²)",
    ylim=(0,1300)
)
ax.legend(ncol=2, loc="upper center")
fig.tight_layout()

# Creating the Animation object
azim_ani = matplotlib.animation.FuncAnimation(
    fig, update_Gplane_azi, 1:length(azi_range), fargs=(lineG,), interval=33.33)

fname = "GHI_$(year)-$(n)_trans_s$(slope)_azim_range.gif"
writer = matplotlib.animation.PillowWriter(fps=30)
fname = "GHI_$(year)-$(n)_trans_s$(slope)_azim_range.mp4"
writer = matplotlib.animation.FFMpegWriter(fps=30, metadata=Dict(
        "title" => "GHI transposition on a plane tilted by $(slope)°, for a range of azimuths",
        "artist" => "Pierre Haessig",
        "date" => "2021"))
azim_ani.save(fname, dpi=216, writer=writer) # 5"*216 = 1080 px

## Show azimuth of the plane

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

function aziplot_setup(ax, plane_color = "tab:red", ra = 1.3)
    # N,E,S,W arrow
    NESW_color = "gray"
    ax.arrow(0, -ra , 0, 2ra, head_width=0.1, length_includes_head=true, color=NESW_color)
    ax.arrow(-ra , 0, 2ra, 0, head_width=0, color=NESW_color)
    ax.text(0, 1.1ra, "S\n0°", ha="center", va="bottom")
    ax.text(1.1ra, 0, "90° W", ha="left", va="center")
    ax.text(0, -1.1ra, "180°\nN", ha="center", va="top")
    ax.text(-1.1ra, 0, "E −90°", ha="right", va="center")

    # Plane orientation plot
    azimuth = 0 #°
    cosγ = cos(deg2rad(azimuth))
    sinγ = sin(deg2rad(azimuth))
    # plane
    plane_line = ax.plot([-cosγ, cosγ], [sinγ, -sinγ], color=plane_color, lw=3)[1]
    # normal arrow
    plane_arrow = ax.arrow(0, 0, sinγ, cosγ, head_width=0.1, length_includes_head=true, color=plane_color)
    plane_arrow_xy0 = plane_arrow.get_xy()

    # Axes setup:
    lim = 1.5ra
    ax.set(
        xlim=(-lim, lim),
        ylim=(-lim, lim),
        aspect = "equal"
    )
    ax.axis("off")

    return plane_arrow, plane_line, plane_arrow_xy0
end

"move plane to azimuth[i]"
function aziplot_update(i, plane_arrow, plane_line, plane_arrow_xy0)
    azimuth = azi_range[i]
    cosγ = cos(deg2rad(azimuth))
    sinγ = sin(deg2rad(azimuth))
    # plane
    plane_line.set_data([-cosγ, cosγ], [sinγ, -sinγ])
    # normal arrow
    rot = matplotlib.transforms.Affine2D().rotate_deg(-azimuth)
    xyr = rot.transform(plane_arrow_xy0)
    plane_arrow.set_xy(xyr)
    return plane_arrow, plane_line
end

plane_arrow, plane_line, plane_arrow_xy0 = aziplot_setup(ax)
aziplot_update(2, plane_arrow, plane_line, plane_arrow_xy0)


## Animation of a range of azimuth + azi_plot ##

slope = 40 #°
fps = 30
azi_range = -90:3:90 # → 2s @ 30 fps
print("azimuth anim $(azi_range[1])° to $(azi_range[end])°: $(length(azi_range)/fps)s @ $(fps) fps\n")

lineG_data = zeros(2, length(tc_range))
lineG_data[1,:] = tc_range

function update_Gplane_azi(i, lineG, ax)
    azimuth = azi_range[i]
    # Compute Gplane
    for (k, tck) in enumerate(tc_range)
        GHIk = GHI_day[k]
        lineG_data[2,k] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
    end
    ax.set_title("Irradiance transposition on a plane\n tilted by $(slope)°, azimuth $(azimuth)°")
    # Update plot
    lineG.set_data(lineG_data)
    return (lineG, )
end

function update_azianim(i, lineG, axG, plane_arrow, plane_line, plane_arrow_xy0)
    lineG = update_Gplane_azi(i, lineG, axG)[1]
    plane_arrow, plane_line = aziplot_update(i, plane_arrow, plane_line, plane_arrow_xy0)
    return (lineG, plane_arrow, plane_line)
end

# Figure creation with null Gplane
fig = figure(figsize=(9, 4.5))

gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

plane_arrow, plane_line = aziplot_setup(ax2)

# GHI plot setup
ax1.plot(tc_range, GHI_day, "o-", label="horizontal (data)")
lineG = ax1.plot(tc_range, lineG_data[2,:], "D-", label="tilted (estimate)")[1]
ax1.grid(true)
ax1.set(
    title="Irradiance transposition on a plane\n tilted by $(slope)°, azimuth XX°",
    xlabel="time (UTC hours)",
    ylabel="Irradiance (W/m²)",
    ylim=(0,1300)
)
ax1.legend(ncol=2, loc="upper center")
fig.tight_layout()

# Creating the Animation object
azim_ani = matplotlib.animation.FuncAnimation(
    fig, update_azianim, 1:length(azi_range),
    fargs=(lineG, ax1, plane_arrow, plane_line, plane_arrow_xy0),
    interval=33.33)

# fname = "GHI_$(year)-$(n)_trans_s$(slope)_azim_range2.gif"
# writer = matplotlib.animation.PillowWriter(fps=30)
fname = "GHI_$(year)-$(n)_trans_s$(slope)_azim_range2.mp4"
writer = matplotlib.animation.FFMpegWriter(fps=30, metadata=Dict(
        "title" => "GHI transposition on a plane tilted by $(slope)°, for a range of azimuths",
        "artist" => "Pierre Haessig",
        "date" => "2021"))
azim_ani.save(fname, dpi=240, writer=writer) # 4.5"*240 = 1080 px


## Animation of a range of slopes ##

azimuth = 0 #°
slo_range = 0:1.5:90 # → 2s @ 30 fps

lineG_data = zeros(2, length(tc_range))
lineG_data[1,:] = tc_range

function update_Gplane_slo(i, lineG)
    slope = slo_range[i]
    # Compute Gplane
    for (k, tck) in enumerate(tc_range)
        GHIk = GHI_day[k]
        lineG_data[2,k] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
    end
    ax.set_title("Irradiance transposition on a plane with azimuth $(azimuth)°,\n tilted by $(slope)°")
    # Update plot
    lineG.set_data(lineG_data)
    return (lineG, )
end

# Figure creation with null Gplane
fig, ax = subplots(figsize=(5,3.5))
ax.plot(tc_range, GHI_day, "o-", label="horizontal (data)")
lineG = ax.plot(tc_range, lineG_data[2,:], "D-", label="tilted (estimate)")[1]
ax.grid(true)
ax.set(
    title="Irradiance transposition on a plane with azimuth $(azimuth)°,\n tilted by XX°",
    xlabel="time (UTC hours)",
    ylabel="Irradiance (W/m²)",
    ylim=(0,1300)
)
ax.legend(ncol=2, loc="upper center")
fig.tight_layout()

# Creating the Animation object
azim_slo = matplotlib.animation.FuncAnimation(
    fig, update_Gplane_slo, 1:length(slo_range), fargs=(lineG,), interval=33.33)

fname = "GHI_$(year)-$(n)_trans_a$(azimuth)_slope_range.mp4"
azim_slo.save(fname, dpi=216, # 5"*216 = 1080 px
    metadata=Dict(
        "title" => "GHI transposition on a plane with azimuth $(azimuth)°, for a range of slopes",
        "artist" => "Pierre Haessig",
        "date" => "2021")
)