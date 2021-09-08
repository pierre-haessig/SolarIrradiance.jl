# SolarIrradiance.jl

The SolarIrradiance Julia package implements computations about solar radiation
on the surface of Earth.
Applications include thermal modeling of buildings and Solar

In particular, it implements several methods from the litterature
to transpose global irradiance data measured on an *horizontal* plane (GHI)
to an irradiance on an *inclined* plane (like a solar panel or the wall of a building).

Authors:

- Pierre Galloo-Beauvais, ENS Rennes
- Pierre Haessig, CentraleSupélec, IETR

This code is made available under a MIT license (see [LICENSE.txt](LICENSE.txt)).

## Usage

### Installation

To install the package, enter the [Pkg REPL](https://pkgdocs.julialang.org/v1/getting-started/)
(pressing `]`) and use the `add` command with the adress of the repository:
```
(@v1.6) pkg> add https://github.com/pierre-haessig/SolarIrradiance.jl
```

Now, it can be used with
```
julia> using SolarIrradiance
```

### Example: Transposing horizontal irradiance data on a tilted plane


Assuming we have the following vector of Global Horizontal Irradiance (GHI) data:
```
julia> GHI_day = [0.0, 0.0, 0.0, 0.0, 0.0, 39.18, 209.0, 374.0, 240.0, 476.0, 379.0, 619.0, 872.01, 664.0, 680.0, 434.0, 311.0, 313.0, 169.0, 20.88, 0.0, 0.0, 0.0, 0.0];
```

You can transform it to an irradiance on an inclined panel with the `global_radiation_tilt` function.
It requires inputs which specify both the GHI measurement to be transposed, the time and location of this measurement, and finally the orientation of the plane for the transposition. These parameters can be set as:

```julia
# Time of the measurement
n = 136 # can be computed with dayofyear(2012,05,15)
tc = 0:23
dt = 1/60 # 1 minute ~ instantaneous (0.0 is orbidden)

# Location (Rennes, France):
lat = 48.117 #°
lon = -1.678 #°

# Panel orientation:
slope = 40 #°
azimuth = 90 #° → west-facing
albedo = 0.5; # in [0,1]
```

Then, we can tranpose each of the 24 values in a loop:

```
Gplane = zeros(length(tc))

for (k, tck) in enumerate(tc)
    GHIk = GHI_day[k]
    Gplane[k] = global_radiation_tilt(GHIk, n, tck, dt, lat, lon, slope, azimuth, albedo)
end
```

For a more detailed introduction, see [doc/Usage.ipynb](doc/Usage.ipynb).