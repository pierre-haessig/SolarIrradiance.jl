# SolarIrradiance.jl
# Copyright (c) 2021 Pierre Galloo-Beauvais, Pierre Haessig
# This code is subject to the terms of the MIT license (see LICENSE.txt)

module SolarIrradiance

export dayofyear, declination, EoT, hour_angle, day_bounds
export cosθavg, cosθZavg, Rb, global_radiation_tilt

include("GHIProcess.jl")
export GHIProcess

using Dates


"""
    declination(n)

Compute the solar declination angle δ of the day `n`, in radians.

δ is the angle between the Sun-Earth direction and the equatorial plane.
It varies between −23.45° on December 21-22 and +23.45° on June 22
and is zero at the equinoxes.
See https://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle.

# Examples

Declination at the summer solstice in the northern hemisphere:
```jldoctest
julia> using Dates

julia> n = dayofyear(2021,06,22)
173

julia> δ = declination(n)
0.40924559967108437

julia> round(rad2deg(δ); digits=3)
23.448
```

Declination at March equinoxe:
```jldoctest
julia> n = dayofyear(2021,03,22)
81

julia> round(declination(n); digits=3)
-0.0
```

"""
declination(n) = deg2rad(23.45) * sin(2π*(284+n)/365)


"""
    EoT(n)

Equation of Time on day `n`, in hours.

The Equation of Time is the difference between the apparent (or true) solar time
and the mean solar time.
See https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time

# Examples

Minimal value is about -14 minutes in mid February:
```jldoctest
julia> EoT(45)*60
-14.26116933144462
```

Maximal value is about +16 minutes in mid November:
```jldoctest
julia> EoT(306)*60
16.388668719846564
```
"""
function EoT(n)
    B= 2π*(n-1)/365
    return 3.82*(0.000075+0.001868*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B))
end


"""
    hour_angle(n, tc, lon)

Compute the hour angle of the Sun ω (in rad) a civil hour `tc` on day `n`, at longitude `lon` (°). 

See also the inverse function `time_hour_angle`.

L'angle horaire du Soleil donne un renseignement sur le temps solaire vrai, il vaut 0 quand le Soleil est à son point culminant de la journée, est négatif avant ce point et positif après ce point.

TODO: add missing time zone input
"""
hour_angle(n, tc, lon) = (tc + lon/15 + EoT(n) - 12)*deg2rad(15)

"""
    time_hour_angle(n, ω, lon)

Compute the civil hour `tc` of an hour angle of the Sun `ω` (in rad) on day `n`,
at longitude `lon` (°). 

It is the inverse function of `tc -> hour_angle(n, tc, lon)`.

TODO: add missing time zone input
"""
time_hour_angle(n, ω, lon) = ω/deg2rad(15) - EoT(n) + 12 - lon/15

"""
    day_bounds(n, slope, azimuth, lat)

Compute the hour angle bounds when the sun is above a panel during day `n`.

Panel has a `slope` (in [0, 90°]) and `azimuth` (in [-90°, 90°])
and is located at latitude `lat` (°).

Return `ωu`, `ωd` (in rad) which are respectively the hour angles of 
sunrise and sunset *over the panel*.
When `slope` is 0, this corresponds to regulat sunrise and sunset (over the ground).
and there is a dedicated method for this, see `day_bounds(n, lat)`.

Mathematically, with θ the angle of incidence of the Sun on the panel,
this means that cosθ=0 at hour angle `ωu` and `ωd` and cosθ ≥ 0 in between.
"""
function day_bounds(n, slope, azimuth, lat)
    beta=deg2rad(slope)
    gamma=deg2rad(azimuth)
    phi=deg2rad(lat)
    delta=declination(n)

    sinδ=sin(delta)
    cosδ=cos(delta)
    sinϕ=sin(phi) 
    cosϕ=cos(phi)
    cosβ=cos(beta)
    sinβ=sin(beta)
    cosγ=cos(gamma)
    
    A = sinδ*sinϕ*cosβ - sinδ*cosϕ*sinβ*cosγ
    B = cosδ*cosϕ*cosβ + cosδ*sinϕ*sinβ*cosγ
    C = cosδ*sinβ*sin(gamma)

    D = C^2+B^2-A^2
    E = A-B

    ωu = 2*atan((-C+sqrt(D))/E)
    ωd = 2*atan((-C-sqrt(D))/E)
    
    return ωu, ωd
end

"""
    day_bounds(n, lat)

Compute the hour angle bounds when the sun is above the horizon during day `n`
at latitude `lat` (°).

Return `ωu`, `ωd` (in rad) which are respectively the hour angles of sunrise and sunset.

Mathematically, with θz the solar zenith angle (angle of incidence of the Sun on the ground),
this means that cosθz=0 at hour angle `ωu` and `ωd` and cosθz ≥ 0 in between.
"""
function day_bounds(n, lat)
    phi=deg2rad(lat)
    
    # Sunset
    ωd = acos(-tan(phi)*tan(declination(n)))
    # Sunrise
    ωu = -ωd
    
    return ωu, ωd
end


"""
    cosθavg(n, tcini, dt, slope, azimuth, lat, lon)

Cosine of θ, the angle of incidence of the sun on a tilted panel.

Panel has orientation (`slope`, `azimuth` °) and is located at (`lon`, `lat` °).
Value is computed on day `n`, averaged on a time interval `dt`
starting at civil hour `tcini`.

See also: [`cosθZavg`](@ref) for the incidence on the _ground_.
"""
function cosθavg(n, tcini, dt, slope, azimuth, lat, lon)
    beta=deg2rad(slope)
    gamma=deg2rad(azimuth)
    phi=deg2rad(lat)
    
    ω1=hour_angle(n, tcini, lon)
    ω2=hour_angle(n, tcini+dt, lon)
    
    ωup, ωdown = day_bounds(n, slope, azimuth, lat)
    
    if ω1<ωup && ω2<ωup
        return 0.0
    elseif ω1>ωdown && ω2>ωdown
        return 0.0
    elseif ω1<ωup && ω2>ωup
        ω1 = ωup
        tdebut = time_hour_angle(n, ωup, lon)
        tfin = tcini+dt
        dt = tfin-tdebut # TODO: check that the change of `dt` is correct
    elseif ω1<ωdown && ω2>ωdown
        ω2 = ωdown
        tdebut = tcini
        tfin = time_hour_angle(n, ωdown, lon)
        dt = tfin-tdebut # TODO: check that the change of `dt` is correct
    end
    
    delta = declination(n)

    sinδ = sin(delta)
    cosδ = cos(delta)
    sinϕ = sin(phi)
    cosϕ = cos(phi)
    cosβ = cos(beta)
    sinβ = sin(beta)
    cosγ = cos(gamma)
    
    a = sinδ*sinϕ*cosβ - sinδ*cosϕ*sinβ*cosγ
    b = cosδ*cosϕ*cosβ + cosδ*sinϕ*sinβ*cosγ
    c = cosδ*sinβ*sin(gamma)
    f = ω2 - ω1
    g = sin(ω2) - sin(ω1)
    h = cos(ω2) - cos(ω1)
    x = (a*f + b*g - c*h)*12/(pi*dt)
    
    return x
end

"""
    cosθZavg(n, tcini, dt, lat, lon)

Cosine of ``θ_Z``, the solar zenith angle, i.e. the angle of incidence of the sun on the ground.

Location is (`lon`, `lat` °).
Value is computed on day `n`, averaged on a time interval `dt`
starting at civil hour `tcini`.

See also: [`cosθavg`](@ref) for the incidence on a _tilted panel_.
"""
function cosθZavg(n, tcini, dt, lat, lon)
    phi=deg2rad(lat)
    
    ω1=hour_angle(n, tcini, lon)
    ω2=hour_angle(n, tcini+dt, lon)
    
    ωup, ωdown = day_bounds(n, lat)
    
    if ω1<ωup && ω2<ωup
        return 0.0
    elseif ω1>ωdown && ω2>ωdown
        return 0.0
    elseif ω1<ωup && ω2>ωup
        ω1 = ωup
        tdebut = time_hour_angle(n, ωup, lon)
        tfin = tcini+dt
        dt = tfin-tdebut # TODO: check that the change of `dt` is correct
    elseif ω1<ωdown && ω2>ωdown
        ω2 = ωdown
        tdebut = tcini
        tfin = time_hour_angle(n, ωdown, lon)
        dt = tfin-tdebut # TODO: check that the change of `dt` is correct
    end
    
    delta = declination(n)
    
    d = cos(phi)*cos(delta)
    e = sin(phi)*sin(delta)
    f = ω2 - ω1
    g = sin(ω2) - sin(ω1)
    x = (d*g + e*f)*12/(pi*dt)
    
    return x
end

"""
    Rb(cosθ, cosθZ, Rbsat=5)

Rb = cosθ/cosθZ with special care of zeros and saturation at `Rbsat`
"""
function Rb(cosθ, cosθZ, Rbsat=5.)
    if cosθ==0. || cosθZ==0.
        return 0.0
    else
        x = cosθ/cosθZ
    end
    
    if x>Rbsat
        return Rbsat
    else
        return x
    end
end

"""
    Rb(n, tcini, dt, slope, azimuth, lat, lon, Rbsat=5)

Cette fonction calcule le terme de redressement de GHI direct pour un panneau PV incliné d'un angle `slope`, orienté d'un angle `azimuth`, situé à une latitude `lat` et une longitude `lon` (°). Elle le calcule pour un jour n donné, sur un intervalle donné de durée dt à partir d'un instant initial tcini.

⚠️ Fonction avec encore quelques problèmes. Des pics apparaissent dans certaines conditions, une saturation a été mise en place afin que ces pics n'impactent pas les calculs utilisant cette fonction. Cette saturation, choisie empiriquement, empêche le terme de redressement de dépasser la valeur de 5.
"""
function Rb(n, tcini, dt, slope, azimuth, lat, lon, Rbsat=5.)
    cosθ = cosθavg(n, tcini, dt, slope, azimuth, lat, lon)
    cosθZ = cosθZavg(n, tcini, dt, lat, lon)
    
    return Rb(cosθ, cosθZ, Rbsat)
end

"""
    global_radiation_tilt(GHI::Real, n, tcini, dt, lat, lon, slope, azimuth, albedo;
                          diffuse_model=:EKD82, transpose_model=:HDKR)

Estimate the global irradiance on a tilted panel from horizontal data `GHI`.

Panel has a orientation (`slope`, `azimuth`) and is located at (`lon`, `lat`).
`GHI` data is for day `n`, averaged on a time interval `dt`,
starting at civil hour `tcini` (UTC).

Reflection from surrounding area is modeled by `albedo` in [0,1].

# Empirical models

Many steps of the irradiance transposition from a horizontal plane to a tilted plane
are deduced from geometry. However, two aspects of the estimation are _empirical_,
so that several models are available:
1. splitting the global horizontal irradiance between direct and diffuse components
2. transposing the diffuse irradiance on a tilted plane

## 1. Diffuse fraction estimation

`diffuse_model` sets the model used to estimate the _fraction of diffuse radiation_
from the clearness index:
- `:EKD82`: Erbs, Klein & Duffie (1982). Estimation of the diffuse radiation fraction for
  hourly, daily and monthly-average global radiation. Solar Energy 28 (4), 293–302.
- `:OH77`: Orgill & Hollands (1977). Correlation equation for hourly diffuse radiation
  on a horizontal surface. Sol. Energy 19, 357–359.
- `:CR79`: Collares-Pereira & Rabl (1979), The average distribution of solar radiation
  —correlations between diffuse and hemispherical and between daily and hourly insolation values.
  Solar Energy 22, 155.

## 2. Diffuse irradiance transposition

`transpose_model` sets the model to estimate the diffuse radiation on a tilted plane.
Models differ on their assumption on how the diffuse irradiance is spatially distributed
on the sky hemisphere.
See (Demain 2013) for an introduction to those models and Duffie & Beckman's book for HDKR.
Models are classified between isotropic and anisotropic.

Isotropic models:
- `:LJ62`: Liu & Jordan (1962). Daily insolation on surfaces tilted towards the equator. ASHRAE, 53:526-41.
- `:Ko86`: Korokanis (1986). On the choice of the angle of tilt for south facing solar collectors in the Athens basin area. Solar Energy, 36:217-25.
- `:Ba02`: Badescu (2002). 3D isotropic approximation for solar diffuse irradiance on tilted
surfaces. Renewable Energy, 26:221-3.

Anisotropic models:
- `:Wi82`: Willmot (1982). On the climatic optimization of the tilt and azimuth of
  flat-plate solar collectors. Solar Energy, 28:205-16.
- `:Bu77`: Bugler (1977). The determination of hourly insolation on an inclined plane
  using a diffuse irradiance model based on hourly measured global horizontal insolation.
  Solar Energy, 19:477-91.
- `:MI83`: Ma & Iqbal (1983) Statistical comparison of models for estimating solar radiation
  on inclined surfaces. Solar Energy 31 (3), 313-317.
- `:Ha79`:  Hay (1979). Study of shortwave radiation on non-horizontal surfaces.
  Canadian Climate Centre, Report No. 79-12, Downsview, Ontario.
- `:TC77`:  Temps & Coulson (1977). Solar radiation incident upon slopes of different orientation.
  Solar Energy, 19:179-84.
- `:Kl79`: Klucher (1979). Evaluation of models to predict insolation on tilted surfaces.
  Solar Energy, 23:111-114.
- `:HDKR`: HDKR model (as named in Duffie & Beckman's book) used by HOMER, proposed by
  Reindl, Beckman & Duffie (1990). Evaluation of Hourly Tilted Surface Radiation Models.
  Solar Energy, 45, 9-17.
"""
function global_radiation_tilt(GHI::Real, n, tcini, dt, lat, lon, slope, azimuth, albedo;
                               diffuse_model=:EKD82, transpose_model=:HDKR)
    delta = declination(n)
    phi = deg2rad(lat)
    beta = deg2rad(slope)
    gamma = deg2rad(azimuth)
    
    # Extraterrestrial normal irradiance
    Gon = 1367*(1+0.033*cos(deg2rad((360*n)/365)))
    
    # Averaged horizontal extraterrestrial irradiance
    ω1 = hour_angle(n, tcini, lon)
    ω2 = hour_angle(n, tcini+dt, lon)
    
    # TODO: use the cosθZavg function here instead:
    # Goavg = Gon*cosθZavg, I believe
    a = cos(phi)*cos(delta) * (sin(ω2) - sin(ω1))
    b = sin(phi)*sin(delta) * (ω2 - ω1)
    Goavg = Gon*(a+b)*12/(pi*dt)
    if Goavg < 0
        Goavg = 0.0
    end
    
    # Clearness index `kt`
    if Goavg > 0
        kt = GHI/Goavg
            if kt > 1
                kt = 1.0
            end
    else
        kt = 0.0
    end
    
    # Diffuse irradiance, estimated from `kt` and GHI using various models
    if diffuse_model == :EKD82 # Erbs et al. 1982
        if kt<=0.22
            GHIdiffus = (1-0.09*kt)*GHI
        elseif kt>0.8
            GHIdiffus = 0.165*GHI
        else
            GHIdiffus = (0.9511-0.1604*kt+4.388*kt^2-16.638*kt^3+12.336*kt^4)*GHI
        end
        
    elseif diffuse_model == :OH77 # Orgill & Hollands 1977
        if kt<=0.35
            GHIdiffus = (1-0.249*kt)*GHI
        elseif kt>0.75
            GHIdiffus = 0.177*GHI
        else
            GHIdiffus = (1.557-1.84*kt)*GHI
        end
        
    elseif diffuse_model == :CR79 # Colarres-Pereira & Rabl 1979
        if kt<=0.17
            GHIdiffus = 0.99*GHI
        elseif kt>=0.80
            GHIdiffus = 0.2*GHI
        elseif kt>0.17 && kt<0.75
            GHIdiffus = (1.188-2.272*kt+9.473*kt^2-21.865*kt^3+14.648*kt^4)*GHI
        else # [0.75, .80[ (not in original publication)
            GHIdiffus = (0.632-0.54*kt)*GHI
        end
    else
        error("unknown diffuse fraction estimation model $diffuse_model")
    end 
    
    # Beam (direct) irradiance:
    GHIdirect=GHI-GHIdiffus
    
    # Direct (beam) radiation on tilted surface
    # averaged cos(θ) and cos(θz)
    cosθ = cosθavg(n, tcini, dt, slope, azimuth, lat, lon)
    cosθZ = cosθZavg(n, tcini, dt, lat, lon)
    RB=Rb(cosθ, cosθZ)
    Gdirect=RB*GHIdirect
    
    # Diffuse radiation on tilted surface
    if transpose_model == :LJ62 # Liu & Jordan (1962)
        Gdiffus=(1+cos(beta))*GHIdiffus/2
        
    elseif transpose_model == :Ko86 # Korokanis (1986)
        Gdiffus=(2+cos(beta))*GHIdiffus/3
        
    elseif transpose_model == :Ba02 # Badescu (2002)
        Gdiffus=(3+cos(2*beta))*GHIdiffus/4
        
    elseif transpose_model == :Wi82 # Willmot (1982)
        if cosθ <= 0
            AW = 0.0
            BW = 1.0115-0.20293*beta-0.080823*beta^2
        else
            AW = Gdirect*RB/(1367*cosθ)
            BW = (1.0115-0.20293*beta-0.080823*beta^2)*(1-(Gdirect/(1367*cosθ)))
        end
        Gdiffus = (AW+BW) * GHIdiffus
        
    elseif transpose_model == :Bu77 # Bugler (1977)
        AB = (1+cos(beta))/2
        if GHIdiffus <= 0
            BB = 0.0
        else
            BB = 0.05*(Gdirect/GHIdiffus)
        end
        if cosθZ <= 0
            CB = cosθ
        else
            CB = cosθ - (1+cos(beta))/(2*cosθZ)
        end  
        Gdiffus = (AB + BB*CB) * GHIdiffus 
        
    elseif transpose_model == :MI83 # Ma & Iqbal (1983)
        AM = (1-kt)*(1+cos(beta))/2
        BM = kt*RB
        Gdiffus = (AM+BM)*GHIdiffus
        
    elseif transpose_model == :Ha79 # Hay (1979)
        if Goavg <= 0
            AH = 1.0
            CH = 0.0
        else
            AH = 1 - (GHIdirect/Goavg)
            CH = RB * (GHIdirect/Goavg)
        end
        BH = (1+cos(beta))/2
        Gdiffus = (AH*BH+CH) * GHIdiffus
        
    elseif transpose_model == :TC77 # Temps & Coulson (1977)
        ATC = (1+cos(beta))/2
        BTC = 1+(cosθ^2)*(sqrt(1-cosθZ^2))^3
        CTC = 1+(sin(beta/2))^3
        Gdiffus = ATC*BTC*CTC*GHIdiffus
        
    elseif transpose_model == :Kl79 # Klucher (1979)
        AK = (1+cos(beta))/2
        if GHI <= 0
            FK = 1.0
        else
            FK = 1 - (GHIdiffus/GHI)^2
        end
        BK = 1 + FK*(cosθ^2)*(sqrt(1-cosθZ^2))^3
        CK = 1 + FK*(sin(beta/2))^3
        Gdiffus = AK*BK*CK * GHIdiffus
        
    elseif transpose_model == :HDKR # HDKR model from Reindl, Beckman & Duffie (1990)
        if Goavg > 0
            Ai = GHIdirect/Goavg
        else
            Ai = 0.0
        end
        if GHI > 0
            f = sqrt(GHIdirect/GHI)
        else
            f = 0.0
        end
        AHDKR = Ai*RB
        BHDKR = 1 - Ai
        CHDKR = (1+cos(beta))/2
        DHDKR = 1 + f*(sin(beta/2))^3
        Gdiffus = (AHDKR + BHDKR*CHDKR*DHDKR) * GHIdiffus
    
    else
        error("unknown diffuse irradiance transposition model $transpose_model")
    end
    
    # Reflected irradiance
    Greflechi = GHI * albedo * (1-cos(beta))/2
    
    # Total irradiance
    Gtot = Gdiffus + Gdirect + Greflechi
    
    # Saturate output
    if Gtot>2000. 
        # TODO: transform fixed threshold into parameter
        Gtot = 2000.
    elseif Gtot<0
        Gtot = 0.
    end
    
    return Gtot
end

"""
    global_radiation_tilt(n, tcini, dt, lat, lon, fonctionGHI, diffuse_model, transpose_model, GHI_data, slope, azimuth, albedo)

Same as first `global_radiation_tilt(GHI::Real...)`, but 

`GHI_data` is the Vector of hourly GHI of a given year, obtained from function `GHIannee`.

Elle prend également en argument une chaîne de caractères `fonctionGHI`, qui permet de choisir l'interpolation des données du GHI. Cette chaîne de caractère peut donc être : "forwardfill", "backwardfill", "centeredfill" et "interpolation".
"""
function global_radiation_tilt(GHI_data::Vector, fonctionGHI, n, tcini, dt, lat, lon, slope, azimuth, albedo;
                               diffuse_model=:EKD82, transpose_model=:HDKR)
    GHI = GHIProcess.GHIavg(GHI_data, n, tcini, dt, fonctionGHI)
    Gtot = global_radiation_tilt(GHI, n, tcini, dt, lat, lon, slope, azimuth, albedo;
                                 diffuse_model=:EKD82, transpose_model=:HDKR)
    return Gtot
end

end # module
