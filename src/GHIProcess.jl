# SolarIrradiance.jl
# Copyright (c) 2021 Pierre Galloo-Beauvais, Pierre Haessig
# This code is subject to the terms of the MIT license (see LICENSE.txt)

module GHIProcess

export GHI_download, GHI_read, sliceday
export forwardfill, backwardfill, centeredfill, interpolation
export GHIavg

using CSV
using DataFrames
using Dates
using HTTP

"""
    GHI_download(year, lat, lon, basename="GHI")

Download of GHI hourly data from PVGIS at location `lat`, `lon` (°) for `year`.

Data is saved in a CSV file with the `basename` prefix (default: `"GHI"`).

⚠️ `lat`, `lon` must be expressed with maximum 3 significant digits,
otherwise the PVGIS API refuses to generate the CSV containing the data

See also: [`GHI_read`](@ref) to load the data.
"""
function GHI_download(year, lat, lon, angle, aspect; basename="GHI" , version = 1)
     url = string("https://re.jrc.ec.europa.eu/api/v5_",version,"/seriescalc?",
        "lat=",lat, "&lon=",lon, "&angle=", angle, "&aspect=", aspect,
        "&startyear=", year, "&endyear=", year,
        "&outputformat=csv&browser=1")
    fname = string(basename, "_", lat, "_", lon, "_SA_", angle, "deg_", aspect, "aspect_",
                   year, "_", year, ".csv")
    
    io = open(fname, "w")
    r = HTTP.request("GET", url, response_stream=io)
    close(io)
    print("PVGIS GHI data downloaded and saved as $fname\n")
    return fname
end

"""
    GHI_read(source, year)

Load the CSV `source` containing PVGIS GHI hourly data for `year`.

Return a Vector of length 8760 or 8784 depending on whether `year` is a leap year.

See also: [`GHI_download`](@ref) to download the data.

# Example:

Download 2012 data for Rennes, France:

```jldoctest
julia> year = 2012; lon = -1.678; lat = 48.117;

julia> fname = GHI_download(year, lat, lon; basename="GHI")
PVGIS GHI data downloaded and saved as GHI_48.117_-1.678_SA_0deg_0deg_2012_2012.csv
"GHI_48.117_-1.678_SA_0deg_0deg_2012_2012.csv"

julia> GHI_year = GHI_read(fname, year)
8784-element Vector{Float64}:
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
 70.0
 97.0
  ⋮
 85.0
 49.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
```
"""
function GHI_read(source, year)
    nhours = 24 * (365 + Dates.isleapyear(year))
    df = CSV.read(source, DataFrame; header=9, select=[2], limit=nhours)
    GHI_year = df[!,1]
    return GHI_year
end

"""
    sliceday(GHI_year,n)

Renvoie une liste de 24 éléments, correspondant à la GHI d'un jour n d'une année. 

La sélection de l'année se fait via l'utilisation préalable de la fonction GHIannee(annee,phi_deg,lambda_deg), qui permet de récupérer la GHI pour toute une année. Il faut donc l'utiliser en amont, et stocker le résultat dans une liste, nommée ici "GHI_year".

#Exemple :

julia> GHI_year = GHIannee(2012,48.117,-1.678)
julia> sliceday(GHI_year,1)

24-element Vector{Any}:
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
 70.0
 97.0
 55.0
 68.0
 83.0
 58.0
 15.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0
  0.0


"""
function sliceday(GHI_year,n)
    GHI_day = GHI_year[1+(n-1)*24:n*24] 
    return GHI_day
end

"""
    forwardfill(lesgjour,heure)

La fonction forwardfill renvoie la valeur correspondant au GHI (moyen sur une minute) à une certaine heure, en réalisant une interpolation type bloqueur d'ordre 0 aux données d'entrées. La valeur de GHI à l'heure demandée correspond à la valeur de la mesure précédant directement cette heure.

La fonction prend en entrée une liste correspondant aux données de GHI d'une journée, sous forme d'une liste de 24 éléments, et une heure, qui peut être un flottant. 

Elle doit donc être précédée de l'utilisation de la fonction GHIannee et de la fonction sliceday.

#Exemple : 

julia>GHI_year=GHIannee(2012,48.117,-1.678)
julia>lesgjour=sliceday(GHI_year,1)
julia>forwardfill(lesgjour,12.5)

68.0


"""
function forwardfill(lesgjour,heure)
    
    if heure>=0 && heure<1+1/6
        return(lesgjour[1])
    elseif heure>=23+1/6 && heure<=24
        return(lesgjour[24])
    end

    for k in [1:22...]
        if heure>=k+1/6 && heure<k+1+1/6
            return(lesgjour[k+1])
        end
    end 
end

"""
    backwardfill(lesgjour,heure)

La fonction backwardfill renvoie la valeur correspondant au GHI (moyen sur une minute) à une certaine heure, en réalisant une interpolation type bloqueur d'ordre 0 inversé aux données d'entrées. La valeur de GHI à l'heure demandée correspond à la valeur de la mesure suivant directement cette heure.

La fonction prend en entrée une liste correspondant aux données de GHI d'une journée, sous forme d'une liste de 24 éléments, et une heure, qui peut être un flottant. 

Elle doit donc être précédée de l'utilisation de la fonction GHIannee et de la fonction sliceday.

#Exemple : 

julia>GHI_year=GHIannee(2012,48.117,-1.678)
julia>lesgjour=sliceday(GHI_year,1)
julia>backwardfill(lesgjour,12.5)

83.0


"""
function backwardfill(lesgjour,heure)
    
    if heure>=0 && heure<=1/6
        return(lesgjour[1])
    elseif heure>=22+1/6 && heure<=24
        return(lesgjour[24])
    end

    for k in [0:21...]
        if heure>k+1/6 && heure<=k+1+1/6
            return(lesgjour[k+2])
        end
    end
end

"""
    centeredfill(lesgjour,heure)

La fonction centeredfill renvoie la valeur correspondant au GHI (moyen sur une minute) à une certaine heure, en réalisant une interpolation type bloqueur centré sur la mesure. La valeur de GHI à l'heure demandée correspond à la valeur de la mesure la plus proche de cette heure.

La fonction prend en entrée une liste correspondant aux données de GHI d'une journée, sous forme d'une liste de 24 éléments, et une heure, qui peut être un flottant. 

Elle doit donc être précédée de l'utilisation de la fonction GHIannee et de la fonction sliceday.

#Exemple : 

julia>GHI_year=GHIannee(2012,48.117,-1.678)
julia>lesgjour=sliceday(GHI_year,1)
julia>centeredfill(lesgjour,12.5)

68.0

julia>GHI_year=GHIannee(2012,48.117,-1.678)
julia>lesgjour=sliceday(GHI_year,1)
julia>centeredfill(lesgjour,12.8)

83.0


"""
function centeredfill(lesgjour,heure)

    if heure>=0 && heure<4/6
        return(lesgjour[1])
    elseif heure>=23+4/6 && heure<=24
        return(lesgjour[24])
    end

    for k in [0:22...]
        if heure>=k+4/6 && heure<k+1+4/6
            return(lesgjour[k+2])
        end
    end
end

"""
    interpolation(lesgjour,heure)

La fonction interpolation renvoie la valeur correspondant au GHI (moyen sur une minute) à une certaine heure, en réalisant une interpolation linéaire sur les données d'entrées. 

La fonction prend en entrée une liste correspondant aux données de GHI d'une journée, sous forme d'une liste de 24 éléments, et une heure, qui peut être un flottant. 

Elle doit donc être précédée de l'utilisation de la fonction GHIannee et de la fonction sliceday.

#Exemple : 

julia>GHI_year=GHIannee(2012,48.117,-1.678)
julia>lesgjour=sliceday(GHI_year,1)
julia>interpolation(lesgjour,12.5)

77.50000000000001


"""
function interpolation(lesgjour,heure)
    if heure>=0 && heure<=1/6
        return(heure*lesgjour[1])
    elseif heure>23+1/6 && heure<=24
        return(lesgjour[24]+(heure-23-1/6)*(lesgjour[24]-lesgjour[23]))
    end 

    for k in [0:22...]
        if heure>k+1/6 && heure<=k+1+1/6
            return(lesgjour[k+1]+(heure-k-1/6)*(lesgjour[k+2]-lesgjour[k+1]))
        end
    end
end

"""
    GHIavg(GHI_year,n,tcini,dt,fonctionGHI)

Cette fonction évalue la valeur moyenne du GHI sur un intervalle de temps donné, un jour donné d'une année donnée, avec le choix d'une interpolation sur les données du jour. 

Cette fonction prend en entrée les données de GHI sur une année, obtenues via la fonction GHIannee.

tcini est l'heure initiale, le début de l'intervalle sur lequel on veut obtenir la valeur moyenne de GHI. Peut être un flottant. Ne peut pas dépasser 24.

dt est la durée de l'intervalle sur lequel on veut obtenir la valeur moyenne de GHI, en heures. Peut être un flottant. Ne peut pas dépasser 24.

fonctionGHI est une chaîne de caractères correspondant à la signature de la fonction d'interpolation choisie pour reconstruire les données de GHI. On peut donc utiliser les chaînes suivantes : "forwardfill", "backwardfill", "centeredfill", "interpolation".

# Exemple :

julia> GHI_year=GHIannee(2012,48.117,-1.678)
julia> GHIavg(GHI_year,1,12.5,0.5,interpolation)

73.0

julia> GHI_year=GHIannee(2012,48.117,-1.678)
julia> GHIavg(GHI_year,1,12.5,0.5,forwardfill)

68.0

"""
function GHIavg(GHI_year,n,tcini,dt,fonctionGHI)
    lesgjour=sliceday(GHI_year,n)
    
    G=0
    compteur=60*dt
    for k in range(tcini,length=convert(Int64,60*dt))
        G+=fonctionGHI(lesgjour,tcini)
    end
    return(G/compteur)
end


end # module