# %%
# M.Stargardt - 10.07.2025
import os
import pickle
import glob
import numpy as np
import geokit as gk


from reskit.default_paths import DEFAULT_PATHS
from reskit.parameters.parameters import OffshoreParameters


# %%


def waterDepthFromLocation(
    latitude,
    longitude,
    waterDepthFolderPath=None,
):
    """
    Returns the water depth (in meters) at a given geographic location.

    Parameters
    ----------
    latitude : float
        Latitude in decimal degrees.
    longitude : float
        Longitude in decimal degrees.
    waterDepthFolderPath : str, optional
        Path to the folder containing water depth .tif files. Loaded from defaults if not specified.

    Returns
    -------
    float or None
        Water depth at the specified location in meters (always positive). Returns None if not found.
    """

    if waterDepthFolderPath is None:
        waterDepthFolderPath = DEFAULT_PATHS.get("waterdepthFile")
        if waterDepthFolderPath is None:
            raise ValueError("No waterDepthFilePath is given. Please add it to default_path.yaml.")


    depthFiles = glob.glob(os.path.join(waterDepthFolderPath, "*.tif"))
    resultDepth = gk.raster.interpolateValues(source=[depthFiles],points=(longitude,latitude))

    return abs(resultDepth) if resultDepth is not None else None


# %% function to calculate the distance to the coastline
# if you want to execute the distance to coastline more often, please separete the loading of the taserband to increase execution time


def distanceToCoastline(latitude, longitude, path=None):
    """
    Computes the distance to the coastline from a given geographic point.

    Parameters
    ----------
    latitude : float
        Latitude in decimal degrees.
    longitude : float
        Longitude in decimal degrees.
    path : str, optional
        File path to the distance-to-coast raster. Loaded from defaults if not specified.

    Returns
    -------
    float or None
        Distance in kilometers, or None if the point is out of bounds or an error occurs.
    """
    if path is None:
        path = DEFAULT_PATHS.get("distancetoCoast")
        if path is None:
            raise ValueError("No distaneFilePath is given. Please add it to default_path.yaml.")
        
    try:
        value = gk.raster.interpolateValues(path, (longitude, latitude))

        return value

    except Exception as e:
        print(f"Error at Lat: {latitude}, Lon: {longitude}: {e}")
    return None


# %%
def calculateOffshoreCapex(
    baseCapex,
    capacity,
    hubHeight,
    waterDepth,
    coastDistance,
    rotorDiam,
    techYear=2050,
    shareTurb=0.449,
    shareFound=0.204,
    shareCable=0.181,
    shareOverhead=0.166,
    maxMonopileDepth=25,
    maxJacketDepth=55,
    baseDepth=17,
    baseDistCoast=27,
    baseCap=None,
    baseHubHeight=None,
    baseRotorDiam=None,
    defaultOffshoreParamsFp=None,
):
    """
    Scales a generic offshore CAPEX value based on water depth and distance to shore by taking capacity, hubheight and rotor diameter of a base case. If no base case is given, a default base case is applied.

    Parameters
    ----------
    baseCapex : float
        Reference CAPEX per kW (cost unit/kW) that should be scaled. base CApex must be given in €/kW to enable correct scaling.
    capacity : float
        Turbine rated capacity in kW.
    hubHeight : float
        Hub height in meters.
    waterDepth : float
        Site-specific water depth in meters.
    coastDistance : float
        Distance from site to nearest coast in kilometers.
    rotorDiam : float
        Rotor diameter in meters.
    techYear : int, optional
        Year of the applied technology, by default 2050.
    shareTurb : float, optional
        Share of turbine cost in total CAPEX in the baseline turbine reference case. Default is 0.449.
    shareFound : float, optional
        Share of foundation costin total CAPEX in the baseline turbine reference case. Default is 0.204.
    shareCable : float, optional
        Share of cable/connection cost in total CAPEX in the baseline turbine reference case. Default is0.181.
    shareOverhead : float, optional
        Share of overhead/miscellaneous costs in total CAPEX in the baseline turbine reference case. Default is 0.166.
    maxMonopileDepth : float, optional
        Maximum depth for monopile foundations, by default 25.
    maxJacketDepth : float, optional
        Maximum depth for jacket foundations, by default 55.
    baseDepth : float, optional
        Reference depth in CAPEX literature, by default 17.
    baseDistCoast : float, optional
        Reference coast distance, by default 27.
    baseCap : float, optional
        Reference turbine capacity. Loaded from CSV if not provided.
    baseHubHeight : float, optional
        Reference hub height. Loaded from CSV if not provided.
    baseRotorDiam : float, optional
        Reference rotor diameter. Loaded from CSV if not provided.
    defaultOffshoreParamsFp : str, optional
        Filepath to offshore turbine parameters CSV.

    Returns
    -------
    float
        Adjusted offshore wind CAPEX per kW for the given configuration. The cost unit is the same as the baseCapex.
    """
    assert np.isclose(
        shareTurb + shareFound + shareCable + shareOverhead, 1.0, rtol=1e-9
    ), "Sum of all cost shares must equal 1"

    assert (
        0 < maxMonopileDepth < 55
    ), "Maximum depth for monopile foundation must be between 0 and 55 m"

    assert (
        55 <= maxJacketDepth < 100
    ), "Maximum depth for jacket foundation must be between 55 and 100 m"

    assert (
        maxMonopileDepth < maxJacketDepth
    ), "Jacket depth must be greater than monopile depth"

    if any(_arg is None for _arg in [baseCap, baseHubHeight, baseRotorDiam, baseCapex]):
        params = OffshoreParameters(fp=defaultOffshoreParamsFp, year=techYear)
    elif not all(_arg is None for _arg in [defaultOffshoreParamsFp, techYear]):
        raise ValueError(
            "techYear and defaultOffshoreParamsFp are expected to be None if "
            "baseCap, baseHubHeight, baseRotorDiam and baseCapex are provided explicitly."
        )

    if baseCap is None:
        baseCap = params.base_capacity
        print("baseCap is taken from overall techno-economic file")
    if baseHubHeight is None:
        baseHubHeight = params.base_hub_height
        print("baseHubHeight is taken from overall techno-economic file")
    if baseRotorDiam is None:
        baseRotorDiam = params.base_rotor_diam
        print("baseRotorDiam is taken from overall techno-economic file")
    if baseCapex is None:
        baseCapex = params.base_capex_per_capacity
        print("inputCapex is taken from overall techno-economic file")


    turbineCostBase = baseCapex * shareTurb
    foundCostBase = baseCapex * shareFound
    cableCostBase = baseCapex * shareCable
    overheadCostBase = baseCapex * shareOverhead

    # Scale turbine cost
    turbineCostNew = onshoreTcc(
        capacity,
        hubHeight,
        rotorDiam,
        gdpEscalator=1,
        bladeMaterialEscalator=1,
        blades=3,
    )
    turbineCostReference = onshoreTcc(
        baseCap,
        baseHubHeight,
        baseRotorDiam,
        gdpEscalator=1,
        bladeMaterialEscalator=1,
        blades=3,
    )
    costRatioTurbine = turbineCostNew / turbineCostReference
    newTurbineCost = turbineCostBase * costRatioTurbine

    # Scale foundation cost
    depthBaseCost = getRatedCostFromWaterDepth(
        baseDepth, maxMonopileDepth, maxJacketDepth
    )
    depthPlantCost = getRatedCostFromWaterDepth(
        waterDepth, maxMonopileDepth, maxJacketDepth
    )
    costRatioFoundation = depthPlantCost / depthBaseCost
    newFoundationCost = foundCostBase * costRatioFoundation

    # Scale cable cost
    cableRatio = getCableCost(coastDistance, capacity) / getCableCost(
        baseDistCoast, baseCap
    )
    newCableCost = cableCostBase * cableRatio

    # Combine all costs
    totalOffshoreCapex = (
        newTurbineCost + newFoundationCost + newCableCost + overheadCostBase
    )

    return totalOffshoreCapex




# %%
def getRatedCostFromWaterDepth(depth, maxMonopileDepth=25, maxJacketDepth=55):
    """
    Estimates the rated cost of offshore wind turbine foundations based on water depth.

    Parameters
    ----------
    depth : float
        Water depth at the installation site (in meters).
    maxMonopileDepth : float, optional
        Threshold depth for monopile foundations, by default 25.
    maxJacketDepth : float, optional
        Threshold depth for jacket foundations, by default 55.

    Returns
    -------
    float
        Rated cost in €_2023/kW.

    References
    ----------
    Rogeau et al. (2023), Renewable and Sustainable Energy Reviews.
    """
    depth = abs(depth)

    if depth < maxMonopileDepth:
        c1, c2, c3 = 181, 552, 370
    elif depth <= maxJacketDepth:
        c1, c2, c3 = 103, -2043, 478
    else:
        c1, c2, c3 = 0, 697, 1223

    return c1 * depth**2 + c2 * depth + c3 * 1000


# %%
def getCableCost(distance, capacity, variableCostFactor=1350, fixedCost=0):
    """
    Calculates the cost for connecting an offshore wind power plant to the coastline.

    Parameters
    ----------
    distance : float
        Distance to coastline in kilometers.
    capacity : float
        Power plant's capacity in kW.
    variableCostFactor : float, optional
        Cost multiplier in €/kW/km, by default 1 350.
    fixedCost : float, optional
        Fixed connection cost in the respective currency unit. Defaults to 0 [€]

    Returns
    -------
    float
        Total cable connection cost in monetary units.

    References
    ----------
    Rogeau et al. (2023), "Review and modeling of offshore wind CAPEX",
    Renewable and Sustainable Energy Reviews, DOI: 10.1016/j.rser.2023.113699
    """
    assert distance > 0, "distance must be larger tan 0"
    assert capacity > 0, " turbine capacity must be larger than 0"
    assert variableCostFactor > 0, "cost factor must be larger tan 0"
    assert fixedCost >= 0, "fixed Cost must be postive or 0"

    

    variableCost = variableCostFactor * distance * capacity
    cableCost = fixedCost + variableCost

    return cableCost


def onshoreTcc(cp, hh, rd, gdpEscalator=None, bladeMaterialEscalator=None, blades=None):
    """
    Calculates the turbine capital cost (TCC) of a 3-blade onshore wind turbine.

    Parameters
    ----------
    cp : float
        Turbine capacity in kW.
    hh : float
        Hub height in meters.
    rd : float
        Rotor diameter in meters.
    gdpEscalator : float, optional
        Labor cost escalator, by default taken from OffshoreParameters.
    bladeMaterialEscalator : float, optional
        Blade material cost escalator, by default taken from OffshoreParameters.
    blades : int, optional
        Number of blades, by default taken from OffshoreParameters.

    Returns
    -------
    float
        Turbine capital cost (TCC) in monetary units.

    References
    ----------
    Fingersh et al. (2006), NREL. https://www.nrel.gov/docs/fy07osti/40566.pdf
    """
    if gdpEscalator is None or bladeMaterialEscalator is None or blades is None:
        offshoreParams = OffshoreParameters()
        gdpEscalator = offshoreParams.gdp_escalator
        bladeMaterialEscalator = offshoreParams.blade_material_escalator
        blades = offshoreParams.blades

    rr = rd / 2
    sa = np.pi * rr * rr

    singleBladeMass = 0.4948 * np.power(rr, 2.53)
    singleBladeCost = (
        (0.4019 * np.power(rr, 3) - 21051) * bladeMaterialEscalator
        + 2.7445 * np.power(rr, 2.5025) * gdpEscalator
    ) * (1 - 0.28)

    hubMass = 0.945 * singleBladeMass + 5680.3
    hubCost = hubMass * 4.25

    pitchSystemCost = 2.28 * (0.2106 * np.power(rd, 2.6578))
    noseConeMass = 18.5 * rd - 520.5
    noseConeCost = noseConeMass * 5.57

    lowSpeedShaftCost = 0.01 * np.power(rd, 2.887)
    bearingMass = (rd * 8 / 600 - 0.033) * 0.0092 * np.power(rd, 2.5)
    bearingCost = 2 * bearingMass * 17.6

    breakCouplingCost = 1.9894 * cp - 0.1141
    generatorCost = cp * 219.33
    electronicsCost = cp * 79
    yawSystemCost = 2 * (0.0339 * np.power(rd, 2.964))
    mainframeMass = 1.228 * np.power(rd, 1.953)
    mainframeCost = 627.28 * np.power(rd, 0.85)
    platformAndRailingMass = 0.125 * mainframeMass
    platformAndRailingCost = platformAndRailingMass * 8.7

    electricalConnectionCost = cp * 40
    hydraulicAndCoolingSystemCost = cp * 12
    nacelleCost = 11.537 * cp + 3849.7
    towerMass = 0.2694 * sa * hh + 1779
    towerCost = towerMass * 1.5

    turbineCapitalCost = (
        singleBladeCost * blades
        + hubCost
        + pitchSystemCost
        + noseConeCost
        + lowSpeedShaftCost
        + bearingCost
        + breakCouplingCost
        + generatorCost
        + electronicsCost
        + yawSystemCost
        + mainframeCost
        + platformAndRailingCost
        + electricalConnectionCost
        + hydraulicAndCoolingSystemCost
        + nacelleCost
        + towerCost
    )

    return turbineCapitalCost
