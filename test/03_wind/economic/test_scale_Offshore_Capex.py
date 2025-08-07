from reskit.wind.economic.scale_Offshore_Capex import *
import numpy as np
import pytest
from reskit.default_paths import DEFAULT_PATHS
from reskit.parameters.parameters import OffshoreParameters

# GPS Coordiantes for location in North Sea
lon = 6.2160
lat = 53.7170


def test_waterDepthFromLocation():
    waterdepth_exact = 22
    depth = waterDepthFromLocation(
        lat, lon, waterDepthFolderPath=DEFAULT_PATHS.get("waterdepthFile")
    )
    assert np.isclose(
        depth, waterdepth_exact
    ), "the waterdepthfile is not working correct"


def test_calculateOffshoreCapex():
    comparedCAPEX = 2943.80819956
    calculatedCAPEX = calculateOffshoreCapex(
        inputCapex=3000,
        capacity=14000,
        hubHeight=150,
        waterDepth=25,
        coastDistance=25,
        rotorDiam=230,
        techYear=2050,
        shareTurb=0.449,
        shareFound=0.204,
        shareCable=0.181,
        shareOverhead=0.166,
        maxMonopileDepth=25,
        maxJacketDepth=55,
        litValueAvgDepth=17,
        litValueAvgDistCoast=27,
        baseCap=13000,
        baseHubHeight=150,
        baseRotorDiam=250,
        defaultOffshoreParamsFp=None,
    )

    assert np.isclose(calculatedCAPEX, comparedCAPEX)


def test_getRatedCostFromWaterDepth():

    test_value = getRatedCostFromWaterDepth(17, 25, 55)
    assert np.isclose(test_value, 431693), "equation is changed"
    assert getRatedCostFromWaterDepth(17) == getRatedCostFromWaterDepth(
        -17
    ), "negative avlues are handled incorrect"


def test_getCableCost():
    short = getCableCost(10, 14000, variableCostFactor=1.35, fixedCost=0)
    long = getCableCost(50, 14000, variableCostFactor=1.35, fixedCost=0)
    small = getCableCost(10, 10000, variableCostFactor=1.35, fixedCost=0)
    large = getCableCost(10, 14000, variableCostFactor=1.35, fixedCost=0)

    assert short < long, "equaiton is wrong in cable cost calculations"
    assert small < large, "equaiton is wrong in cable cost calculations"


def test_getFoundationCost():

    c1 = getFoundationCost(
        capacity=10000,
        applicationType="ac",  # AC substation offshore
        waterDepth=56,  # floating water depth
        foundationType=None,  # no given foundation
        maxJacketDepth=55,
        convention="RogeauEtAl2023",  # Rogeau et al
    )
    assert c1 == 461856.0

    c2 = getFoundationCost(
        capacity=10000,
        applicationType="dc",  # DC substation offshore
        waterDepth=56,  # floating water depth
        foundationType="jacket",  # jacket given but too deep -> warning, no error
        maxJacketDepth=55,
        convention="RogeauEtAl2023",  # Rogeau et al
    )
    assert c2 == 679784.0

    c3 = getFoundationCost(
        capacity=10000,
        applicationType="electrolysis",  # central offshore electrolysis
        waterDepth=55,  # jacket water depth
        foundationType="floating",  # jacket would have been possibel, too, but floating allowed
        maxJacketDepth=55,
        convention="RogeauEtAl2023",  # Rogeau et al
    )
    assert c3 == 1553080.0

    # TEST MUST-FAIL CASES
    with pytest.raises(Exception):
        getFoundationCost(
            capacity=10000,
            applicationType="does_not_exist",  # must fail
            waterDepth=55,
            foundationType=None,
            maxJacketDepth=55,
            convention="RogeauEtAl2023",
        )

    with pytest.raises(Exception):
        getFoundationCost(
            capacity=10000,
            applicationType="AC",
            waterDepth=50,
            foundationType="does_not_exist",  # must fail
            maxJacketDepth=55,
            convention="RogeauEtAl2023",
        )

    with pytest.raises(Exception):
        getFoundationCost(
            capacity=10000,
            applicationType="AC",
            waterDepth=-1,  # must fail
            foundationType=None,
            maxJacketDepth=55,
            convention="RogeauEtAl2023",
        )


def test_getConverterStationCost():
    # test onshore AC substation
    c1 = getConverterStationCost(
        capacity=10000,
        waterDepth=None,  # onshore
        voltageType="ac",
        maxJacketDepth=55,
        convention="RogeauEtAl2023",
    )

    assert c1 == 231875.0

    # test offshore DC substation
    c2 = getConverterStationCost(
        capacity=10000,
        waterDepth=55,  # jacket depth
        voltageType="dc",
        maxJacketDepth=55,
        convention="RogeauEtAl2023",
    )

    assert c2 == 1713505.0

    # TEST MUST-FAIL CASES

    with pytest.raises(Exception):
        getConverterStationCost(
            capacity=10000,
            waterDepth=55,  # jacket depth
            voltageType="does_not_exist",  # must fail
            maxJacketDepth=55,
            convention="RogeauEtAl2023",
        )
