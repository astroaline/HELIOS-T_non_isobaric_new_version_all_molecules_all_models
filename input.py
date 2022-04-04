import numpy as np
import os

## Constants ##

kboltz = 1.38064852e-16    # Boltzmann's constant
amu = 1.660539040e-24      # atomic mass unit
gamma = 0.57721
rjup = 7.1492e9            # equatorial radius of Jupiter
rsun = 6.9566e10           # solar radius
rearth = 6.378e8            # earth radius
pressure_probed = 1e-2      # probed pressure in bars
# pressure_cia = 1e-2         # pressure for cia in bars
# m = 2.4*amu                 # assummed hydrogen-dominated atmosphere
m_water = 18.0*amu          # mean molecular mass of any molecules you want to consider
m_cyanide = 27.0*amu
m_ammonia = 17.0*amu
m_methane = 16.0*amu
m_carbon_monoxide = 28.0*amu


## Planet Data ##

planet_name = 'TRAPPIST-1d'

g = 474
g_uperr = 51
g_loerr = 47
g_uncertainty = (g_uperr+g_loerr)/2
rstar = 0.121
rstar_uncertainty = 0.003
r0 = 0.05402   # Rp = [0.068,0.072]
r0_uncertainty = 0.02   # R0 = [0.03402,0.07402], R0 = [0.6668,0.7868]

wavelength_bins = np.array([1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65]) # must be 1 longer than transit_depth and transit_depth_error
transit_depth = np.array([0.3723,0.3992,0.3987,0.4163,0.4234,0.38420000000000004,0.4106,0.38680000000000003,0.40650000000000003,0.3922])
transit_depth_error = np.array([0.027200000000000002,0.024300000000000002,0.0229,0.028,0.0176,0.0205,0.0164,0.018600000000000002,0.0201,0.0189])

pmin = 1e-6



## Retrieval info ##

model_name = 'greycloud'

molecules = ['1H2-16O__POKAZATEL_e2b']  # list of molecules (determines which opacity tables are loaded)
parameters = ["T", "log_xh2o", "R0", "log_P_cloudtop", "Rstar", "G"]   # parameters you wish to retrieve (MUST MATCH MOLECULES)
res = 2         # resolution used for opacities
live = 1000     # live points used in nested sampling
wavenumber=True     # True if opacity given in terms of wavenumber, False if wavelength

priors = {"T": [2700, 200], "log_xh2o": [13,-13], "log_xch4": [13,-13], "log_xco": [13,-13],
          "log_P0": [4,-1], "R0": [2*r0_uncertainty, r0-r0_uncertainty], "log_tau_ref": [7,-5], "Q0": [99,1], "a": [10,3],
          "log_r_c": [6,-9], "log_p_cia": [3,-3], "log_P_cloudtop": [5,-4], "log_cloud_depth": [2,0],
          "Rstar": [2*rstar_uncertainty,rstar-rstar_uncertainty],
          "G": [2*g_uncertainty,g-g_uncertainty], "line": [5,0]} # priors for all possible parameters



## info for all possible parameters ##
molecular_name_dict = {'1H2-16O__POKAZATEL_e2b': 'water', '12C-1H4__YT10to10_e2b': 'methane', '12C-16O__Li2015_e2b': 'carbon_monoxide'}  # dictionary list of all possible molecules and corresponding names
molecular_abundance_dict = {'1H2-16O__POKAZATEL_e2b': 'log_xh2o', '12C-1H4__YT10to10_e2b': 'log_xch4', '12C-16O__Li2015_e2b': 'log_xco'}  # dictionary list of all possible molecules and corresponding abundance names

parameter_dict = {"T": 1000, "log_xh2o": "Off", "log_xch4": "Off", "log_xco": "Off", "R0": r0,
                  "Rstar": rstar, "log_P0": 1, "log_tau_ref": "Off", "Q0": "Off", "a": "Off", "log_r_c": "Off", "log_p_cia": -2,
                  "log_P_cloudtop": "Off", "log_cloud_depth": "Off", "G": g, "line": "Off"}    # default parameter values used if not retrieved

molecular_mass_dict = {'1H2-16O__POKAZATEL_e2b': m_water, '12C-1H4__YT10to10_e2b': m_methane, '12C-16O__Li2015_e2b': m_carbon_monoxide}   # dictionary of molecules and their mean molecular masses
temperature_array = np.r_[50:700:50, 700:1500:100, 1500:3100:200]
#temperature_array = np.array([1500, 1700, 1900, 2100])
#temp_dict = {'01': temperature_array[9:], '12C-1H4__YT10to10_e2b': temperature_array[9:], '12C-16O__HITEMP2010_e2b': temperature_array}   # temperature values for corresponding opacity tables
temperature_array_cia = np.r_[200:3025:25]          # temperature array for CIA table
opacity_path = os.environ['HOME'] + "/Desktop/PhD/OPACITIES/"  # path to opacity binary files
cia_path = os.environ['HOME'] + "/Desktop/PhD/HITRAN/"      # path to CIA files
