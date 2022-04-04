from input import *
import load_files
import numpy as np
import struct
import pdb


def load_opacity(temperature, pressure, molecule):
    res = 2     # resolution for the opacities
    step_size = int(res/0.01)

    wavenumber_min = int(1e4/wavelength_bins[-1])
    wavenumber_max = int(1e4/wavelength_bins[0])

    index_min = int((wavenumber_min)/res)
    if res == 2:
        index_max = int((wavenumber_max)/res) - 1
    else:
        index_max = int((wavenumber_max)/res)

    temp_str = str(temperature).zfill(5)     # temperature as in opacity filename
    pressure_load = int(np.log10(pressure) * 100)

    if pressure_load < 0:
        pressure_str = 'n' + str(abs(pressure_load)).rjust(3, '0')  # pressure as in opacity filename
    else:
        pressure_str = 'p' + str(abs(pressure_load)).rjust(3, '0')

    wavenumber_dict = {'1H2-16O__POKAZATEL_e2b': '42000', '12C-1H4__YT10to10_e2b': '13000', '12C-16O__Li2015_e2b': '22000'}

    filename = molecule + '/Out_00000_' + wavenumber_dict[molecule] + '_' + temp_str + '_' + pressure_str + '.bin'

    data = []
    with open(opacity_path + filename, "rb") as f:
        byte = f.read(4)
        while byte:
            data.extend(struct.unpack("f", byte))
            byte = f.read(4)

    x_full = np.r_[0:42000:0.01]

    x_full = x_full[index_min * step_size:index_max * step_size:step_size]
    data = np.array(data[index_min * step_size:index_max * step_size:step_size])

    if len(data) < len(x_full):
        print('padding opacities...')
        data = np.pad(data, (0, len(x_full)-len(data)), 'constant')

    return data, x_full



def load_cia(x_full):

    sigma_h2h2 = load_files.load_sigma('H2', 'H2', x_full)
    sigma_h2he = load_files.load_sigma('H2', 'He', x_full)

    sigma_cia = sigma_h2h2 + 0.1*sigma_h2he # Because we assume xhe = 0.1*xh2

    return sigma_cia


def tau(p0):
    # Compute tau for all pressures

    pressure_array = 10 ** np.array([-8, -7.66, -7.33, -7, -6.66, -6.33, -6, -5.66, -5.33, -5, -4.66, -4.33, -4, -3.66,
                                     -3.33, -3, -2.66, -2.33, -2, -1.66, -1.33, -1, -0.66, -0.33, 0, 0.33, 0.66, 1.0])

    pressure_array_pmin = pressure_array[np.where(pressure_array == pmin)[0][0]:]  # remove everything below pmin

    pressure_array_pmin = pressure_array_pmin * 1e6  # convert to cgs
    p0 = p0 * 1e6

    wavenumber_min = int(1e4/wavelength_bins[-1])
    wavenumber_max = int(1e4/wavelength_bins[0])

    opacity_line_length = int((wavenumber_max - wavenumber_min) / res)
    if planet_name == 'GJ1214b':
        opacity_line_length = int((wavenumber_max - wavenumber_min) / res) - 1
    if planet_name == 'WASP-19b':
        opacity_line_length = int((wavenumber_max - wavenumber_min) / res) - 1
    if (opacity_line_length % 2) == 0:
        opacity_line_length = int((wavenumber_max - wavenumber_min) / res)
    else:
        opacity_line_length = int((wavenumber_max - wavenumber_min) / res) - 1

    integral_dict = {}


    for molecule in molecules:

        integral_grid_molecule = np.zeros((len(temperature_array), len(pressure_array_pmin), opacity_line_length))

        # we will integrate over pressure, for each temperature, for all wavelengths

        # Load integrands for all pressures
        for i, t in enumerate(temperature_array):

            integrand_grid_molecule = np.zeros((len(pressure_array_pmin), opacity_line_length))  # This will be the integrand for water

            for j, p in enumerate(pressure_array_pmin):
                p = p / 1e6  # Need bars to read in opacity file
                opacity, x_full = load_opacity(t, p, molecule)  # load water opacity for this temperature and pressure
                # print(x_full)
                integrand_grid_molecule[j] = opacity * molecular_mass_dict[molecule] / np.sqrt(np.log(p0 / p))  # compute sigma/sqrt(ln(P0/P))

            # Integrate for each pressure p, from pmin to p
            for j, p in enumerate(pressure_array_pmin):
                pressure_sliced = pressure_array_pmin[:j + 1]  # pass in pressure values and integrand values for all pressures
                integrand_grid_sliced = integrand_grid_molecule[:j + 1]  # below p, above pmin

                integral_value = np.trapz(integrand_grid_sliced, pressure_sliced,
                                          axis=0)  # calculate integral using trapezoid approximation

                integral_grid_molecule[i, j] = integral_value


        integral_dict[molecule] = integral_grid_molecule



    # now we do the same for CIA

    integral_cia_grid = np.zeros((len(temperature_array_cia), len(pressure_array_pmin), opacity_line_length))

    sigma_cia = load_cia(x_full)

    for i, t in enumerate(temperature_array_cia):

        integrand_cia_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))

        for j, p in enumerate(pressure_array_pmin):

            ntot = p/kboltz/t   # total number density in cgs
            integrand_cia_grid[j] = sigma_cia[i]*ntot / np.sqrt(np.log(p0 / p))


        for j, p in enumerate(pressure_array_pmin):
            pressure_sliced = pressure_array_pmin[:j + 1]  # pass in pressure values and integrand values for all pressures
            integrand_cia_grid_sliced = integrand_cia_grid[:j + 1]  # below p, above pmin

            integral_value = np.trapz(integrand_cia_grid_sliced, pressure_sliced,
                                      axis=0)  # calculate integral using trapezoid approximation

            integral_cia_grid[i, j] = integral_value


    # And Rayleigh scattering, but this one is independent of temperature

    sigma_rayleigh = 8.4909e-45 * (x_full ** 4)

    integral_rayleigh_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))
    integrand_rayleigh_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))

    for j, p in enumerate(pressure_array_pmin):

        integrand_rayleigh_grid[j] = sigma_rayleigh / np.sqrt(np.log(p0 / p))


    for j, p in enumerate(pressure_array_pmin):
        pressure_sliced = pressure_array_pmin[:j + 1]  # pass in pressure values and integrand values for all pressures
        integrand_rayleigh_grid_sliced = integrand_rayleigh_grid[:j + 1]  # below p, above pmin

        integral_value = np.trapz(integrand_rayleigh_grid_sliced, pressure_sliced,
                                  axis=0)  # calculate integral using trapezoid approximation

        integral_rayleigh_grid[j] = integral_value


    # Here is a dictionary of integral grids, each with different associated temperature arrays !!!

    integral_dict['cia'] = integral_cia_grid
    integral_dict['rayleigh'] = integral_rayleigh_grid

    return pressure_array_pmin, integral_dict, x_full

