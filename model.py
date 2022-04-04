from input import *
from load_files import *
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pdb


class Model:

    def __init__(self, len_x, x_full, bin_indices, parameter_dict, integral_dict):

        self.len_x = len_x
        self.x_full = x_full
        self.bin_indices = bin_indices
        self.parameter_dict = parameter_dict
        # self.opacity_grid = opacity_grid
        self.integral_dict = integral_dict


    def molecular_opacity(self, my_temp, opacity_table):

        fn = RegularGridInterpolator((self.x_full, temperature_array), opacity_table)
        pt = (self.x_full, my_temp)
        y = fn(pt)
        return y


    def cia_cross_section(self, my_temp, cia_table):

        fn = RegularGridInterpolator((temperature_array_cia, self.x_full), cia_table)
        pt = (my_temp, self.x_full)
        y = fn(pt)
        return y


    def integral_grid_interpolate(self, my_temp, temp_arr, pressure_arr, integral_grid):

        # pdb.set_trace()

        fn = RegularGridInterpolator((temp_arr, pressure_arr, self.x_full), integral_grid) # Interpolate on 3D grid
        pts = np.meshgrid(np.array([my_temp]), pressure_arr, self.x_full)
        flat_pts = np.array([m.flatten() for m in pts])     # Some complicated fix for 3D
        out_array = fn(flat_pts.T)
        # y = out_array.reshape(*pts[0].shape)
        y = out_array.reshape(len(pressure_arr), len(self.x_full))


        return y


    def transit_depth(self):
        ## calculates transit depth ##

        my_temp = self.parameter_dict["T"]
        R0 = self.parameter_dict["R0"]*rjup
        P0 = (10**self.parameter_dict["log_P0"])
        Rstar = self.parameter_dict["Rstar"]*rsun
        G = self.parameter_dict["G"]

        epsilon = 0.000001
        p0 = (P0 + epsilon)*1e6   # add some tiny value to p0 to avoid infinities in integration

        pressure_array = 10 ** np.array(
            [-8, -7.66, -7.33, -7, -6.66, -6.33, -6, -5.66, -5.33, -5, -4.66, -4.33, -4, -3.66,
             -3.33, -3, -2.66, -2.33, -2, -1.66, -1.33, -1, -0.66, -0.33, 0, 0.33, 0.66, 1.0])

        pressure_values = pressure_array[np.where(pressure_array == pmin)[0][0]:]  # remove everything below pmin
        pressure_values = pressure_values*1e6   # convert to cgs

        mass_fraction = []
        molecular_mass = []
        integral_molecules = []

        for molecule in molecules:
            abundance_name = molecular_abundance_dict[molecule]
            mass_fraction.append([10**self.parameter_dict[abundance_name]])
            molecular_mass.append([molecular_mass_dict[molecule]])
            abundance_name = molecular_abundance_dict[molecule]
            integral_mol = self.integral_grid_interpolate(my_temp, temperature_array, pressure_values, self.integral_dict[molecule]) * (10 ** self.parameter_dict[abundance_name])
            integral_molecules.append(integral_mol)

        integral_molecules_sum = np.sum(np.array(integral_molecules), axis=0)

        mass_fraction = np.array(mass_fraction)
        molecular_mass = np.array(molecular_mass)

        xh2 = (1 - np.sum(mass_fraction))/1.1   # calculate abundance of H2
        if 'm' not in globals():                # set mean molecular weight if not given in input
            m = 2.4*xh2*amu + np.sum(mass_fraction*molecular_mass)

            
        scale_height = kboltz*my_temp/m/G

        factor = np.sqrt(2 * scale_height * R0) / (kboltz * my_temp)


        integral_grid_cia = self.integral_grid_interpolate(my_temp, temperature_array_cia, pressure_values, self.integral_dict['cia'])
        integral_grid_rayleigh = self.integral_dict['rayleigh']

        if self.parameter_dict["log_P_cloudtop"] == "Off":  # Cloud free

            cloud_tau = np.zeros((len(pressure_values),1))

        elif self.parameter_dict["Q0"] == "Off":        # Grey cloud

            cloudtop_pressure = (10**self.parameter_dict["log_P_cloudtop"])*1e6
            cloud_tau = np.zeros((len(pressure_values),1))
            cloud_tau[pressure_values > cloudtop_pressure,0] = np.inf

        else:       # Non-grey cloud

            cloudtop_pressure = (10**self.parameter_dict["log_P_cloudtop"])*1e6
            bc = 10**self.parameter_dict["log_cloud_depth"]
            cloudbottom_pressure = bc*cloudtop_pressure

            Q0 = self.parameter_dict["Q0"]
            a = self.parameter_dict["a"]
            rc = 10**self.parameter_dict["log_r_c"]
            tau_ref = 10**self.parameter_dict["log_tau_ref"]
            x_lambda = 2 * np.pi * rc * self.x_full
            x_lambda_ref = 2 * np.pi * rc * 1e4
            tau_nongrey = tau_ref * (Q0 * x_lambda_ref ** (-a) + x_lambda_ref ** 0.2) / (Q0 * x_lambda ** (-a) + x_lambda ** 0.2)

            cloud_tau = np.zeros((len(pressure_values),len(self.x_full)))
            cloud_tau[(pressure_values > cloudtop_pressure) & (pressure_values < cloudbottom_pressure),:] = tau_nongrey


        tau_values = factor * (integral_molecules_sum + integral_grid_cia*xh2 + integral_grid_rayleigh*xh2) + cloud_tau

        h_values = np.zeros(len(tau_values[0]))

        for i in range(len(tau_values[0])):        # This should be 1458

            new_integrand = np.zeros(len(pressure_values))
 
            for j in range(len(pressure_values)):       # This should be 22

                new_integrand[j] = (1 - np.exp(-tau_values[j,i]))/pressure_values[j]*(R0 + scale_height*np.log(p0/pressure_values[j]))
                # print(new_integrand)


            integral_value = np.trapz(new_integrand, pressure_values)     # This should be a scalar
            # print(integral_value)
           
            h_values[i] = (scale_height/R0)*integral_value


        r = R0 + h_values
        
        result = 100.0 * (r / Rstar) ** 2  # return percentage transit depth

        # with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_3]/save_transit_results_'+model_name+'.txt', 'w') as save_transit_results:
        #     save_transit_results.write(str(np.mean(result)))

        return result
        



    def binned_model(self):
    ## calculates average transit depth in given bins ##

        if self.parameter_dict['line'] == 'Off':
            y_full = self.transit_depth()
            y_mean = np.zeros(self.len_x)
            for i in range(self.len_x):
                j = int(self.bin_indices[i])
                k = int(self.bin_indices[i + 1])
                y_mean[i] = np.mean(y_full[k:j])  # bin transit depth      
        else:
            y_mean = np.full(self.len_x, self.parameter_dict['line'])
            
        return y_mean
