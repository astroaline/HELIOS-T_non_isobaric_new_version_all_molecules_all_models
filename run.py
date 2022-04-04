import os

planet_name = input('Enter the name of the planet: ')


with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_all_models]/planets/'+planet_name+'/input_non_isobaric_cloudfree.py', 'r') as f:
    data = f.read()

with open('input.py', 'w') as f:
    f.write(data)

os.system('python helios-t.py')



with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_all_models]/planets/'+planet_name+'/input_non_isobaric_flat_line.py', 'r') as f:
    data = f.read()

with open('input.py', 'w') as f:
    f.write(data)

os.system('python helios-t.py')



with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_all_models]/planets/'+planet_name+'/input_non_isobaric_greycloud.py', 'r') as f:
    data = f.read()

with open('input.py', 'w') as f:
    f.write(data)

os.system('python helios-t.py')



with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_all_models]/planets/'+planet_name+'/input_non_isobaric_greycloud_CH4_CO.py', 'r') as f:
    data = f.read()

with open('input.py', 'w') as f:
    f.write(data)

os.system('python helios-t.py')



with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_all_models]/planets/'+planet_name+'/input_non_isobaric_non_greycloud.py', 'r') as f:
    data = f.read()

with open('input.py', 'w') as f:
    f.write(data)

os.system('python helios-t.py')



with open('/home/aline/Desktop/PhD/HELIOS-T-master_[non_isobaric_new_version_all_molecules_all_models]/planets/'+planet_name+'/input_non_isobaric_non_greycloud_CH4_CO.py', 'r') as f:
    data = f.read()

with open('input.py', 'w') as f:
    f.write(data)

os.system('python helios-t.py')

