import pandas as pd 
#dataframe1 = pd.read_excel('Umist5.xlsx')
# The purpose of this code is to read an excel spreadsheet of the umist database 
# and convert it into Julia catalyst code, the output of this python script
# is a giant string array where each element of that array is a line of julia code. 


all = pd.read_excel(r'/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/UMIST/UmistALL.xlsx')
rows = all.shape[0]
cols = all.shape[1]
col_labels = ['Ref. No.', 'Type', 'R1', 'R2', 'P1', 'P2', 'P3', 'P4', 'NE', 'Alpha', 'Gamma', 'Tl', 'Tu', 'ST', 'ACC', 'unknown', 'REF']

#print(all.at[0, 'P1'])

final_print = []
### 1: Packages ###
packages_print = ['### Packages ###',
                  'using Catalyst',
                  'using DifferentialEquations',
                  'using Plots',
                  'using Sundials, LSODA',
                  'using OrdinaryDiffEq',
                  'using Symbolics',
                  'using DiffEqDevTools',
                  'using ODEInterface, ODEInterfaceDiffEq',
                  'using ModelingToolkit', 
                  '', 
                  '']



### 2: Timespan ###
tspan_print = ['### Timespan ###', 
               'seconds_per_year = 3600 * 24 * 365',
               'tspan = (0.0, 30000 * seconds_per_year) # ~30 thousand yrs', 
               '', 
               '']

### 3: Initial Conditions ###
in_cond_print = ['### Initial Conditions ###', 'u0 = [']
ordered_spec_list = []
spec = ''
u0_line = ''

# this one's special because Julia will order the initial conditions sepcies
# on a first come first serve basis, so python needs to scan the reactions
# and find which one comes first, second, third, and so on.

# This will loop through all the species cells and create an ordered string array
# which will hopefully be the correct order that Julia will have
for i in range(rows):
    for j in range(6):
        col_name = col_labels[j+2]
        spec = all.at[i,col_name]
        if pd.isna(spec):
            break
        elif spec == 'PHOTON' or spec == 'CRP' or spec == 'CRPHOT':
            ordered_spec_list = ordered_spec_list
        elif spec not in ordered_spec_list:
            ordered_spec_list.append(spec)

# This creates the initial conditon julia code 
for s in range(len(ordered_spec_list)):
    if s == len(ordered_spec_list)-1:
        u0_line = '      0.0]   # ' + str(s+1) + ': ' + ordered_spec_list[s]
    else:
        u0_line = '      0.0,   # ' + str(s+1) + ': ' + ordered_spec_list[s]
    in_cond_print.append(u0_line)
in_cond_print.extend(['', ''])



### 4: Parameters ###
T = 10
Av = 2
n_H = 611
cr_ion_rate = 6e-18 # The cosmic-ray ionisation rates listed here are normalised to a total rate for electron production from cosmic-ray ionisation (primarily from H2 and He in dark clouds) of ζ0 = 1.36 × 10−17s−1. Rates for both direct cosmic-ray ionisation and cosmic-ray-induced photoreactions can be scaled to other choices of the ionisation rate, ζ,by multiplying the appropriate rate coefficients by ζ/ζ0.
omega = .5 # omega is the dust-grain albedo in the far ultraviolet (typically 0.4–0.6 at 150 nm).
# I don't know much about omega

react_type = ''
rate = ''
params_line = ''

params_print = ['### Parameters ###', 
                'T = 10',
                'omega = 0.5', 
                'Av = 2'
                'params = Dict(', 
                '         :T => T,', 
                '         :Av => Av,', 
                '         :n_H => 611,', 
                '         :cr_ion_rate => 6e-18,', 
                '         :omega => omega,']

# This will find the rate of the ith reaction, i.e., it will plug in the alpha, beta, gamma, values into the correct formula for the appropriate rezction type
for i in range(rows):
    react_type = all.at[i, 'Type']
    alpha = all.at[i, 'Alpha']
    beta = all.at[i, 'Beta']
    gamma = all.at[i, 'Gamma']

    # For Cosmic Ray Ionization (Type = CP) rate is given by eq. (2)
    if react_type == 'CP':
        rate = str(alpha)  # s^-1

    # Cosmic Ray Induced Photoionizations (Type = CR) are given by eq. (3)
    elif react_type == 'CR':
        rate = str(alpha) + ' * (T/300)^(' + str(beta) + ') * (' + str(gamma) + '/(1-omega))'  # s^-1  

    # Insterstellar Photoreactions (Type = PH) are given by eq. (4)
    elif react_type == 'PH':
        rate = str(alpha) + ' * exp(-(' + str(gamma) + '*Av))'  # s^-1

    # Two body reactions are given by eq. (1)
    else:
        rate = str(alpha) + ' * (T/300)^(' + str(beta) + ') * exp(-(' + str(gamma) + '/T))'  # cm^3 s^-1

    # Take out the comma on the last line
    if i == rows-1:
        params_line = '         :k' + str(i+1) + ' => ' + rate + ' # Reaction rate number ' + str(i+1) + ' with type ' + react_type
    else:
        params_line = '         :k' + str(i+1) + ' => ' + rate + ', ' + ' # Reaction rate number ' + str(i+1) + ' with type ' + react_type
    
    params_print.append(params_line)

params_print.append('         )')
params_print.append('')
params_print.append('')


### 4.5 Tempurature range flags ###
temp_cond_print = ['### Tempurature range flags ###']
T_low_bounds = 'T_low_bounds = ['
T_upp_bounds = 'T_upp_bounds = ['
for i in range(rows):
    T_low_bounds = T_low_bounds + str(all.at[i, 'Tl']) + ' '
    T_upp_bounds = T_upp_bounds + str(all.at[i, 'Tu']) + ' '

T_low_bounds = T_low_bounds + ']'
T_upp_bounds = T_upp_bounds + ']'

temp_cond_print.append(T_low_bounds)
temp_cond_print.append(T_upp_bounds)
temp_cond_print.append('for i = 1:length(T_low_bounds)')
temp_cond_print.append('    if T < T_low_bounds[i]') 
temp_cond_print.append('        print("\\nTempurature below the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) ')
temp_cond_print.append('    elseif T > T_upp_bounds[i]') 
temp_cond_print.append('        print("\\nTempurature higher than the tempurature range [", T_low_bounds[i], ", ", T_upp_bounds[i], "] for reaction number ", i ) ')
temp_cond_print.append('    end')
temp_cond_print.append('end')
temp_cond_print.append('print("\\n")')
temp_cond_print.append('print("\\n")')
temp_cond_print.append('')
temp_cond_print.append('')




# 4586 reactions have tempurature ranges [10, 41000]
# 1060 reactions have tempurature ranges [10, 300]
# 2 have range [200, 700]
# 21 have range [10, 2500]
# 1 has range [2720, 5190]
# 1 has range [300, 3160]
# 1 has range [383, 41000]
# 1 has range [1250, 1600]
# 35 have range [300, 500]

# total ~5669 (out of 6173)


# for parameters, we can make every single reaction rate a parameter, (is this redundant storage?)
# conditional parameters will need to have custom if else statements, 
# this means we are already looking at a good 32,000 lines of julia code.
# Is there a better way to do this?

# Nevermind, consider this:
# UMIST consists of 6173 gas phase reactions, 6153 have only 1 tempurature range, 
# leaving only 20 reactions with 2 tempurature ranges 
# I am going to customize the 20 reactions with two temperature ranges





### 5: Network Admin things ###
net_admin_print = ['### Network Admin things ###', '@variables t ']
spec_macro_line = '@species '
params_macro_line = '@parameters T Av n_H cr_ion_rate omega '

for i in range(len(ordered_spec_list)):
    spec_macro_line = spec_macro_line + str(ordered_spec_list[i]) + '(t) '
net_admin_print.append(spec_macro_line)

for i in range(rows):
    params_macro_line = params_macro_line + 'k' + str(i+1) + ' '
net_admin_print.append(params_macro_line)

net_admin_print.append('')
net_admin_print.append('')


### 6: Reaction Equations ###
react_eqs_print = ['### Reaction Equations ###', 'reaction_equations = [']
ith_reaction = ''
spec = ''
R_labels = ['R1', 'R2']
P_labels = ['P1', 'P2', 'P3', 'P4']

for i in range(rows):
    ith_reaction = '    (@reaction n_H * k' + str(i+1) + ', '
    # Reactants
    for j in range(2):
        spec = all.at[i, R_labels[j]]
        if pd.isna(spec) or spec == 'PHOTON' or spec == 'CRP' or spec == 'CRPHOT':
            break
        elif j == 1:
            ith_reaction = ith_reaction + ' + ' + spec
        else:
            ith_reaction = ith_reaction + spec

    ith_reaction = ith_reaction + ' --> '

    # Products
    for j in range(4):
        spec = all.at[i, P_labels[j]]
        if pd.isna(spec) or spec == 'PHOTON' or spec == 'CRP' or spec == 'CRPHOT':
            ith_reaction = ith_reaction[:-3]
            break
        elif j == 3:
            ith_reaction = ith_reaction + spec
        else:
            ith_reaction = ith_reaction + spec + ' + '


    if i == rows-1:
        ith_reaction = ith_reaction + ' )] '
    else:
        ith_reaction = ith_reaction + ' ), '

    react_eqs_print.append(ith_reaction)

react_eqs_print.append('')
react_eqs_print.append('')

### 7: Turn the Network into a system of ODEs ###
convert_print = ['### Turn the Network into a system of ODEs ###']
convert_print.append('@named system = ReactionSystem(reaction_equations, t)')
convert_print.append('sys = convert(ODESystem, complete(system))')
convert_print.append('ssys = structural_simplify(sys)')
convert_print.append('prob = ODEProblem(ssys, u0, tspan, params)')
convert_print.append('sol = solve(prob, Rodas4())')
convert_print.append('')
convert_print.append('')

### 7.5: Timing ###
time_print = ['### Timing ###']
time_print.append('print("Time to convert:")')
time_print.append('@time convert(ODESystem, complete(system))')
time_print.append('print("Time to simplify:")')
time_print.append('@time structural_simplify(sys)')
time_print.append('print("Time to create the simplified problem:")')
time_print.append('@time ODEProblem(ssys, u0, tspan, params)')
time_print.append('print("Time to solve the simplified 1000 reaction system with Rodas4(): ")')
time_print.append('@time solve(prob, Rodas4());')
time_print.append('')
time_print.append('')

### 8: Plot ###
plot_print = ['### Plotting ###']
plot_line1 = ''
plot_line2 = ''

for i in range(len(ordered_spec_list)):
    plot_line1 = '# Species number ' + str(i+1) + ': ' + ordered_spec_list[i]
    plot_line2 = 'plot(sol, idxs = (0,' + str(i+1) + '), lw = 3, lc = "blue", title = "Umist: Abundance of ' + ordered_spec_list[i] + '")'
    plot_print.append(plot_line1)
    plot_print.append(plot_line2)
    plot_print.append('')

plot_print.append('')
plot_print.append('')


### 9: Print it all ###
final_print.extend(packages_print)
final_print.extend(tspan_print)
final_print.extend(in_cond_print)
final_print.extend(params_print)
final_print.extend(temp_cond_print)
final_print.extend(net_admin_print)
final_print.extend(react_eqs_print)
final_print.extend(convert_print)
final_print.extend(time_print)
final_print.extend(plot_print)


len_fin = len(final_print)

for line in final_print:
    print(line)

print(all)



#file_path = "/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/UMIST/Umist.txt"  # Replace with your desired file path
#with open(file_path, "w") as file:
#    for i in range(len_fin):
#       file.writelines(final_print[i] + '\n')

# only reactions with two reacants shouls be multiplied by n_H
# Two problems: replace e- with e in excel, 
# do the same for all the pluses and minuses for all the species
# Also, this program currently neglects the 20 reactions with two tempurature ranges instead of one

