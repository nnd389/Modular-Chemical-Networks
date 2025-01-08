import pandas as pd 
import numpy as np
all = pd.read_excel(r'/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/UMIST/Umist_and_Glover.xlsx')

col_labels = ['UI', 'URa', 'UT', 'URe', 'Tl', 'Tu', 'Column7', 'GI', 'GRe']
line = ''

UR1 = []
UR2 = []
UP1 = []
UP2 = []
UP3 = []
UP4 = []
for i in range(6173):
    UR1.append(all.at[i, 'UR1'])
    UR2.append(all.at[i, 'UR2'])
    UP1.append(all.at[i, 'UP1'])
    UP2.append(all.at[i, 'UP2'])
    UP3.append(all.at[i, 'UP3'])
    UP4.append(all.at[i, 'UP4'])

GR1 = []
GR2 = []
GR3 = []
GP1 = []
GP2 = []
GP3 = []
for i in range(218):
    GR1.append(all.at[i, 'GR1'])
    GR2.append(all.at[i, 'GR2'])
    GR3.append(all.at[i, 'GR3'])
    GP1.append(all.at[i, 'GP1'])
    GP2.append(all.at[i, 'GP2'])
    GP3.append(all.at[i, 'GP3'])

gspec = ''
label = ''
index_list1 = []


count = 0

glover_reaction = []
umist_reaction = []
final_print = ['### Parameters ###', 'T = 10','omega = 0.5','Av = 2','params = Dict(','         :T => T,','         :Av => Av,','         :n_H => 611,','         :cr_ion_rate => 6e-18,','         :omega => omega,',]
T_low_print = ['T_low_bounds = [']
T_upp_print = ['T_upp_bounds = [']


match_index_list = []
exit_code = 0


#change to range 218
for i in range(218):
    match_index_list = []
    flag = 0
    count = 0
    count_fin = 0
    exit_code = 0
    line = ''
    # Check is there are three reactants, if there is, then the reaction is definitely not in the UMIST database
    if pd.isna(GR3[i]) == False:
        line = '         # :k' + str(i+1) + ' => ' + 'Three reactants found, no match for Glover reaction ' + str(i+1) + str(all.at[i, 'GRe'])
        T_low_print.append('                -1   # No match for Glover reaction ' + str(i+1) + ' So no temp range')
        T_upp_print.append('                1e10 # No match for Glover reaction ' + str(i+1) + ' So no temp range')
        final_print.append(line)
    else:
        for row_UR1, element_UR1 in enumerate(UR1):
            
            # Look for two body forward reactant matches
            if GR1[i] == element_UR1:
                if GR2[i] == UR2[row_UR1]:
                    index_list1.append(row_UR1)
                    match_index_list.append(row_UR1)
                    count = count + 1
            # Look for two body backward reactant matches
            elif GR2[i] == element_UR1:
                if GR1[i] == UR2[row_UR1]:
                    index_list1.append(row_UR1)
                    match_index_list.append(row_UR1)
                    count = count + 1

            # Look for single body matches
            if (GR1[i] == element_UR1) and (i >=165):
                index_list1.append(row_UR1)
                match_index_list.append(row_UR1)
                count = count + 1

        
        # create glover individual product set
        glover_reaction = []
        glover_reaction.append(all.at[i, 'GP1'])
        glover_reaction.append(all.at[i, 'GP2'])
        glover_reaction.append(all.at[i, 'GP3'])
        # clean up nan values
        for j in range(len(glover_reaction)):
            if np.nan in glover_reaction:
                glover_reaction.remove(np.nan)

        # Check if there are no matches:
        if len(match_index_list) == 0:
            line = '         # :k' + str(i+1) + ' => ' + " Glover " + str(i+1) + str(all.at[i, 'GRe']) + " had " + str(count) + " product matches"
            T_low_print.append('                -1   # No match for Glover reaction ' + str(i+1) + ' So no temp range')
            T_upp_print.append('                1e10 # No match for Glover reaction ' + str(i+1) + ' So no temp range')
            final_print.append(line)

        # create umist individual product set
        for k in match_index_list:
            umist_reaction = []
            umist_reaction.append(all.at[k, 'UP1'])
            umist_reaction.append(all.at[k, 'UP2'])
            umist_reaction.append(all.at[k, 'UP3'])
            # clean up nan values and clean up PHOTON,CRP,CRPHOT values
            for j in range(len(umist_reaction)):
                if np.nan in umist_reaction:
                    umist_reaction.remove(np.nan)
                elif 'PHOTON' in umist_reaction:
                    umist_reaction.remove('PHOTON')
                elif 'CRP' in umist_reaction:
                    umist_reaction.remove('CRP')
                elif 'CRPHOT' in umist_reaction:
                    umist_reaction.remove('CRPHOT')
            # Check if the products match
            if len(glover_reaction) != len(umist_reaction):
                exit_code = 1
                
            elif set(glover_reaction) == set(umist_reaction):
                line = '         :k' + str(i+1) + ' => ' + str(all.at[k, 'URa']) + ",     # Found Match for Glover " + str(i+1) + str(all.at[i, 'GRe']) + " and Umist " + str(k+1) + str(all.at[k, 'URe']+ 'with temp range [' + str(all.at[k, 'Tl']) +  ', ' + str(all.at[k, 'Tu']) + ']')
                T_low_print.append('                ' + str(all.at[k, 'Tl']) + '   # Glover reaction ' + str(i+1) + ' matches with Umist ' + str(k+1) + ' with temp range [' + str(all.at[k, 'Tl']) +  ', ' + str(all.at[k, 'Tu']) + ']')
                T_upp_print.append('                ' + str(all.at[k, 'Tu']) + '   # Glover reaction ' + str(i+1) + ' matches with Umist ' + str(k+1) + ' with temp range [' + str(all.at[k, 'Tl']) +  ', ' + str(all.at[k, 'Tu']) + ']')
                final_print.append(line)
                #break
                
            
    #check reaction 110
    if exit_code == 1:
        line = '         # :k' + str(i+1) + ' => ' + 'NO match for Glover reaction ' + str(i+1) + str(all.at[i, 'GRe'])
        T_low_print.append('                -1   # No match for Glover reaction ' + str(i+1) + ' So no temp range')
        T_upp_print.append('                1e10 # No match for Glover reaction ' + str(i+1) + ' So no temp range')
        final_print.append(line)


    final_print.append('')
    T_low_print.append('')
    T_upp_print.append('')
    

final_print.append('         )')
T_low_print.append('                ]')
T_upp_print.append('                ]')
final_print.extend(T_low_print)
final_print.extend(T_upp_print)


for x in final_print:
    print(x)


len_fin = len(final_print)
file_path = "/Users/kneenaugh/Desktop/Git/AstroChemNetwork/scripts/UMIST/glover_params.txt"  # Replace with your desired file path
with open(file_path, "w") as file:
    for i in range(len_fin):
        file.writelines(final_print[i] + '\n')


#params_glover_print.append(line)
    


#for x in params_glover_print:
#    print(x)






















'''
g_param_print = ['T = 10', 'omega = 0.5', 'Av = 2', 'params = Dict(', '         :T => T,', '         :Av => Av,', '         :n_H => 611,', '         :cr_ion_rate => 6e-18,', '         :omega => omega,' ]
line = ''
params_macro_line = '@parameters T Av n_H cr_ion_rate omega '
for i in range(218):
    line = '         :k' + str(i+1) + ' => '
    params_macro_line = params_macro_line + 'k' + str(i+1) + ' '
    g_param_print.append(line)

for x in g_param_print:
    print(x)

print(params_macro_line)
'''