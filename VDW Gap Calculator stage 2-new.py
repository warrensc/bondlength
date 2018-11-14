import os
import glob
import math
from collections import OrderedDict

read_path = os.path.join(os.sep,  'Users', 'Admin', 'Desktop',  'jemsData', 'Materials genome',  'P1_cif_out')
write_path = os.path.join(os.sep,  'Users', 'Admin', 'Desktop', 'jemsData',  'Materials genome',  'BondLength-17')

cut_off_angle_is_a_constant = 1 # 0 is no, 1 is yes
cut_off_angle = 60
cut_off_angle_end = 75
angle_distance_lower_bound = 0.56 #specifies upper bound of cut_off_angle
angle_distance_upper_bound = 0.60 #specifies lower bound of cut_off_angle_end
examine_all_bonds = 1  # 0 is no, 1 is yes
target_atom_1 = 'O'
target_atom_2 = 'Ca'

def assure_path_exists(name_of_path):
    if not os.path.exists(name_of_path):
        os.makedirs(name_of_path)

assure_path_exists(write_path)

cif_counter = len(glob.glob1(read_path,"*.cif"))
counter = 0

z_count = 0
os.chdir(write_path)

while z_count < 80:
    z_name = 'vdw_' + str(z_count) + '.txt'
    out = open(z_name,  'w+')
    out.close()
    z_count += 1

element = ['', 'Intentionally empty', 'H',  'He',     'Li',    'Be',     'B',       'C',      'N',      'O',      'F',      'Ne',     'Na',   'Mg',     'Al',    'Si',     'P',      'S',      'Cl',    'Ar',        'K',      'Ca',     'Sc',     'Ti',    'V',       'Cr',    'Mn',     'Fe',    'Co',    'Ni',    'Cu',     'Zn',    'Ga',     'Ge',    'As',     'Se',    'Br',    'Kr',      'Rb',    'Sr',     'Y',       'Zr',    'Nb',    'Mo',    'Tc',      'Ru',    'Rh',   'Pd',      'Ag',     'Cd',    'In',     'Sn',    'Sb',    'Te',     'I',       'Xe',     'Cs',      'Ba',    'La',     'Ce',     'Pr',      'Nd',     'Pm',    'Sm',      'Eu',       'Gd',      'Tb',     'Dy',      'Ho',     'Er',      'Tm',      'Yb',      'Lu',     'Hf',     'Ta',     'W',     'Re',    'Os',     'Ir',     'Pt',    'Au',    'Hg',     'Tl',     'Pb',    'Bi',     'Po',   'At',    'Rn',    'Ac',    'Th',    'Pa',        'U',        'Np',     'Pu',       'Am',    'Cm',     'Bk',    'Cf',  'Es']
radius = [  '1', '1',                             '1.2', '1.43', '2.12', '1.98',  '1.91', '1.77', '1.66', '1.50', '1.46', '1.58', '2.4',  '2.51', '2.25', '2.19', '1.90', '1.89', '1.82',  '1.83',   '2.73', '2.62', '2.58',  '2.46', '2.42', '2.45', '2.45', '2.44', '2.40', '2.40', '2.38', '2.39', '2.32', '2.29', '2.05', '1.82', '1.86', '2.25',  '3.21', '2.84', '2.75', '2.52', '2.56', '2.45', '2.44', '2.44', '2.44', '2.15', '2.53',  '2.49', '2.43', '2.42', '2.47', '1.99', '2.04',  '2.06',  '3.48',  '3.03', '2.98', '2.88',  '2.92',  '2.95',  '2.92',  '2.90',  '2.87',  '2.83',   '2.79',  '2.87',  '2.81',  '2.83',   '2.79',  '2.80',  '2.74',  '2.63', '2.53', '2.57', '2.49', '2.48', '2.41', '2.29', '2.32', '2.45', '2.47',  '2.60', '2.54',   '2.4',   '2.2',   '2',   '2.8',  '2.93',  '2.88',   '2.71',   '2.82',  '2.81',  '2.83',  '3.05',  '3.4',  '3.05',  '2.7']

while counter < cif_counter:
    #1 Read in file.
    new_counter = str(counter).zfill(5)
    counter += 1
    os.chdir(read_path)
    f = open(new_counter + '.cif',  'r')
    f_lines = f.readlines()
    f.close()
    
    if examine_all_bonds == 1:
        atom = 1
    else:
        atom = 0

###########################################    
    for line in f_lines:
        if line.startswith(target_atom_1):
            for line in f_lines:
                if line.startswith(target_atom_2):
                    atom = 1
########################################## 
    
    if atom == 0:
        continue
    else:
        #1a  Create variable that has total number of atoms, "atom_count"
        atom_count = int(f_lines[0])
        if atom_count > 20:
            continue
        else:
            #1b  Create cell_positions[0],[1],[2] = a, b, c
            #1c  Create cell_angles[0],[1],[2]= alpha, beta, gamma
            cell_positions = f_lines[1].split(' ')
            cell_angles = f_lines[2].split(' ')

            s = 0
            while s <= 2:
                cell_positions[s] = float(cell_positions[s])
                cell_angles[s] = float(cell_angles[s])
                s+=1 

            #1d  Create atom list.  atom_id[0]...atom_id[atoms].  [0]=empty, starting really on [1]
            #1e  Create position list.
            #  atom_fx[0]....atom_fx[atom_count]  ...  with atom_fx[0] = 1.... thus starting with real data at x[1]
            #  atom_fy[0]....atom_fy[atom_count]   ... with atom_fy[0] = 1
            #  atom_fz[0]...atom_fz[atom_count]    ...  with atom_fz[0] = 1
            s = 3
            s_max = s + atom_count
            atom_id = [str('')]*2000
            atom_fx = [float(1)]*2000
            atom_fy = [float(1)]*2000
            atom_fz = [float(1)]*2000
            atom_id[0] = str('Intentionally empty')
            while s < s_max:
                parts = f_lines[s].split(' ')
                atom_id[s-2]=parts[0]
                atom_fx[s-2]=float(parts[1])
                atom_fy[s-2]=float(parts[2])
                atom_fz[s-2]=float(parts[3])
                s += 1
            #1f Create full atom list for ALL atoms.  all_atom_id = 0 to 27*atom_count - 1    -->   index [n] below
            all_atom_id = [str('')]*2000*27
            s=0
            s_max = 27
            while s<s_max:
                s2 = 0
                s2_max = atom_count
                while s2 < s2_max:
                    s3 = s2 + 27 * s
                    all_atom_id[s3] = atom_id[s2+1]
                    s2 += 1
                s += 1
            #2  Calculate unit cell in cartesian format  ( this is the tranlation matrix, multiplied by 1 ).  Outputs:  x[0], y[0], z[0].  One index for each atom.
            x = [float(0)]*2000
            y = [float(0)]*2000
            z = [float(0)]*2000
            
            cell_xx = cell_positions[0]
            cell_xy = 0
            cell_xz = 0 
            cell_yx = cell_positions[1] * math.cos(math.radians((cell_angles[2])))
            cell_yy = cell_positions[1] * math.sin(math.radians((cell_angles[2])))
            cell_yz = 0
            cell_zx = cell_positions[2] * math.cos(math.radians((cell_angles[1])))
            cell_zy = cell_positions[2] * (math.cos(math.radians((cell_angles[0])))  - math.cos(math.radians((cell_angles[1]))) * math.cos(math.radians((cell_angles[2])))  ) / (math.sin(math.radians(cell_angles[2])) ) 
            cell_zz = cell_positions[0] *  cell_positions[1]  *  cell_positions[2]* (1-(math.cos(math.radians((cell_angles[0]))))**2 - (math.cos(math.radians((cell_angles[1]))))**2-(math.cos(math.radians((cell_angles[2]))))**2 + 2* math.cos(math.radians((cell_angles[0]))) * math.cos(math.radians((cell_angles[1]))) * math.cos(math.radians((cell_angles[2]))))**(1/2)/(cell_positions[0] * cell_positions[1] * math.sin(math.radians((cell_angles[2])))  )
            
            x[0]    =    cell_positions[0] * atom_fx[0]  +  ( cell_positions[1] * math.cos(math.radians((cell_angles[2])))  ) * atom_fy[0]    +    ( cell_positions[2] * math.cos(math.radians((cell_angles[1]))) ) * atom_fz[0] 
            y[0]    =    ( cell_positions[1] * math.sin(math.radians((cell_angles[2])))) * atom_fy[0]  +  ( cell_positions[2] * (math.cos(math.radians((cell_angles[0])))  - math.cos(math.radians((cell_angles[1]))) * math.cos(math.radians((cell_angles[2])))  ) / (math.sin(math.radians(cell_angles[2])) ) ) * atom_fz[0]
            z[0]    =    ( cell_positions[0] *  cell_positions[1]  *  cell_positions[2]* (1-(math.cos(math.radians((cell_angles[0]))))**2 - (math.cos(math.radians((cell_angles[1]))))**2-(math.cos(math.radians((cell_angles[2]))))**2 + 2* math.cos(math.radians((cell_angles[0]))) * math.cos(math.radians((cell_angles[1]))) * math.cos(math.radians((cell_angles[2]))))**(1/2)/(cell_positions[0] * cell_positions[1] * math.sin(math.radians((cell_angles[2])))  )) * atom_fz[0]
            
            #3  Calculate atom positions in cartesian format.  Outputs:  x[1], y[1], z[1] - x[atom_count], y[atom_count], z[atom_count]
            n = 1
            n_max = atom_count
            while n <= n_max:
                x[n]    =    cell_positions[0] * atom_fx[n]  +  ( cell_positions[1] * math.cos(math.radians((cell_angles[2])))  ) * atom_fy[n]    +    ( cell_positions[2] * math.cos(math.radians((cell_angles[1]))) ) * atom_fz[n] 
                y[n]    =    ( cell_positions[1] * math.sin(math.radians((cell_angles[2])))) * atom_fy[n]  +  ( cell_positions[2] * (math.cos(math.radians((cell_angles[0])))  - math.cos(math.radians((cell_angles[1]))) * math.cos(math.radians((cell_angles[2])))  ) / (math.sin(math.radians(cell_angles[2])) ) ) * atom_fz[n]
                z[n]    =    ( cell_positions[0] *  cell_positions[1]  *  cell_positions[2]* (1-(math.cos(math.radians((cell_angles[0]))))**2 - (math.cos(math.radians((cell_angles[1]))))**2-(math.cos(math.radians((cell_angles[2]))))**2 + 2* math.cos(math.radians((cell_angles[0]))) * math.cos(math.radians((cell_angles[1]))) * math.cos(math.radians((cell_angles[2]))))**(1/2)/(cell_positions[0] * cell_positions[1] * math.sin(math.radians((cell_angles[2])))  )) * atom_fz[n]
                n += 1
            
            #4    Calculate atom positions in all 27 cells in cartesian format.  Outputs:  shift_x[n], shift_y[n], shift_z[n]   NOTE:   data starts on [0], not [1] like above
            #  FORMULA:  index =  ( shift * atom_count ) + atom_no - 1  
                                  ## shift = 0 to 26
                                  ## n  =  atom_count (fixed #)
                                  ## atom_no = 1 to atom_count
            shift_x = [float(0)]*2000*27
            shift_y = [float(0)]*2000*27
            shift_z = [float(0)]*2000*27
            shift = 0
            shift_max = 27
            while shift < shift_max:
                n = 0
                n_max = atom_count
                while n < n_max:
                    if shift == 0:  
                        shift_x[atom_count*shift+n] = x[n+1]
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 1:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 2:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 3:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_yx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 4:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_yx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy + cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 5:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_yx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy - cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 6:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_yx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 7:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_yx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy + cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 8: 
                        shift_x[atom_count*shift+n] = x[n+1] - cell_yx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy - cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 9:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 10:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 11:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 12:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx + cell_yx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 13:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx + cell_yx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 14:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx + cell_yx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 15:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx - cell_yx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 16:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx - cell_yx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy + cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 17:
                        shift_x[atom_count*shift+n] = x[n+1] + cell_xx - cell_yx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy - cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 18:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 19:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 20:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1]
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 21:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx + cell_yx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 22:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx + cell_yx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy + cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 23:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx + cell_yx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] + cell_yy - cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz
                    if shift == 24:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx - cell_yx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy
                        shift_z[atom_count*shift+n] = z[n+1]
                    if shift == 25:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx - cell_yx + cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy + cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] + cell_zz
                    if shift == 26:
                        shift_x[atom_count*shift+n] = x[n+1] - cell_xx - cell_yx - cell_zx
                        shift_y[atom_count*shift+n] = y[n+1] - cell_yy - cell_zy
                        shift_z[atom_count*shift+n] = z[n+1] - cell_zz       
                    n += 1
                shift+=1

            #5    Calculate % of vdW distances.  For each atom in unit cell, calculate distance to every other atom (whether in unit cell or not).  This creates a 2D matrix of the two indices.
            #  initialize distance_matrix
            distance_matrix = [[0]*(atom_count*27) for m in range(atom_count)]     #NOTE:   the index on left (n) here is on right (n) below.  m = x,y,z;  n = shift_x,y,z
            
            element_1_matrix = [[0]*(atom_count*27) for m in range(atom_count)]
            element_2_matrix = [[0]*(atom_count*27) for m in range(atom_count)]
            
            m=0
            m_max = atom_count
            n=0
            n_max = atom_count*27
            while m < m_max:
                while n < n_max:
                    n_new = n%atom_count+1
                    distance_matrix[m][n] =  ((shift_x[n] - x[m+1])**2 + (shift_y[n] - y[m+1])**2  + (shift_z[n] - z[m+1])**2)**(1/2) / (float(radius[element.index(atom_id[m+1])]) + float(radius[element.index(atom_id[n_new])]))
                    element_1_matrix[m][n] = element.index(atom_id[m+1])
                    element_2_matrix[m][n] = element.index(atom_id[n_new])
                    n += 1
                n=0
                m += 1
            #####################################
            
            m=0
            
            #decide how many distances to check.  54 in general, unless just one atom in unit cell
            stop_at = 54
            if atom_count == 1:
                stop_at = 27
            
            index_list = []    #creates a place to save indices [n]
            
            output_distances = str('')
            
            while m<m_max:
                distance_list = distance_matrix[m]  #extract row m from matrix 
                element_1_list = element_1_matrix[m]
                element_2_list = element_2_matrix[m]
                sorted_distance_list = sorted(distance_list) #organize distances from low to high
                sorted_distance_list = sorted_distance_list[:stop_at] #reduce size of matrix to 27 or 54 items
                sorted_distance_list_save = sorted_distance_list
                sorted_distance_list = list(OrderedDict.fromkeys(sorted_distance_list)) #remove duplicates
                
                #we now have a list of the shorted unique distances from atom m
                # need to find indices [n] for those nearest neighbors
                                
                for x in sorted_distance_list:
                    for i, j in enumerate(distance_list):
                        if j == x:
                            index_list.extend([i])
                
                # at this point, we now have "index_list", which contains 54 (or more) indices of [n] for a given [m]
                # this can be more than 54 if the final search returns several atoms of the same distance
                # for simplicity, I want to trim the list to 54 items
                
                index_list = index_list[:54]
                
                # now we can use the values of [n] to get the coordinates from shift_x, shift_y, and shift_z
                
                index_list_counter = 0
                ################
                
                
                
                
                #################
                
                x_checker = []
                y_checker = []
                z_checker = []
                element_1_checker = []
                element_2_checker = []
                
                element_1_checker_list = []
                element_2_checker_list = []
                
                angles = []
                angle_min = []
                
                limit = 54
                
                if atom_count == 1:
                    limit = 27
                
                while index_list_counter < limit:
                    if index_list_counter == 0:
                        x_checker.extend([shift_x[index_list[index_list_counter]]])
                        y_checker.extend([shift_y[index_list[index_list_counter]]])
                        z_checker.extend([shift_z[index_list[index_list_counter]]])
                    if index_list_counter > 0:
                        x_checker.extend([shift_x[index_list[index_list_counter]]-x_checker[0]])
                        y_checker.extend([shift_y[index_list[index_list_counter]]-y_checker[0]])
                        z_checker.extend([shift_z[index_list[index_list_counter]]-z_checker[0]])
                    
                    #this builds a list of the 54 or 27 elements that are nearest to the atom in question.  this list is organized from shortest to farthest.
                    element_1_checker.extend([element[element_1_list[index_list[index_list_counter]]]])
                    element_2_checker.extend([element[element_2_list[index_list[index_list_counter]]]])
                    
                    # now, in each x_ y_, z_checker, the [1]-[53] spots are filled with corrected coordinates... i.e., referenced to an origin of 0,0,0
                    
                    angle_checker = 1
                    
                    if index_list_counter > 1:  #we need at least three points to calculate angles
                        
                        while angle_checker < index_list_counter:  #compare new point to each proceeding point
                            angles.extend([math.degrees(math.acos(round((((x_checker[index_list_counter])*(x_checker[angle_checker])+(y_checker[index_list_counter])*(y_checker[angle_checker])+(z_checker[index_list_counter])*(z_checker[angle_checker]))/((math.sqrt(x_checker[index_list_counter]**2+y_checker[index_list_counter]**2+z_checker[index_list_counter]**2))*(math.sqrt(x_checker[angle_checker]**2+y_checker[angle_checker]**2+z_checker[angle_checker]**2)))), 6)))])
                            angle_checker += 1
                                
                        angle_min.extend([min(float(s) for s in angles)]) #stores only the smallest angle for new bond (angle btw new bond and all previous bonds)
                
                    index_list_counter += 1

                #we now output the angles for each bond, recording the smallest angle to neighbors.  let's pick a threshold value for the angle to include as a bond: ____ degrees
                # if angle in angle_min is less than ____ degrees, then record the index and acquire the distance
                #note, we have 52 entries in each angle_min.  to translate back to sorted_distance_list_save 
                
                angle_min_list = [] #this contains an index of atoms that meet the angle criterion
                
                element_1_checker_list.append(element_1_checker[0])
                element_2_checker_list.append(element_2_checker[1])
                
                
                #here, I want to implement an angle cut-off that depends on distance.  I will a variable, cut_off_angle_checker to depend on the distance.  The distance is checked from

                cut_off_angle_calculated = []

                for i, j in enumerate(sorted_distance_list_save):
                    #need to plug in distance and output a cut-off angle
                    if j < angle_distance_lower_bound:
                        cut_off_angle_calculated.extend([cut_off_angle])
                    if angle_distance_lower_bound <= j < angle_distance_upper_bound:
                        cut_off_angle_calculated.extend([(j-angle_distance_lower_bound)*(cut_off_angle_end-cut_off_angle)/(angle_distance_upper_bound-angle_distance_lower_bound)])
                    if j >= angle_distance_upper_bound:
                        cut_off_angle_calculated.extend([cut_off_angle_end])

                for i, j in enumerate(angle_min):
                    if j > cut_off_angle_calculated[i]:
                        angle_min_list.extend([i])
                        element_1_checker_list.extend([element_1_checker[i+1]]) #record element if the angle is kept
                        element_2_checker_list.extend([element_2_checker[i+1]]) #record element if the angle is kept
                
                angle_min_list.insert(0, -1) #add
                angle_min_list = [x+2 for x in angle_min_list]

                final_distance_list = []
                
                for x in angle_min_list:
                    if examine_all_bonds == 1:
                        final_distance_list.extend([sorted_distance_list_save[x]])
                    else:
                        if target_atom_1 == element_1_checker_list[x-1]:
                            if target_atom_2 == element_2_checker_list[x-1]:
                                final_distance_list.extend([sorted_distance_list_save[x]])
                        elif target_atom_1 == element_2_checker_list[x-1]:
                            if target_atom_2 == element_1_checker_list[x-1]:
                                final_distance_list.extend([sorted_distance_list_save[x]])
                    
                    
                
                output_distances = output_distances + ' ' + " ".join(str(x) for x in final_distance_list)

                
                del angle_min
                angle_min = []
                
                #now that we have 3 lists of x, y, & z coordinates, we need to calculate angles.
                
                
                
                del index_list
                index_list = []
                m+=1
                #get index from sorted_distance_list
                

            
            # (2) check bond distance and angle against each earlier neighbor
            # (3) If distance / angle is longer than threshold, do not include bond.  Otherwise, record distance.
            
            #####################################
            #6   Write out results
            

            
            z_count = 0
            m=1
            m_max = atom_count
            n=0
            n_max = atom_count*27
            
            text_out = str('')
            os.chdir(write_path)
            
            text_out = text_out + output_distances + ' '
            

            if counter < 875:
                out = open('vdw_0.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 1750:
                out = open('vdw_1.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 2625:
                out = open('vdw_2.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 3500:
                out = open('vdw_3.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 4375:
                out = open('vdw_4.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 5250:
                out = open('vdw_5.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 6125:
                out = open('vdw_6.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 7000:
                out = open('vdw_7.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 7875:
                out = open('vdw_8.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 8750:
                out = open('vdw_9.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 9625:
                out = open('vdw_10.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 10500:
                out = open('vdw_11.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 11375:
                out = open('vdw_12.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 12250:
                out = open('vdw_13.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 13125:
                out = open('vdw_14.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 14000:
                out = open('vdw_15.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 14875:
                out = open('vdw_16.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 15750:
                out = open('vdw_17.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 16625:
                out = open('vdw_18.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 17500:
                out = open('vdw_19.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 18375:
                out = open('vdw_20.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 19250:
                out = open('vdw_21.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 20125:
                out = open('vdw_22.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 21000:
                out = open('vdw_23.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 21875:
                out = open('vdw_24.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 22750:
                out = open('vdw_25.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 23625:
                out = open('vdw_26.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 24500:
                out = open('vdw_27.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 25375:
                out = open('vdw_28.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 26250:
                out = open('vdw_29.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 27125:
                out = open('vdw_30.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 28000:
                out = open('vdw_31.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 28875:
                out = open('vdw_32.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 29750:
                out = open('vdw_33.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 30625:
                out = open('vdw_34.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 31500:
                out = open('vdw_35.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 32375:
                out = open('vdw_36.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 33250:
                out = open('vdw_37.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 34125:
                out = open('vdw_38.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 35000:
                out = open('vdw_39.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 35875:
                out = open('vdw_40.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 36750:
                out = open('vdw_41.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 37625:
                out = open('vdw_42.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 38500:
                out = open('vdw_43.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 39375:
                out = open('vdw_44.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 40250:
                out = open('vdw_45.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 41125:
                out = open('vdw_46.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 42000:
                out = open('vdw_47.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 42875:
                out = open('vdw_48.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 43750:
                out = open('vdw_49.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 44625:
                out = open('vdw_50.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 45500:
                out = open('vdw_51.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 46375:
                out = open('vdw_52.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 47250:
                out = open('vdw_53.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 48125:
                out = open('vdw_54.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 49000:
                out = open('vdw_55.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 49875:
                out = open('vdw_56.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 50750:
                out = open('vdw_57.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 51625:
                out = open('vdw_58.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 52500:
                out = open('vdw_59.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 53375:
                out = open('vdw_60.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 54250:
                out = open('vdw_61.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 55125:
                out = open('vdw_62.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 56000:
                out = open('vdw_63.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 56875:
                out = open('vdw_64.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 57750:
                out = open('vdw_65.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 58625:
                out = open('vdw_66.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 59500:
                out = open('vdw_67.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 60375:
                out = open('vdw_68.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 61250:
                out = open('vdw_69.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 62125:
                out = open('vdw_70.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 63000:
                out = open('vdw_71.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 63875:
                out = open('vdw_72.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 64750:
                out = open('vdw_73.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 65625:
                out = open('vdw_74.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 66500:
                out = open('vdw_75.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 67375:
                out = open('vdw_76.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 68250:
                out = open('vdw_77.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 69125:
                out = open('vdw_78.txt',  'a')
                out.write(text_out)
                out.close()
            elif counter < 70000:
                out = open('vdw_79.txt',  'a')
                out.write(text_out)
                out.close()
            text_out = str('')
            n=0
            m += 1
        print(str(round((counter/cif_counter*100), 2)) + '% done')
        os.chdir(read_path)
