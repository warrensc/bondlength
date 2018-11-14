import os
import re

start_path = os.path.join(os.sep,  'Sync', 'IDrive-Sync', 'jemsData',  'Materials genome',  'P1_cif')
out_path = os.path.join(os.sep,  'Sync', 'IDrive-Sync', 'jemsData',  'Materials genome',  'P1_cif_out')

os.chdir (start_path)

a = float()
b = float()
c = float()
alpha = float()
beta = float()
gamma = float()
atom=[str('')]*2000
x=[float(0)]*2000
y=[float(0)]*2000
z=[float(0)]*2000
last_n = int()

loop_line = int()

for file in os.listdir(start_path):
    n=0
    os.chdir(start_path)
    f = open(file, 'r')
    for num,  line in enumerate(f, 0):
        if line.startswith('_cell_length_a'):
            line_all = str(re.sub(' +',' ', line))
            line_parts = line_all.split(' ')
            a = float(line_parts[1]) 
        if line.startswith('_cell_length_b'):
            line_all = str(re.sub(' +',' ', line))
            line_parts = line_all.split(' ')
            b = float(line_parts[1])
        if line.startswith('_cell_length_c'):
            line_all = str(re.sub(' +',' ', line))
            line_parts = line_all.split(' ')
            c = float(line_parts[1])
        if line.startswith('_cell_angle_alpha'):
            line_all = str(re.sub(' +',' ', line))
            line_parts = line_all.split(' ')
            alpha = float(line_parts[1]) 
        if line.startswith('_cell_angle_beta'):
            line_all = str(re.sub(' +',' ', line))
            line_parts = line_all.split(' ')
            beta = float(line_parts[1])
        if line.startswith('_cell_angle_gamma'):
            line_all = str(re.sub(' +',' ', line))
            line_parts = line_all.split(' ')
            gamma = float(line_parts[1])
        if line.startswith('loop_'):
            loop_line = num
        if num > loop_line+2:
            if not line.startswith('_'):
                if line[:1].isalpha():
                    line_all = str(re.sub(' +',' ', line))
                    line_parts = line_all.split(' ')
                    atom[n]=line_parts[1]
                    x[n]=line_parts[2]
                    y[n]=line_parts[3]
                    z[n]=line_parts[4]
                    last_n = n
                    n+=1
    n=0
    # we  now have a, b, c // alpha, beta, gamma // atom, x y z
    text = str(last_n+1) + '\n' + str(a) + ' ' + str(b) + ' ' + str(c) + '\n'
    text = text + str(alpha) + ' ' + str(beta) + ' ' + str(gamma) + '\n'
    while n <= last_n:
        text = text + atom[n] + ' ' + x[n] + ' ' + y[n] + ' ' + z[n]
        if n != last_n:
            text = text + '\n'
        n+=1
    n=0
    os.chdir (out_path)
    out = open(file,  'w+')
    out.write(text)
    out.close()
