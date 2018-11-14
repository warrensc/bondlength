import os
import numpy as np
import glob

read_path = os.path.join(os.sep,  'Users', 'Admin', 'Desktop', 'jemsData',  'Materials genome',  'BondLength-16')

os.chdir(read_path)
file_counter = len(glob.glob1(read_path,"*.yes_csv"))
file_no = 1

f = np.genfromtxt('histogram_0.yes_csv')
f = np.reshape(f, (-1, 1))

while file_no < file_counter:
    name = 'histogram_' + str(file_no) + '.yes_csv'
    g = np.genfromtxt(name)
    g = np.reshape(g, (-1, 1))
    f = np.hstack((f, g))
    print(str(np.shape(f)))
    file_no += 1
    
print('here 1')

f = np.delete(f, 0, axis=1)

print('here 2')

sum = np.sum(f,  axis=1)
sum = np.reshape(sum, (-1, 1))

print('here 3')

f = np.genfromtxt('histogram_0.yes_csv')
f = np.reshape(f, (-1, 1))
f = np.hstack((f, sum))

np.savetxt('bond_length_histogram.csv', f, delimiter=',')

