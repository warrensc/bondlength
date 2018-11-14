import os
import numpy as np
import glob

read_path = os.path.join(os.sep,  'Users', 'Admin', 'Desktop', 'jemsData',  'Materials genome',  'BondLength-16')

os.chdir(read_path)
file_counter = len(glob.glob1(read_path,"*.txt"))
file_no = 0

lower_range = 0
upper_range = 3
step_size = 0.01
bin_lower = lower_range + (step_size/2)
bin_upper = upper_range - (step_size/2)

bins = []*999
bins = np.arange(bin_lower, bin_upper, step_size)
bins=np.reshape(bins, (-1, 1))
print(str(np.shape(bins)))
bins = np.hsplit(bins, 1)
print(str(np.shape(bins[0])))

hist_data = [0]*999
hist_data[0]=bins[0]

np.savetxt('histogram_0.yes_csv', bins[0])

while file_no < file_counter:
    file_name = 'vdw_' + str(file_no) + '.txt'
    print('loading ' + file_name)
    a = np.genfromtxt(file_name, dtype=np.float, delimiter=' ')
    print('read in ' + file_name)

    a = np.reshape(a, (1, -1))
    print('reshaped ' + file_name)

    a = a[~np.isnan(a)]
    print('removed NaN')

    hist,bins_out = np.histogram(a, bins=np.arange(lower_range,upper_range,step_size))
    print('did histogram')

    hist =np.reshape(hist, (-1, 1))
    print('reshaped histogram')
    hist = np.hsplit(hist, 1)

    np.savetxt('histogram_' + str(file_no + 1) + '.yes_csv', hist[0])
    file_no += 1
    
print('hstacked bins & hist')

#we now have hist_data[0], hist_data[1], and hist_data[2].  We want to concatenate these into one array.
#out = np.hstack(hist_data)


print('done.')
