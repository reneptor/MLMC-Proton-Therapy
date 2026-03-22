import numpy as np
import os

array = np.zeros(shape=(100, 100, 500), dtype=np.int16)

output_dir = os.path.join('..', 'input')  # This creates the path '../input'
file_path = os.path.join(output_dir, 'phantom.dat')
array.tofile(file_path)
