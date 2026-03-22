import numpy as np
import os
import SimpleITK as sitk
import matplotlib.pyplot as plt

CT_path = os.path.join('..', 'CT_phantom', 'CT.mhd')
CT_image = sitk.ReadImage(CT_path)
CT_array = sitk.GetArrayFromImage(CT_image)

print("CT_array shape:", CT_array.shape) # (192, 191, 155)

CT_flatten = CT_array.flatten()

print(np.max(CT_flatten)) # 3071
print(np.min(CT_flatten)) # -1024

# CT_filtered = CT_flatten[CT_flatten > 1000]

# print(CT_array.dtype)


# plt.hist(CT_filtered, bins=200)
# plt.xlabel('Pixel Intensity')
# plt.savefig("CT_histogram.png")


# output_dir = os.path.join('..', 'input')  # This creates the path '../input'
# file_path = os.path.join(output_dir, 'CT.dat')
# CT_array.tofile(file_path)