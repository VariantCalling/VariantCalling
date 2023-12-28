from VCDataGenerator import generate_data_for_noise_reduction, generate_random_clones, simulate_read
import time
import numpy as np
import matplotlib.pyplot as plt
import VariantCallingMini as vc

# OOP Python
dg = vc.VariantCallingDataMini()
start = time.time()
noisy, clean, _ = dg.generate_data_for_noise_reduction(1000, 100, 0.01, 0.05)
print("OOP Python time taken: ", time.time() - start)

plt.imshow(noisy[0], cmap='jet')
plt.show()


plt.imshow(clean[0], cmap='jet')
plt.show()


# Cython
start = time.time()
noisy, clean = generate_data_for_noise_reduction(1000, alignment_error_prob=0.01, sequencing_error_prob=0.05)
print("Cython time taken: ", time.time() - start)

plt.imshow(noisy[0], cmap='jet')
plt.show()


plt.imshow(clean[0], cmap='jet')
plt.show()