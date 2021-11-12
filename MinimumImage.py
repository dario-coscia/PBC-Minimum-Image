# ============ ============= ============= ============= =============
# Python program to calculate minimum distance, using minimum image 
# convention for 2D (1,1,1) lattice. The ideces (m,n,l) refers to the
# Miller indeces used in crystallography. The algorithm projects the 
# Euclidean distance between the two particles into the two basis vectors, 
# then classical minimum image is applied. 
#
# To run the code simply uncomment the last part and type:
#
#       python3 MinimumImage.py
#
# If everything is correct you should get as output:
#
#       Minimum distance between the two particles: 0.0
#
# ============ ============= ============= ============= =============

import numpy as np

def MinimumImage(x, y, L1, L2):
    ''' 
    Function to calculate minimum distance between two numpy
        vectors x, y, in a 2D square lattice with basis vector
        two numpy vectors L1 and L2.
    '''
    
    dr = x - y
    
    beta = (L1[0]*dr[1]-L1[1]*dr[0])/(L1[0]*L2[1]-L1[1]*L2[0])
    alpha = (L2[1]*dr[0]-L2[0]*dr[1])/(L1[0]*L2[1]-L1[1]*L2[0])
    
    dr1 = alpha*L1
    dr2 = beta*L2 
    
    mod_1 = L1*round(np.linalg.norm(dr1)/np.linalg.norm(L1))
    mod_2 = L2*round(np.linalg.norm(dr2)/np.linalg.norm(L2))
    
    if np.linalg.norm(dr1 - mod_1) <= np.linalg.norm(dr1 + mod_1):
        dr1 = dr1 - mod_1
    else:
        dr1 = dr1 + mod_1   

    if np.linalg.norm(dr2 - mod_2) <= np.linalg.norm(dr2 + mod_2):
        dr2 = dr2 - mod_2
    else:
        dr2 = dr2 + mod_2

    return np.linalg.norm(dr1 + dr2)

# --- Example of program (uncomment to run it) --- #

# x = np.array([0, 0])
# y = np.array([np.sqrt(3), 0])
# L1 = np.array([np.sqrt(3)*0.5, 0.5])
# L2 = np.array([np.sqrt(3)*0.5, -0.5])
# sep = MinimumImage(x, y, L1, L2)
# print("Minimum distance between the two particles: ",sep)
