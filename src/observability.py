import numpy as np
from scipy.linalg import expm, inv, null_space, eig
from sklearn.preprocessing import normalize



def is_system_observable(A, H, domain="TD"):
    is_observable = False
    O = H
    n = A.shape[0]
    for i in range(1, n):
        O = np.vstack((O, np.dot(H, np.linalg.matrix_power(A, i))))

    if np.linalg.matrix_rank(O) == n:
        is_observable = True
    return is_observable, O

def recommended_sv_measurements(A, H, domain="TD"):
    # this function returns the index of state-variable (in X) related to which measurements should be added
    is_observable, O = is_system_observable(A, H)
    # normalize observability matrix
    On = normalize(np.asarray(O.T), axis=1, norm='l1')
    # OOt
    OOt = On @ On.T
    eigenvalues, eigenvectors = np.linalg.eig(OOt)
    # extract the index of the minimum eigenvalue (absolute)
    min_eigenvalue_index = np.argmin(np.abs(eigenvalues))
    # extract the eigenvector corresponding to the minimum eigenvalue
    min_eigenvector = eigenvectors[:, min_eigenvalue_index]
    # find the index which has the highest value in this eigenvector
    max_eigenvector_index = np.argmax(np.abs(min_eigenvector))

    print("top_3_meas_sv = ", np.argsort(np.abs(min_eigenvector))[-3:][::-1])

    return max_eigenvector_index

def observability_degree(A, H, system_type='Continuous'):
    is_observable, O = is_system_observable(A, H)
    On = normalize(np.asarray(O.T), axis=1, norm='l1')
    # OOt
    OOt = On @ On.T
    eigenvalues, eigenvectors = np.linalg.eig(OOt)
    # extract the index of the minimum eigenvalue (absolute)
    ratio = np.min(np.abs(eigenvalues)) / np.max(np.abs(eigenvalues))
    return ratio


def is_system_detectable(A, H, system_type = 'Continuous'):
    """
    Check if the system is detectable using Kalman decomposition.
    
    Parameters:
    A (ndarray): State matrix of the system.
    C (ndarray): Output matrix of the system.
    
    Returns:
    bool: True if the system is detectable, False otherwise.
    """

    is_detectable = False
    eigenvalues_td, eigenvectors_td = np.linalg.eig(A)
    is_observable, observability_matrix = is_system_observable(A, H)
    # Number of states
    n = A.shape[0]
    obs_num = np.linalg.matrix_rank(observability_matrix)
    unobs_num = n - obs_num
    
    # Find the null space of the observability matrix
    unobservable_subspace = null_space(observability_matrix)

    #print("null_subspace: ", unobservable_subspace, np.shape(unobservable_subspace))

    # If the null space is empty, the system is observable and hence detectable
    if unobservable_subspace.size == 0:
        is_detectable = True
    
    # Transform A to the unobservable subspace basis
    A_unobs = unobservable_subspace.T @ A @ unobservable_subspace
    

    # Check eigenvalues of the unobservable subspace
    eigenvalues, _ = eig(A_unobs)
    
    # Compute the row space (observable subspace)
    q, r = np.linalg.qr(observability_matrix.T, mode='complete')

    
    # Combine observable and unobservable basis to form a new state-space basis
    #T = np.hstack((observable_subspace, unobservable_subspace))
    T = q

    # Transform A and C matrices to the new basis
    T_inv = np.linalg.inv(T)
    A_transformed = T_inv @ A @ T
    C_transformed = H @ T
    
    A_unobs1 = A_transformed[-unobs_num:, -unobs_num:]
    eig_A_unobs1, _ = eig(A_unobs1)
    

    # Check stability of unobservable modes 
    for eigenvalue in eig_A_unobs1:
        if system_type == 'Continuous':
            if np.real(eigenvalue) >= 0:  # For continuous-time systems
                print("positive eig value: ", eigenvalue)
                is_detectable = False  # Unstable mode is unobservable, so system is not detectable
            else: 
                is_detectable = True
        elif system_type == 'Discrete':
            if np.abs(np.real(eigenvalue)) > 1:  # For discrete-time systems; if real part lies outside of unit circle
                is_detectable = False  # Unstable mode is unobservable, so system is not detectable
            else: 
                is_detectable = True

    
    return is_detectable  # All unobservable modes are stable, so system is detectable


