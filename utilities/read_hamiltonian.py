#!/usr/bin/python3

import numpy as np
class HamiltonianMatrix:
  """ Class to read and manipulate a Hamiltonian matrix, as written by AMBiT"""
  def __init__(self, hamfile):
    """ Reads the binary file with filename hamfile and stores the contained Hamiltonian matrix.

        Initialiser takes in the filename of the binary file containing the Hamiltonian matrix and 
        stores the data it contains in an ndarry. The Hamiltonian matrix is symmetrical, with only 
        the lower triangle stored in the binary file, so we need to do some algebra to get the full
        matrix.
    """

    with open(hamfile, 'rb') as ifp:
      # First read the size of the matrix
      self.N = np.fromfile(ifp, count = 1, dtype = np.int32)[0]

      # Now initialise an NxN matrix full of zeros to hold our half-matrix
      self.H = np.zeros((self.N, self.N), dtype = np.float64)

      # Time to read the matrix elements from the binary file. The elements are laid out consecutively:
      #
      # H[0,0], H[1,0], H[1,1], H[2,0], ... H[N,N]
      #
      # So, to get the data for the i'th row, we want to read in N-i floats and store them in the slice
      # H[i, i:N].
      for ii in range(self.N):
        rowSize = ii + 1
        row = np.fromfile(ifp, count = rowSize, dtype = np.float64)
        self.H[ii, 0:ii+1] = row


      # Finally, recover the full matrix, store the diagonal and we're done
      self.H = self.H + self.H.transpose() - np.diagflat(self.H.diagonal())
      self.diag = self.H.diagonal()

  def solve(self):
    """ Solves for the eigenvalues/vectors of H and returns all energies and CI expansion coefficients.
    """
    Es, Cs = np.linalg.eigh(self.H)
    return(Es, Cs)

  def print_solutions(self, numSolns = 0, get_coeffs = False):
    """ Prints the energy and leading configuration solutions for the first numSolns CI solutions.

    This function prints the energy levels and either leading configuration percentages (the default)
    or CI expansion coefficients/amplitudes for the atom described by H. This 
    function treats numSolns similarly to AMBiT, in that it prints all solutions if numSolns is
    zero, negative or not supplied.
    """

    # Get the eigensolutions via numpy's special solver for Hermitian matrices
    Es, Cs = self.solve()
    confs = Cs**2

    if(numSolns <= 0):
      numToPrint = Es.size
    else:
      if(numSolns > Es.size):
        print("Warning: found fewer than the requested {} solutions. Printing all {} instead".format(numSolns, Es.size))
      
      numToPrint = min(Es.size, numSolns)

    for ii in range(numToPrint):
      print("===\nSolution {}:".format(ii))
      if(get_coeffs == True):
        print("E = {}\nLeading Configs:\n{}\n".format(Es[ii], Cs[:,ii]))
      else:
        print("E = {}\nBasis Coefficients:\n{}\n".format(Es[ii], confs[:,ii]))

  def crop_matrix(self, newN):
    """ Returns a newN X newN matrix by cropping the current Hamiltonian matrix."""
    return(self.H[0:newN, 0:newN])

