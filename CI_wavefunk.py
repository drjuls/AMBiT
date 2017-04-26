#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plot
import read_basis as rb
import read_hamiltonian as rh

class HamSolver:

  def __init__(self, fname):
    # First, grab the Hamiltonian matrix file we want to read in
    print("Reading Hamiltonian matrix from file: {}...\n".format(fname))
    self.ham = rh.HamiltonianMatrix(fname)

    # Now solve to get the eigenvectors and values of the Hamiltonian
    self.Es, self.Cs = self.ham.solve()

  def read_potential(self, fname):
    """ Reads the direct HF potential from the specified file and constructs a HF potential, with
    centrifugal barrier term."""

    self.V = []
    l = 3

    # Loop over the file and read its contents into our prepared lists. The file should be structured
    # like:
    #
    # r[i]  f[i]  dfdr[ii]
    #
    # And has no header. Note that we don't need to store r since that's already part of the basis
    # orbital data
    with open(fname, 'r') as fp:
      for line in fp:
        data = line.split()
        r = float(data[0])

        self.V.append(-float(data[1]) + l*(l+1)/(2*r**2))

  def get_rel_coeffs(self):

    # Get the 4f coefficient first, since that's the simplest (1CSF)
    eta = self.Cs[-1,0]
    #eta = -self.Cs[-1,0]

    # And do the 5f, which will be the first 4 elements in the lowest eigenvector
    fst4 = self.Cs[0:4,0]
    epsilon = np.sqrt(sum(fst4**2))
    # Now run through the remaining chunks of the eigenvector and print out those contributions
    deltas = []
    for ii in range(4, self.Cs[:,0].size - 1, 4):
      n = 5 + ii//4
      coeffs = self.Cs[ii:ii+4, 0] 

      deltas.append(np.sqrt(sum(coeffs**2)))

    # Finally, return an ndarray containing all of the relativistic config coefficients
    return(np.concatenate(([eta], [epsilon], deltas)))


  def print_rel_coeffs(self):
    # Factor out the leading coefficients epsilon and deltas to get the coefficients for each CSF
    # Get the 4f coefficient first, since that's the simplest (1CSF)
    print("Ground-state coefficients:")
    eta = self.Cs[-1,0]
    print("C_4f = {}\n===".format(eta))

    # And do the 5f, which will be the first 4 elements in the lowest eigenvector
    fst4 = self.Cs[0:4,0]
    epsilon = np.sqrt(sum(fst4**2))
    CSF_5f = fst4/epsilon
    print("C_5f = {}".format(epsilon))
    print("5f CSFs = {}\n===".format(CSF_5f))

    # Now run through the remaining chunks of the eigenvector and print out those contributions
    deltas = []
    for ii in range(4, self.Cs[:,0].size - 1, 4):
      n = 5 + ii//4
      coeffs = self.Cs[ii:ii+4, 0] 

      deltas.append(np.sqrt(sum(coeffs**2)))
      CSF_nf = coeffs/deltas[-1]
      print("C_{}f = {}".format(n, deltas[-1]))
      print("{}f CSFs = {}\n===".format(n, CSF_nf))

    # Sanity check: all the coeffecients' mod-squares should sum to 1
    total = eta**2 + epsilon**2 + sum([d**2 for d in deltas])
    print("\nTotal probability = {}".format(total))
  

  def build_wavefunk(self, basisFile):

    # First, grab the coefficients for each relativistic configuration in our solution (e.g. 4f-2 6s2)
    relCoeffs = self.get_rel_coeffs()

    # Make sure to calculate the maximum F state in our basis
    maxF = 3 + len(relCoeffs)

    
    # Now initialise a basis orbital object so we can use these coefficients to build up a ground-state
    # like wavefunction
    basisOrbs = rb.OrbitalMap(basisFile, silent = True) 

    # Now we want to run through an build up a wavefunction for the solution using the CI expansion 
    # coefficients. There's only one CSF per nonrelativistic configuration and each config has the form 
    # 4f-2 6s2 nf1 for n > 1 (note: the first condition only holds because we hacked AMBiT to ignore 
    # 4f*-1; J = 5/2 contributions). The first two terms don't change between configurations, so we're 
    # able to simply build up a ``4f-like'' wavefunction via a sum over the nf basis orbitals.
    labels = [str(n) + 'f' for n in range(4, maxF + 1)] # List of all labels in the CI expansion
    coeffMap = dict(zip(labels, relCoeffs)) # Map a CSF label to its coefficient

    
    wavefunkf = np.zeros(basisOrbs.latticeSize, dtype = np.float64)
    wavefunkg = np.zeros(basisOrbs.latticeSize, dtype = np.float64)
    for label in labels:
      orb = basisOrbs[label]
      wavefunkf[0:orb.size] += coeffMap[label] * orb.f
      wavefunkg[0:orb.size] += coeffMap[label] * orb.g
    
    # Now we want to plot the constructed 4f wavefunction, as well as the 4f orbital and the ``5f''
    # correction from InjectOrbitals
    plot.plot(basisOrbs.r, wavefunkf,
        label="Aritificial CI - Max {}f: Psi(r)".format(maxF))
    plot.plot(basisOrbs.r[0:basisOrbs["4f"].size], basisOrbs["4f"].f,
        label="HF: 4f(r)")
    plot.plot(basisOrbs.r[0:basisOrbs["5f"].size], basisOrbs["5f"].f,
        label="HF: 5f(r)")
    plot.plot(basisOrbs.r[0:basisOrbs["6f"].size], basisOrbs["6f"].f,
        label="HF: 6f(r)")

    # Also plot the direct potential (if we've read it)
    if(self.V != None):
      plot.plot(basisOrbs.r, self.V, label="V(direct) + centrifugal", linestyle="--", color=(0,0,0))

    plot.xlabel("r (a.u.)")
    plot.ylabel("wavefunction")
    plot.legend(loc="lower right")
    ax = plot.gca()
    ax.set_xlim([1e-1, 30])
    ax.set_xscale('log')
    ax.set_ylim([-7.5, 2.0])
    plot.show() 
   
  def find_turning_points(self, E, tol = 1e-4):
    """ Finds the classical turning point of the direct HF potential (i.e. where E - V(direct) < tol). 
    
    Returns a list containing the indices of the turning points. Since we're dealing with floating point
    arithmetic, we can only find the points that are within some specified tolerance of the true turning
    point. Currently this defaults to 1e-4, but may need to be tweaked on a per-atom basis."""
    # Use numpy.where() to get the indices of the potential array where V = E. This has an
    # extraordinarily messy syntax, but it's probably better than just looping over the whole damn thing
    # manually

    indices = np.where(abs(np.array(self.V) - E) < tol)[0]
    return(list(indices))

