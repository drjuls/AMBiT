#!/usr/bin/python3
""" Reads an (binary) AMBiT basis data file and analyse the basis orbitals, similar to AMBiT.m. 

This module is primarily meant to be used interactively (i.e. import it and play around with the 
orbitals in the interpreter), and contains the following classes:

Orbital - Class representing all the information for a single basis orbital. It contains the upper 
and lower components of the basis spinor (i.e. f and g) and their derivatives, as well as its energy,
principal quantum number and relativistic quantum number, kappa. Each orbital may contain fewer (but 
never more) points than the global lattice.

OrbitalMap - Basically an overgrown dictionary to hold and manage all the stored basis orbitals. The 
orbitals are indexed by strings containing their spectroscopic label (e.g. "4f+"); individual
Orbital objects stored by an instance of OrbitalMap can be directly accessed via the overloaded 
[] operator. The class also contains the global lattice properties: the number of points in the 
lattice, as well as ndarrays containing the r and dr points. Finally, this class also contains 
functions to analyse and visualise the orbital functions.

Finally, this module makes heavy use of the numpy and matplotlib libraries for data analysis and 
plotting, so some familiarity with these is recommended to get the most out of this module.

"""
import numpy as np
import matplotlib.pyplot as plot


class Orbital:
  """Class representing all the information for a single basis orbital. 
  
  This class contains the upper and lower components of the basis spinor (i.e. f and g) and their 
  derivatives, as well as its energy, principal quantum number and relativistic quantum 
  number, kappa. Each orbital may contain fewer (but never more) points than the global lattice. 

  This class contains to following functions:

  __init__ - initialises an orbital function with the quantum numbers, energy and size given
             by the input arguments.
  set_spinor - Sets the ndarrays containing the radial components (f, g, dfdr, dgdr) of the 
               orbital from an already open binary file. This is meant only to be called from  
               the OrbitalMap  initialiser.
  get_density - Returns an nd array containing the radial probability density of the orbital
  """
  def __init__(self, pqn, kappa, E, size):

    # Set the essential parameters of the orbital
    self.pqn = int(pqn)
    self.kappa = int(kappa)
    self.energy = float(E)
    self.size = int(size)

    # Finally, give the orbital a label
    specLabels = ['s', 'p', 'd', 'f', 'g', 'h']

    # Start with the principal quantum number, then add the angular momentum in
    # spectroscopic notation
    self.label = str(self.pqn)
    if(self.kappa >= 0):
      self.label += specLabels[self.kappa]
    else:
      self.label += specLabels[-self.kappa - 1]

    # Finally, append a '+' to indicate if kappa < 0
    if(self.kappa < -1):
      self.label += '+'

  def set_spinor(self, infile):
    """ Reads the values of f, g, dfdr and dgdr from a given (binary) input file. 
    
        Only self.size elements are stored, which may be less than the global lattice size
        
    """

    self.f = np.fromfile(infile, count = self.size, dtype = np.float64)
    self.g = np.fromfile(infile, count = self.size, dtype = np.float64)
    self.dfdr = np.fromfile(infile, count = self.size, dtype = np.float64)
    self.dgdr = np.fromfile(infile, count = self.size, dtype = np.float64)

  def get_density(self):
    """ Returns an ndarray containing the probability density of the orbital"""
    return(self.f**2 + self.g**2)

class OrbitalMap:
  """ Class to hold and manage all stored basis orbitals.

  This is basically an overgrown dictionary to hold and manage all the stored basis orbitals. The 
  orbitals are indexed by strings containing their spectroscopic label (e.g. "4f+"); individual
  Orbital objects stored by an instance of OrbitalMap can be directly accessed via the overloaded 
  [] operator. The class also contains the global lattice properties: the number of points in the 
  lattice, as well as ndarrays containing the r and dr points. Finally, this class also contains 
  functions to 
  analyse and visualise the orbital functions.

  This class contains the following functions:

  __init__ - Reads the binary orbtials file (basisFileName) and creates a map of all the stored 
             orbitals. This function is roughly equivalent to AmbitReadAll in AMBiT.m.
  
  __getitem__ - Overloads the [] operator to return a stored orbital, indexed by its spectroscopic
                label.
  list_orbitals - Prints the spectroscopic label and energy of all stored orbitals.

  plot_orbitals - Plots the f component of a list of orbitals. Prints a warning and returns nothing
                  if any of the supplied orbitals do not exist.
  plot-densities - Plots the radial probability density of a list of orbitals. Prints a warning and
                   returns nothing if any of the supplied orbitals do not exist.

  get_overlap - Returns the overlap between two stored orbitals.

  """
  def __init__(self, basisFileName, silent = False):
    """ Reads the binary orbtials file (basisFileName) and creates a map of all the stored orbitals
    
    Initialising an instance of orbital map reads the AMBiT basis file with the filename 
    basisFileName and constructs a map of all the stored orbitals, indexed by their spectroscopic
    label (e.g. 4f+).

    It also reads the total number of orbitals, as well as the global lattice properties, which 
    are stored in the ndarrays r and dr. Note that basisFileName should include the relative 
    path to the file, which shouldn't be already open when this constructor is invoked.

    This function is roughly equivalent to the AmbitReadAll function in AMBiT.m.

    """

    with open(basisFileName, 'rb') as ifp:
      # First, grab the header line, which contains the lattice parameters with
      # the structure [double, double, double, int]
      header_t = np.dtype([('beta', np.float64), ('h', np.float64), 
          ('rmin', np.float64), ('latticeSize', np.int32)])

      header = np.fromfile(ifp, count = 1, dtype = header_t)
      self.latticeSize = int(header['latticeSize'])

      # Next, read the lattice information (arrays containing r and dr)
      self.r = np.fromfile(ifp, count=self.latticeSize, dtype=np.float64)
      self.dr = np.fromfile(ifp, count=self.latticeSize, dtype=np.float64)

      # And grab the total number of orbitals
      self.numOrbitals = np.fromfile(ifp, count=1, dtype = np.int32)[0]

      # Now we can finally go through and read the actual basis orbitals, which
      # are organised in the file like this (n.b. size is the size of the
      # orbital and may be less than the size of the lattice):
      #
      # pqn, energy, kappa, size
      # f[0], f[1], ..., f[size]
      # g[0], g[1], ..., g[size]
      # dfdr[0], ... dfdr[size]
      # dgdr[0], ... dgdr[size]
      #
      # Store the orbitals in a dict with tuples of (pqn, kappa) as keys
      self.__orbitals = {}
      for _ in range(self.numOrbitals):
        header_t = np.dtype([('pqn', np.int32), ('E', np.float64), 
            ('kappa', np.int32), ('size', np.int32)])
        header = np.fromfile(ifp, count = 1, dtype = header_t)
        orb = Orbital(header['pqn'], header['kappa'], header['E'],
            header['size'])
        orb.set_spinor(ifp)
        
        # Generate a label in spectroscopic notation so it's easier to identify
        # the levels
        pqn = int(header['pqn'])
        kappa = int(header['kappa'])
        self.__orbitals[orb.label] = orb

      # Finally, print out the number of orbitals we successfully read (if we're in verbose mode)
      if(not silent):
        self.list_orbitals()
 
  def __getitem__(self, label):
    """ Return the orbital object indexed by the given spectroscopic label (e.g. 1s, 4d+)."""
    return(self.__orbitals[label])

  def __setitem__(self, label, orb):
    """ Stores the given orbital object using the supplied key/label. """
    self.__orbitals[label] = orb

  def list_orbitals(self):
    """ Prints the spectroscopic label and energy of all stored orbitals."""
    print("{} basis orbitals stored:\n#Orbital  Energy".format(self.numOrbitals))

    for pair in sorted(self.__orbitals.items(), 
        key = lambda pair: pair[1].pqn):
      print("{}  {}".format(pair[0], pair[1].energy))
 
  def plot_orbitals(self, *args):
    """Prints the orbitals corresponding to a given list of spectroscopic labels. """
    
    # Exit if the user passed no arguments (for some reason)
    if(len(args) == 0):
      print("No orbitals given to plot_orbitals")
      return(None)

    # Now loop over the arguments and plot the orbitals
    for label in args:
      try:
        target = self.__orbitals[label]
      except KeyError:
        print("Error: could not find orbital {}".format(label))
        return(None)

      lattice = self.r[:target.size]

      plot.plot(lattice, target.f, label='{}: f(r)'.format(label)) 
    
    plot.xlabel('r (atomic units)')
    plot.ylabel("Radial Wavefunction")
    plot.legend(loc="upper right")
    plot.show()

  def plot_densities(self, *args):
    """ Prints the probability density of a given list of orbitals."""
    # Exit if the user passed no arguments (for some reason)
    if(len(args) == 0):
      print("No orbitals given to plot_densities")
      return(None)

    # Now loop over the arguments and plot the orbitals
    for label in args:
      try:
        target = self.__orbitals[label]
      except KeyError:
        print("Error: could not find orbital {}".format(label))
        return(None)

      lattice = self.r[:target.size]

      plot.plot(lattice, target.get_density(), label='{}: P(r)'.format(label)) 
      
    plot.xlabel('r (atomic units)')
    plot.ylabel("Probability density")
    plot.legend(loc="upper right")
    plot.show()

  def get_overlap(self, leftLabel, rightLabel):
    """ Returns the overlap between the two supplied orbitals."""

    # First, grab the two orbitals corresponding to the labels and grab the
    # size of the smallest one
    leftOrb = self.__orbitals[leftLabel]
    rightOrb = self.__orbitals[rightLabel]

    resultSize = min(leftOrb.size, rightOrb.size)
    integrand = (leftOrb.f[0:resultSize]*rightOrb.f[0:resultSize] 
        + leftOrb.g[0:resultSize]*rightOrb.g[0:resultSize])

    resultF = np.trapz(integrand, self.r[0:resultSize])
    return(resultF)

