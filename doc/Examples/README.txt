AMBiT example calculations
==========================

This directory contains a set of small, runnable input files chosen to span the
features users typically need first. Each example is in its own subdirectory and
documented in-file with header comments explaining what it demonstrates.

Run any example with:

    ambit <subdir>/<atom>.input

For the larger ones (CrII, CuIII) you will want OpenMP enabled at build time and
a few cores available; the smaller ones (HeI, NaI) finish quickly on a laptop.


Picking the right example for what you want to learn
----------------------------------------------------

If you want to see ...                            ... look at

  the very simplest CI possible (2 electrons,     HeI/He.input
  no MBPT, V^0 starting potential, E1 line
  strengths from pure CI)

  a single valence electron above a closed core,  NaI/Na.input
  Brueckner MBPT, hyperfine constants with
  RPA core polarisation (V^{N-1})

  two valence electrons above a closed core,      MgI/Mg.input
  V^{N-2} starting potential, CI + MBPT, multiple
  even/odd J^pi solved together

  an open d-shell with several leading            CrII/Cr+.input
  configurations (3d^5, 3d^4 4s, 3d^4 4p), the
  Emu-CI 'SmallSide' construction, M1 and E2
  transitions between named levels

  hole excitations from a closed shell (3d^-1     CuIII/Cu.input
  rather than 3d^9), FrozenCore, mixed
  ElectronExcitations + HoleExcitations


Feature-by-feature index
------------------------

  Starting potential
    V^0 (bare nucleus, all electrons in CI) ...... HeI
    V^{N-1} (alkali, one valence electron) ....... NaI
    V^{N-2} (two valence electrons) .............. MgI
    Open d-shell V^{N} ........................... CrII
    Closed shell + holes ......................... CuIII

  CI features
    Pure CI, no MBPT ............................. HeI
    Multiple leading configurations .............. CrII
    Emu CI ([CI/SmallSide]) ...................... CrII
    Hole excitations + FrozenCore ................ CuIII
    Both parities solved together ................ HeI, NaI, MgI

  MBPT features
    Diagrammatic (no Brueckner) .................. MgI, CrII, CuIII
    Brueckner orbitals (--brueckner) ............. NaI

  Transitions
    E1 with RPA .................................. NaI
    E1 from pure CI (no RPA) ..................... HeI
    E2 + M1 between named levels ................. CrII
    Magnetic-dipole hyperfine (HFS1) with RPA .... NaI


What is intentionally NOT covered here
--------------------------------------

These examples deliberately stay small. The user guide (doc/UserGuide/) covers
the larger topics that don't fit a quick example:

  * MPI / ScaLAPACK parallelisation of the CI matrix
  * Dielectronic recombination ([DR] section)
  * Multi-run parameter sweeps (now handled by shell scripts or HPC array jobs,
    see Section "Multiple runs" of the user guide)
  * Lorentz-invariance, NMS/SMS isotope shift, Yukawa, and field-shift operators

For the exhaustive list of input options, consult Chapter "User input options"
of doc/UserGuide/user_guide.tex.
