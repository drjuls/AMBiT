(* ::Package:: *)

BeginPackage["AMBiT`"];


AmbitReadLattice::usage = "AmbitReadLattice[\"filename.basis\"] returns list of points r used by AMBiT.";

AmbitReadBasis::usage = "AmbitReadBasis[\"filename.basis\"] returns {energies, orbitals}. These are Association objects with keys {pqn, kappa}.
e.g. energies[{1,-1}] returns energy of 1s orbital. orbitals[key]={f,g} where f,g are interpolated functions.";

AmbitReadDiscreteBasis::usage = "AmbitReadDiscreteBasis[\"filename.basis\"] returns {energies, orbitals}. These are Association objects with keys {pqn, kappa}.
e.g. energies[{1,-1}] returns energy of 1s orbital. orbitals[key]={f,g} where f,g are values at the points r.";

AmbitReadHamiltonian::usage = "AmbitReadHamiltonian[\"filename.matrix\"] returns symmetric real matrix."


AmbitMatrixElement::usage = "AmbitMatrixElement[orb1,orb2,V] = \!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(0\)], \(\[Infinity]\)]\)(\!\(\*SubscriptBox[\(f\), \(1\)]\)\!\(\*SubscriptBox[\(f\), \(2\)]\)+\!\(\*SubscriptBox[\(g\), \(1\)]\)\!\(\*SubscriptBox[\(g\), \(2\)]\))V(r)dr where {\!\(\*SubscriptBox[\(f\), \(1\)]\),\!\(\*SubscriptBox[\(g\), \(1\)]\)}=orb1 etc.";

AmbitHartreeY::usage = "AmbitHartreeY[k,{fa,ga},{fb,gb}] returns the Hartree function \!\(\*SubscriptBox[SuperscriptBox[\(Y\), \(k\)], \(a, b\)]\)=\[Integral]\!\(\*FractionBox[SuperscriptBox[SubscriptBox[\(r\), \(<\)], \(k\)], SuperscriptBox[SubscriptBox[\(r\), \(>\)], \(k + 1\)]]\)\!\(\*SubscriptBox[\(\[Psi]\), \(a\)]\)(r')\!\(\*SubscriptBox[\(\[Psi]\), \(b\)]\)(r')dr'. k is an natural number; {f,g} are interpolated functions (orbitals)."

AmbitRadialIntegral::usage = "AmbitRadialIntegral[k,{fa,ga},{fb,gb},{fc,gc},{fd,gd}] returns \!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(0\)], \(\[Infinity]\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(Y\), \(k\)], \(b, d\)]\)(r)\!\(\*SubscriptBox[\(\[Psi]\), \(a\)]\)(r)\!\(\*SubscriptBox[\(\[Psi]\), \(c\)]\)(r)dr"


Begin["`Private`"];

AmbitReadAll[filename_]:=Module[
{basisFileStream,orbitals,energies,beta,h,rmin,latticeSize,r,dr,numOrbitals},
basisFileStream = OpenRead[filename,BinaryFormat->True];
(* Read lattice parameters and lattice. *)
{beta,h,rmin,latticeSize}=BinaryRead[basisFileStream,{"Real64","Real64","Real64","Integer32"}];
r=BinaryReadList[basisFileStream,"Real64",latticeSize];
dr=BinaryReadList[basisFileStream,"Real64",latticeSize];
numOrbitals=BinaryRead[basisFileStream,"Integer32"];
orbitals=Association[];
energies=Association[];
(* Read all orbitals *)
For[i=0,i<numOrbitals,i++,
Block[{pqn,energy,kappa,size,f,g,dfdr,dgdr,padding},
{pqn,energy,kappa,size}=BinaryRead[basisFileStream,{"Integer32","Real64","Integer32","Integer32"}];
padding=ConstantArray[0.,latticeSize-size];
f=BinaryReadList[basisFileStream,"Real64",size]~Join~padding;
g=BinaryReadList[basisFileStream,"Real64",size]~Join~padding;
dfdr=BinaryReadList[basisFileStream,"Real64",size]~Join~padding;
dgdr=BinaryReadList[basisFileStream,"Real64",size]~Join~padding;
AppendTo[energies,{pqn,kappa}->energy];
AppendTo[orbitals,{pqn,kappa}->{f,g}];
]];
Close[basisFileStream];
Return[{beta,h,rmin,latticeSize,r,dr,energies,orbitals}];
]

AmbitReadLattice[filename_]:=Module[{beta,h,rmin,latticeSize,r,dr,energies,orbitals},
{beta,h,rmin,latticeSize,r,dr,energies,orbitals}=AmbitReadAll[filename];
Return[r];
]

AmbitReadDiscreteBasis[filename_]:=Module[{beta,h,rmin,latticeSize,r,dr,energies,orbitals},
{beta,h,rmin,latticeSize,r,dr,energies,orbitals}=AmbitReadAll[filename];
Return[{energies,orbitals}];
]

AmbitReadBasis[filename_]:=Module[
{beta,h,rmin,latticeSize,r,dr,energies,discreteOrbitals,F,G,interpolatedOrbitals},
{beta,h,rmin,latticeSize,r,dr,energies,discreteOrbitals}=AmbitReadAll[filename];
interpolatedOrbitals={Interpolation[{r,#[[1]]}\[Transpose]],Interpolation[{r,#[[2]]}\[Transpose]]}&/@discreteOrbitals;
Return[{energies,interpolatedOrbitals}];
]

AmbitReadHamiltonian[filename_]:=Module[{matrixFileStream,dim,Nsmall,hamiltonian},
matrixFileStream=OpenRead[filename,BinaryFormat->True];
(* First two integers are size *)
dim=BinaryRead[matrixFileStream,"Integer32"];
Nsmall=BinaryRead[matrixFileStream,"Integer32"];
(* Read lower triangle *)
hamiltonian=Table[BinaryReadList[matrixFileStream,"Real64",Min[i+1,Nsmall]],{i,0,dim-1}];
Close[matrixFileStream];
hamiltonian=PadRight[hamiltonian,{dim,dim}];
(* Return symmetric square matrix *)
Return[hamiltonian+hamiltonian\[Transpose]-DiagonalMatrix[Diagonal[hamiltonian]]];
]


AmbitMatrixElement[orb1_,orb2_,V_]:=Module[{f1,g1,f2,g2,rsmall,rlarge},
{f1,g1}=orb1;
{f2,g2}=orb2;
{{rsmall,rlarge}}=f1[[1]];
NIntegrate[(f1[r] f2[r]+g1[r]g2[r])V[r],{r,0,rlarge}]
]

AmbitHartreeY[k_,{fa_,ga_},{fb_,gb_}]:=Module[{rsmall,rlarge,density,theG,theY,H},
{{rsmall,rlarge}}=fa[[1]];
density[x_]:=fa[x] fb[x] + ga[x] gb[x];
theG = NDSolve[{G'[r]+(k+1)/r G[r]-density[r]/r==0,G[rsmall]==0},G,{r,rsmall,rlarge}];
theY = NDSolve[{Y'[r]-k/r Y[r]+(2k+1)/r (G[r]/.First[theG])==0,Y[rlarge]==(G[rlarge]/.First[theG])},Y,{r,rlarge,rsmall}];
Return[Y/.First[theY]];
]

AmbitRadialIntegral[k_,{fa_,ga_},{fb_,gb_},{fc_,gc_},{fd_,gd_}]:=Module[{rsmall,rlarge,Y},
{{rsmall,rlarge}}=fa[[1]];
Y=AmbitHartreeY[k,{fb,gb},{fd,gd}];
NIntegrate[Y[r](fa[r] fc[r] + ga[r] gc[r]),{r,rsmall,rlarge}]
]


End[];
EndPackage[];
