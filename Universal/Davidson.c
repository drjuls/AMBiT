/* acpz.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "blaswrap.h"

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b18 = -1.;
static doublereal c_b38 = 0.;
static real c_b51 = (float)-1.;
static doublereal c_b81 = 1.;

/* ACPZDVDSON.  A DAVIDSON PROGRAM FOR FINDING A FEW SELECTED EXTREME */
/* 1   EIGENPAIRS OF A LARGE, SPARSE, REAL, SYMMETRIC MATRIX. */
/* 2   A. STATHOPOULOS, C.F. FISCHER. */
/* REF. IN COMP. PHYS. COMMUN. 79 (1994) 268 */
/* ======================================================================= */
/* Subroutine */ int dvdson_(op, n, lim, diag, ilow, ihigh, iselec, niv, 
	mblock, crite, critc, critr, ortho, maxiter, work, iwrsz, iwork, 
	iiwsz, hiend, nloops, nmv, ierr)
/* Subroutine */ int (*op) ();
integer *n, *lim;
doublereal *diag;
integer *ilow, *ihigh, *iselec, *niv, *mblock;
doublereal *crite, *critc, *critr, *ortho;
integer *maxiter;
doublereal *work;
integer *iwrsz, *iwork, *iiwsz;
logical *hiend;
integer *nloops, *nmv, *ierr;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer neig, iicv, nume, i__;
    extern /* Subroutine */ int dscal_();
    static integer isvec;
    extern /* Subroutine */ int dcopy_(), setup_();
    static integer iscra1, iscra2, iscra3, is, ibasis, itemps, istart;
    extern /* Subroutine */ int dvdrvr_();
    static integer iab, ieigval, ioldval;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static cilist io___7 = { 0, 6, 0, 0, 0 };


/* ======================================================================= */

/*       Author: Andreas Stathopoulos, Charlotte F. Fischer */

/*       Computer Science Department */
/*       Vanderbilt University */
/*       Nashville, TN 37212 */
/*       andreas@vuse.vanderbilt.edu */
/*       cff@vuse.vanderbilt.edu                       DECEMBER 1993 */

/*       Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer */

/*       DVDSON is a Fortran77 program that finds a few selected */
/*       eigenvalues and their eigenvectors at either end of spectrum of */
/*       a large, symmetric (and usually sparse) matrix, denoted as A. */
/*       The matrix A is only referenced indirectly through the user */
/*       supplied routine OP which implements a block matrix-vector */
/*       operation(see below). Either the range of the eigenvalues wanted */
/*       or an array of the indices of selected ones can be specified. */
/*       DVDSON is a front-end routine for setting up arrays, and initial */
/*       guess (calling SETUP). It also performs detailed error checking. */
/*       DVDRVR is the driver routine that implements a version of the */
/*       Davidson algorithm. The characteristics of this version are: */
/*        o  All arrays used by the program are stored in MEMORY. */
/*        o  BLOCK method (many vectors may be targeted per iteration.) */
/*        o  Eigenvectors are targeted in an optimum way without */
/*           the need to compute all unconverged residuals, */
/*        o  It REORTHOGONILIZES the basis in case of orthogonality loss. */
/*        o  Finds HIGHEST eigenpairs by using the negative of the A. */
/*        o  Finds SELECTED eigenpairs specified by the user. */
/*        o  It accepts INITIAL eigenvector ESTIMATES or it can */
/*           CREATE INITIAL ESTIMATES from the diagonal elements. */
/*        o  It uses a USER SUPPLIED block matrix-vector operation, OP. */
/*           Depending on the implementation, OP can operate in either */
/*           memory or on disc, and for either sparse or dense matrix. */
/*        o  The user can provide STOPPING CRITERIA for eigenvalues, */
/*           and residuals. The user can also CONTROL reorthogonalization */
/*            and block size. */
/*        o  On exit INFORMATION is given about the convergence status */
/*           of eigenpairs and the number of loops and OP operations. */

/*       The program consists of the following routines: */
/*       DVDSON, SETUP, DVDRVR, ADDABS, TSTSEL, */
/*       MULTBC, OVFLOW, NEWVEC, ORTHNRM. */
/*       It also calls some basic BLAS routines: */
/*       DCOPY, DSCAL, DDOT, DAXPY, IDAMAX, DGEMV, DINIT */
/*       For solving the small eigenproblem, the routine DSPEVX from */
/*       LAPACK is used. DSPEVX is obtainable from NETLIB, together */
/*       with a series of subroutines that it calls. */

/*       All the routines have IMPLICIT DOUBLE PRECISION(A-H,O-Z) */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*  (Important to the following is the concept of NUME, the distance of */
/*   the index of the eigenpair wanted which is farthest from the */
/*   extremes,i.e., */
/*      if  lowest  eigepairs i1<i2<...<ik are wanted, NUME=ik */
/*      if highest eigenpairs i1<i2<...<ik are wanted, NUME=N-i1+1 */
/*   where i1,...,ik are the indices of the wanted eigenpairs. */
/*   Obviously, NUME.GE.(No. of EiGenpairs wanted). ) */
/*   on entry */
/*   ------- */
/*   OP          User supplied routine with calling sequence OP(N,M,B,C). */
/*               B and C are N x M matrices and C stores the result AxB. */
/*               It should be declared external in the main program. */
/*   N           Order of the matrix. */
/*   LIM         The upper limit on the dimension of the expanding basis. */
/*               NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed */
/*               only for LIM=NUME=N. The choice of LIM depends on the */
/*               available workspace (see below). If the space is */
/*               available it is preferable to have a large LIM, but not */
/*               larger than NUME$+$40. */
/*   DIAG        Array of size N with the diagonal elements of the */
/*               matrix A. */
/*   ILOW        The index of the lowest eigepair to be computed. If */
/*               (ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs */
/*               to be computed should be contained in array ISELEC. */
/*               (Modified on exit). */
/*   IHIGH       The index of the highest eigenpair to be computed. */
/*               Considered ONLY when ILOW is in the range */
/*               (0.LT.ILOW.LE.N). (Modified on exit). */
/*   ISELEC      Array of size LIM holding the user specified indices */
/*               for the eigenpairs to be computed. Considered only when */
/*               (ILOW.LE.0).or.(ILOW.GT.N). The indices are read from */
/*               the first position until a non positive integer is met. */
/*                  Example: if N=500, ILOW=0, and ISELEC(1)=495, */
/*                  ISELEC(2)=497, ISELEC(3)=-1, the program will find */
/*                  2 of the highest eigenpairs, pairs 495 and 497. */
/*               Any order of indices is acceptable (Modified on exit). */
/*   NIV         Number of Initial Vector estimates provided by the user. */
/*               If NIV is in the range:  (NUME).LE.(NIV).LE.(LIM), */
/*               the first NIV columns of size N of WORK should contain */
/*               the estimates (see below). In all other cases of NIV, */
/*               the program generates initial estimates. */
/*   MBLOCK      Number of vectors to be targeted in each iteration. */
/*               1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold. */
/*               Large block size reduces the number of iterations */
/*               (matrix acceses) but increases the matrix-vector */
/*               multiplies. It should be used when the matrix accese */
/*               is expensive (disc, recomputed or distributed). */
/*   CRITE       Convergence threshold for eigenvalues. */
/*               If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted */
/*               eigenvalues, convergence is signaled. */
/*   CRITC       Convergence threshold for the coefficients of the last */
/*               added basis vector(s). If all of those corresponding to */
/*               unconverged eigenpairs are less than CRITC convergence */
/*               is signaled. */
/*   CRITR       Convergence threshold for residual vector norms. If */
/*               all the residual norms ||Ax_i-l_ix_i|| of the targeted */
/*               x_i are less than CRITR convergence is signaled. */
/*               If ANY of the criteria are satisfied the algorithm stops */
/*   ORTHO       The threshold over which loss of orthogonality is */
/*               assumed. Usually ORTHO.LE.CRITR*10 but the process can */
/*               be skipped by setting ORTHO to a large number(eg,1.D+3). */
/*   MAXITER     Upper bound on the number of iterations of the */
/*               algorithm. When MAXITER is exceeded the algorithm stops. */
/*               A typical MAXITER can be MAX(200,NUME*40), but it can */
/*               be increased as needed. */
/*   WORK        Real array of size IWRSZ. Used for both input and output */
/*               If NIV is in ((NUME).LE.(NIV).LE.(LIM)), on input, WORK */
/*               must have the NIV initial estimates. These NIV N-element */
/*               vectors start from WORK(1) and continue one after the */
/*               other. They must form an orthonormal basis. */
/*   IWRSZ       The size of the real workspace. It must be at least as */
/*               large as: */

/*                       2*N*LIM + LIM*LIM + (NUME+10)*LIM + NUME */

/*   IWORK       Integer work array of size IIWSZ. Used as scrath array */
/*               for indices and for use in the LAPACK routines. */
/*   IIWSZ       The size of the integer workspace. It must be at least */
/*               as large as: */
/*                                    6*LIM + NUME */

/*               If LIM or NUME needs to be increased, the space should */
/*               also be increased accordingly. For given IWRSZ and */
/*               IIWSZ one can calculate how big a problem one can */
/*               solve (LIM,NUME). */

/*   on exit */
/*   ------- */
/*   WORK(1)     The first NUME*N locations contain the approximations to */
/*               the NUME extreme eigenvectors. If the lowest eigenpairs */
/*               are required, (HIEND=false), eigenvectors appear in */
/*               ascending order, otherwise (HIEND=false), they appear in */
/*               descending order. If only some are requested, the order */
/*               is the above one for all the NUME extreme eigenvectors, */
/*               but convergence has been reached only for the selected */
/*               ones. The rest are the current approximations to the */
/*               non-selected eigenvectors. */
/*   WORK(NUME*N+1) */
/*               The next NUME locations contain the approximations to */
/*               the NUME extreme eigenvalues, corresponding to the above */
/*               NUME eigenvectors. The same ordering and convergence */
/*               status applies here as well. */
/*   WORK(NUME*N+NUME+1) */
/*               The next NUME locations contain the corresponding values */
/*               of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of */
/*               the last step of the algorithm. */
/*   WORK(NUME*N+NUME+NUME+1) */
/*               The next NUME locations contain the corresponding */
/*               residual norms of the NUME above eigenvectors, of the */
/*               last step. */
/*   HIEND       Logical. If .true. on exit the highest eigenpairs are */
/*               found in descending order. Otherwise, the lowest */
/*               eigenpairs are arranged in ascending order. */
/*   NLOOPS      The number of iterations it took to reach convergence. */
/*               This is also the number of matrix references. */
/*   NMV         The number of Matrix-vector(M-V) multiplies. Each matrix */
/*               reference can have up to size(block) M-V multiplies. */
/*   IERR        An integer denoting the completions status: */
/*               IERR = 0        denotes normal completion. */
/*               IERR = -k       denotes error in DSPEVX (k eigenpairs */
/*                               not converged) */
/*               0<IERR<=2048    denotes some inconsistency as follows: */
/*        If (INT( MOD(IERR,  2)/1  ) N < LIM */
/*        If (INT( MOD(IERR,  4)/2  ) LIM < 1 */
/*        If (INT( MOD(IERR,  8)/4  ) ISELEC(1)<1, and no range specified */
/*        If (INT( MOD(IERR, 16)/8  ) IHIGH > N (in range or ISELEC) */
/*        If (INT( MOD(IERR, 32)/16 ) IHIGH < ILOW (Invalid range) */
/*        If (INT( MOD(IERR, 64)/32 ) NEIG >= LIM (Too many wanted) */
/*        If (INT( MOD(IERR,128)/64 ) Probable duplication in ISELEC */
/*        If (INT( MOD(IERR,256)/128) NUME >= LIM (max eigen very far) */
/*        If (INT( MOD(IERR,512)/256) MBLOCK is out of bounds */
/*        If (INT( MOD(IERR,1024)/512) IWRSZ or IIWSZ is not enough */
/*        If (INT( MOD(IERR,2048)/1024) Orthogonalization Failed */
/*        If (INT( MOD(IERR,4096)/2048) NLOOPS > MAXITER */

/*               The program will also print an informative message to */
/*               the standard output when NIV is not proper but it will */
/*               continue by picking initial estimates internally. */
/* ----------------------------------------------------------------------- */

/* Checking user input errors, and setting up the problem to solve. */

    /* Parameter adjustments */
    --diag;
    --iselec;
    --work;
    --iwork;

    /* Function Body */
    *ierr = 0;
    if (*lim > *n) {
	++(*ierr);
    }
    if (*lim <= 0) {
	*ierr += 2;
    }
    *hiend = FALSE_;
    if (*ilow <= 0 || *ilow > *n) {
/*          ..Look for user choice of eigenpairs in ISELEC */
	if (iselec[1] <= 0) {
/*             ..Nothing is given in ISELEC */
	    *ierr += 4;
	} else {
/*             ..Find number of eigenpairs wanted, and their */
/*             ..min/max indices */
	    neig = 1;
	    *ilow = iselec[1];
	    *ihigh = iselec[1];
	    i__1 = *lim;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (iselec[i__] <= 0) {
		    goto L20;
		}
/* Computing MIN */
		i__2 = *ilow, i__3 = iselec[i__];
		*ilow = min(i__2,i__3);
/* Computing MAX */
		i__2 = *ihigh, i__3 = iselec[i__];
		*ihigh = max(i__2,i__3);
		++neig;
/* L10: */
	    }
/*             ..Check if a very large index is asked for */
L20:
	    if (*ihigh > *n) {
		*ierr += 8;
	    }
	}
    } else {
/*          ..Look for a range between ILOW and IHIGH */
/*          ..Invalid range. IHIGH>N */
	if (*ihigh > *n) {
	    *ierr += 8;
	}
	neig = *ihigh - *ilow + 1;
/*          ..Invalid range. IHIGH<ILOW */
	if (neig <= 0) {
	    *ierr += 16;
	}
	if (neig > *lim) {
/*             ..Not enough Basis space. Increase LIM or decrease NEIG */
	    *ierr += 32;
	} else {
/*             ..Fill in the ISELEC with the required indices */
	    i__1 = neig;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
		iselec[i__] = *ilow + i__ - 1;
	    }
	}
    }
    if (*ierr != 0) {
	return 0;
    }
    nume = *ihigh;
/*       ..Identify if few of the highest eigenpairs are wanted. */
    if (*ilow + *ihigh - 1 > *n) {
	*hiend = TRUE_;
	nume = *n - *ilow + 1;
/*          ..Change the problem to a minimum eipenpairs one */
/*          ..by picking the corresponding eigenpairs on the */
/*          ..opposite side of the spectrum. */
	i__1 = neig;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	    iselec[i__] = *n - iselec[i__] + 1;
	}
    }
/*       ..duplications in ISELEC */
    if (neig > nume) {
	*ierr += 64;
    }
/*       ..Not enough Basis space. Increase LIM or decrease NUME */
    if (nume > *lim || nume == *lim && nume != *n) {
	*ierr += 128;
    }
/*       ..Size of Block out of bounds */
    if (*mblock < 1 || *mblock > neig) {
	*ierr += 256;
    }
/*       ..Check for enough workspace for Dvdson */
    if (*iwrsz < *lim * ((*n << 1) + *lim + (nume + 10)) + nume || *iiwsz < *
	    lim * 6 + nume) {
	*ierr += 512;
    }
    if (*ierr != 0) {
	return 0;
    }
    if (*niv > *lim) {
/*          ..Check number of initial estimates NIV is lower than LIM. */
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, "WARNING: Too many initial estimates.?", (ftnlen)
		37);
	e_wsle();
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, "The routine will pick the appropriate number", (
		ftnlen)44);
	e_wsle();
    } else if (*niv < nume && *niv > 0) {
/*          ..check if enough initial estimates. */
/*          ..(NIV<1 => program chooses) */
	s_wsle(&io___6);
	do_lio(&c__9, &c__1, "WARNING: Not enough initial estimates", (ftnlen)
		37);
	e_wsle();
	s_wsle(&io___7);
	do_lio(&c__9, &c__1, "The routine will pick the appropriate number", (
		ftnlen)44);
	e_wsle();
    }

/* Assigning space for the real work arrays */

    ibasis = 1;
    ieigval = ibasis + *n * *lim;
    iab = ieigval + *lim;
    is = iab + *n * *lim;
    itemps = is + *lim * (*lim + 1) / 2;
    isvec = itemps + *lim * (*lim + 1) / 2;
    iscra1 = isvec + *lim * nume;
    ioldval = iscra1 + (*lim << 3);

/* Assigning space for the integer work arrays */

    iscra2 = 1;
    iscra3 = iscra2 + *lim * 5;
    iicv = iscra3 + *lim;
    if (*hiend) {
	dscal_(n, &c_b18, &diag[1], &c__1);
    }
    istart = *niv;
    setup_(op, n, lim, &nume, hiend, &diag[1], &iwork[1], &work[ibasis], &
	    work[iab], &work[is], &istart);
    *nloops = 1;
    *nmv = istart;
    dvdrvr_(op, n, hiend, lim, mblock, &diag[1], &nume, &istart, &neig, &
	    iselec[1], crite, critc, critr, ortho, maxiter, &work[ieigval], &
	    work[ibasis], &work[iab], &work[is], &work[itemps], &work[isvec], 
	    &work[iscra1], &iwork[iscra2], &iwork[iscra3], &iwork[iicv], &
	    work[ioldval], nloops, nmv, ierr);
    if (*hiend) {
	dscal_(n, &c_b18, &diag[1], &c__1);
	dscal_(&nume, &c_b18, &work[ieigval], &c__1);
    }

/* -Copy the eigenvalues after the eigenvectors */
/* -Next, copy the difference of eigenvalues between the last two steps */
/* -Next, copy the residuals for the first NUME estimates */

    dcopy_(&nume, &work[ieigval], &c__1, &work[ibasis + *n * nume], &c__1);
    dcopy_(&nume, &work[ioldval], &c__1, &work[ibasis + (*n + 1) * nume], &
	    c__1);
    dcopy_(&nume, &work[iscra1], &c__1, &work[ibasis + (*n + 2) * nume], &
	    c__1);
/* L100: */
    return 0;
} /* dvdson_ */

/* ======================================================================= */
/* Subroutine */ int setup_(op, n, lim, nume, hiend, diag, minelem, basis, ab,
	 s, niv)
/* Subroutine */ int (*op) ();
integer *n, *lim, *nume;
logical *hiend;
doublereal *diag;
integer *minelem;
doublereal *basis, *ab, *s;
integer *niv;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer imin, i__, j;
    extern /* Subroutine */ int dinit_();
    static integer kpass;
    extern /* Subroutine */ int addabs_();

/* ======================================================================= */
/*       Subroutine for setting up (i) the initial BASIS if not provided, */
/*       (ii) the product of the matrix A with the Basis into matrix AB, */
/*       and (iii) the small matrix S=B^TAB. If no initial estimates are */
/*       available, the BASIS =(e_i1,e_i2,...,e_iNUME), where i1,i2,..., */
/*       iNUME are the indices of the NUME lowest diagonal elements, and */
/*       e_i the i-th unit vector. (ii) and (iii) are handled by ADDABS. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*   on entry */
/*   -------- */
/*   OP          The block matrix-vector operation, passed to ADDABS */
/*   N           the order of the matrix A */
/*   LIM         The limit on the size of the expanding Basis */
/*   NUME        Largest index of the wanted eigenvalues. */
/*   HIEND       Logical. True only if the highest eigenpairs are needed. */
/*   DIAG        Array of size N with the diagonal elements of A */
/*   MINELEM     Array keeping the indices of the NUME lowest diagonals. */

/*   on exit */
/*   ------- */
/*   BASIS       The starting basis. */
/*   AB, S       The starting D=AB, and small matrix S=B^TAB */
/*   NIV         The starting dimension of BASIS. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --diag;
    --s;
    --ab;
    --basis;
    --minelem;

    /* Function Body */
    if (*niv > *lim || *niv < *nume) {

/*          ..Initial estimates are not available. Give as estimates unit */
/*          ..vectors corresponding to the NUME minimum diagonal elements */
/*          ..First find the indices of these NUME elements (in MINELEM). */
/*          ..Array AB is used temporarily as a scratch array. */

	dinit_(n, &c_b18, &ab[1], &c__1);
	i__1 = *nume;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*             ..imin= the first not gotten elem( NUME<=N ) */
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
/* L20: */
		if (ab[j] < 0.) {
		    goto L30;
		}
	    }
L30:
	    imin = j;
	    i__2 = *n;
	    for (j = imin + 1; j <= i__2; ++j) {
/* L40: */
		if (ab[j] < 0. && diag[j] < diag[imin]) {
		    imin = j;
		}
	    }
	    minelem[i__] = imin;
	    ab[imin] = 1.;
/* L10: */
	}

/*          ..Build the Basis. B_i=e_(MINELEM(i)) */

	i__1 = *n * *lim;
	dinit_(&i__1, &c_b38, &basis[1], &c__1);
	i__1 = *nume;
	for (j = 1; j <= i__1; ++j) {
	    i__ = (j - 1) * *n + minelem[j];
	    basis[i__] = 1.;
/* L50: */
	}
	*niv = *nume;
    }

/* Find the matrix AB by matrix-vector multiplies, as well as the */
/* small matrix S = B^TAB. */

    kpass = 0;
    addabs_(op, n, lim, hiend, &kpass, niv, &basis[1], &ab[1], &s[1]);
    return 0;
} /* setup_ */

/* ======================================================================= */
/* Subroutine */ int dvdrvr_(op, n, hiend, lim, mblock, diag, nume, niv, neig,
	 iselec, crite, critc, critr, ortho, maxiter, eigval, basis, ab, s, 
	temps, svec, scra1, iscra2, incv, icv, oldval, nloops, nmv, ierr)
/* Subroutine */ int (*op) ();
integer *n;
logical *hiend;
integer *lim, *mblock;
doublereal *diag;
integer *nume, *niv, *neig, *iselec;
doublereal *crite, *critc, *critr, *ortho;
integer *maxiter;
doublereal *eigval, *basis, *ab, *s, *temps, *svec, *scra1;
integer *iscra2, *incv, *icv;
doublereal *oldval;
integer *nloops, *nmv, *ierr;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static logical done;
    extern doublereal ddot_();
    static integer info, nncv, i__;
    extern /* Subroutine */ int dcopy_();
    static integer kpass;
    extern /* Subroutine */ int daxpy_();
    static logical first;
    extern /* Subroutine */ int addabs_(), multbc_(), newvec_();
    static integer nfound;
    extern /* Subroutine */ int dspevx_(), ovflow_();
    extern logical tstsel_();
    static logical restart;

/* ======================================================================= */
/*       called by DVDSON */

/*       Driver routine implementing Davidson's main loop. On entry it */
/*       is given the Basis, the work matrix D=AB and the small symmetric */
/*       matrix to be solved, S=B^TAB (as found by SETUP). In each step */
/*       the small problem is solved by calling DSPEVX. */
/*       TSTSEL tests for eigenvalue convergence and selects the next */
/*       pairs to be considered for targeting (as a block). */
/*       NEWVEC computes the new vectors (block) to be added in the */
/*       expanding basis, and tests for residual convergence. */
/*       ADDABS is the critical step of matrix multiplication. The new */
/*       vectors of D are found Dnew=ABnew, and the new small problem S, */
/*       is calculated. The algorithm is repeated. */
/*       In case of a large expanding basis (KPASS=LIM) the Basis, AB, */
/*       SVEC and S are collapsed. */
/*       At the end the current eigenvector estimates are computed as */
/*       well as the residuals and eigenvalue differences. */

/*       Subroutines called: */
/*       DSPEVX, MULTBC, TSTSEL, OVFLOW, NEWVEC, ADDABS, */
/*       DCOPY, DDOT, DAXPY */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*   on entry */
/*   ------- */

/*   OP          The user specified block-matrix-vector routine */
/*   N           The order of the matrix A */
/*   HIEND       Logical. True only if the highest eigenpairs are needed. */
/*   LIM         The limit on the size of the expanding Basis */
/*   MBLOCK      Number of vectors to be targeted in each iteration. */
/*   DIAG        Array of size N with the diagonal elements of A */
/*   NUME        The largest index of the eigenvalues wanted. */
/*   NIV         Starting dimension of expanding basis. */
/*   NEIG        Number of eigenvalues wanted. */
/*   ISELEC      Array containg the indices of those NEIG eigenpairs. */
/*   CRITE       Convergence thresholds for eigenvalues, coefficients */
/*   CRITC,CRITR and residuals. */
/*   BASIS       Array with the basis vectors. */
/*   AB          Array with the vectors D=AB */
/*   S           Array keeping the symmetric matrix of the small problem. */
/*   TEMPS       scratch array */
/*   SVEC        Array for holding the eigenvectors of S */
/*   SCRA1       Srcatch array used by DSPEVX. */
/*   ISCRA2      Integer Srcatch array used by DSPEVX. */
/*   INCV        Srcatch array used in DSPEVX. Also used in TSTSEL and */
/*               NEWVEC where it holds the Indices of uNConVerged pairs */
/*   ICV         It contains "1" to the locations of ConVerged eigenpairs */
/*   OLDVAL      Array keeping the previous' step eigenvalue estimates. */

/*   on exit */
/*   ------- */

/*   EIGVAL      Array containing the NUME lowest eigenvalues of the */
/*               the matrix A (or -A if the highest are sought). */
/*   Basis       On exit Basis stores the NUME corresponding eigenvectors */
/*   OLDVAL      On exit it stores the final differences of eigenvalues. */
/*   SCRA1       On exit it stores the NUME corresponding residuals. */
/*   NLOOPS      Number of loops taken by the algorithm */
/*   NMV         Number of matrix-vector products performed. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --diag;
    --incv;
    --iscra2;
    --scra1;
    --temps;
    --s;
    --ab;
    --basis;
    --eigval;
    --oldval;
    --icv;
    --svec;
    --iselec;

    /* Function Body */
    i__1 = *nume;
    for (i__ = 1; i__ <= i__1; ++i__) {
	eigval[i__] = 1e30;
/* L5: */
	icv[i__] = 0;
    }
    first = TRUE_;
    kpass = *niv;
    nncv = kpass;
L10:
/*       (iterations for kpass=NUME,LIM) */

/* Diagonalize the matrix S. Find only the NUME smallest eigenpairs */

    dcopy_(nume, &eigval[1], &c__1, &oldval[1], &c__1);
    i__1 = kpass * (kpass + 1) / 2;
    dcopy_(&i__1, &s[1], &c__1, &temps[1], &c__1);
    dspevx_("Vectors also", "In a range", "Upper triangular", &kpass, &temps[
	    1], &c_b51, &c_b51, &c__1, nume, &c_b38, &nfound, &eigval[1], &
	    svec[1], &kpass, &scra1[1], &iscra2[1], &incv[1], &info, (ftnlen)
	    12, (ftnlen)10, (ftnlen)16);
    *ierr = -abs(info);
    if (*ierr != 0) {
	goto L60;
    }

/* TeST for convergence on the absolute difference of eigenvalues between */
/* successive steps. Also SELect the unconverged eigenpairs and sort them */
/* by the largest magnitude in the last added NNCV rows of Svec. */

    done = tstsel_(&kpass, nume, neig, &iselec[1], &svec[1], &eigval[1], &icv[
	    1], crite, critc, &scra1[1], &iscra2[1], &oldval[1], &nncv, &incv[
	    1]);
    if (done || kpass >= *n) {
	goto L30;
    }
    if (kpass == *lim) {
/* Maximum size for expanding basis. Collapse basis, D, and S, Svec */
/* Consider the basis vectors found in TSTSEL for the newvec. */

	multbc_(n, lim, nume, &svec[1], &scra1[1], &basis[1]);
	multbc_(n, lim, nume, &svec[1], &scra1[1], &ab[1]);
	ovflow_(nume, lim, &s[1], &svec[1], &eigval[1]);
	kpass = *nume;
    }

/* Compute and add the new vectors. NNCV is set to the number of new */
/* vectors that have not converged. If none, DONE=true, exit. */

    newvec_(n, nume, lim, mblock, &kpass, critr, ortho, &nncv, &incv[1], &
	    diag[1], &svec[1], &eigval[1], &ab[1], &basis[1], &icv[1], &
	    restart, &done);
/*          ..An infinite loop is avoided since after a collapsing Svec=I */
/*          ..=> Res=Di-lBi which is just computed and it is orthogonal. */
/*          ..The following is to prevent an improbable infinite loop. */
    if (! restart) {
	first = TRUE_;
    } else if (first) {
	first = FALSE_;
	i__1 = kpass + nncv;
	multbc_(n, &i__1, nume, &svec[1], &scra1[1], &basis[1]);
	i__1 = kpass + nncv;
	multbc_(n, &i__1, nume, &svec[1], &scra1[1], &ab[1]);
	i__1 = kpass + nncv;
	ovflow_(nume, &i__1, &s[1], &svec[1], &eigval[1]);
	kpass = *nume;
	goto L10;
    } else {
	*ierr += 1024;
	goto L30;
    }
    if (done) {
	goto L30;
    }

/* Add new columns in D and S, from the NNCV new vectors. */

    addabs_(op, n, lim, hiend, &kpass, &nncv, &basis[1], &ab[1], &s[1]);
    *nmv += nncv;
    kpass += nncv;
    ++(*nloops);
    if (*nloops <= *maxiter) {
	goto L10;
    }
    *ierr += 2048;
    --(*nloops);
    kpass -= nncv;
L30:

/* Calculate final results. EIGVAL contains the eigenvalues, BASIS the */
/* eigenvectors, OLDVAL the eigenvalue differences, and SCRA1 residuals. */

    i__1 = *nume;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	oldval[i__] = (d__1 = oldval[i__] - eigval[i__], abs(d__1));
    }
    multbc_(n, &kpass, nume, &svec[1], &scra1[1], &basis[1]);
    multbc_(n, &kpass, nume, &svec[1], &scra1[1], &ab[1]);

/* i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi */
/* temporarily stored in AB(NUME*N+1) */

    i__1 = *nume;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dcopy_(n, &ab[(i__ - 1) * *n + 1], &c__1, &ab[*nume * *n + 1], &c__1);
	d__1 = -eigval[i__];
	daxpy_(n, &d__1, &basis[(i__ - 1) * *n + 1], &c__1, &ab[*nume * *n + 
		1], &c__1);
	scra1[i__] = ddot_(n, &ab[*nume * *n + 1], &c__1, &ab[*nume * *n + 1],
		 &c__1);
	scra1[i__] = sqrt(scra1[i__]);
/* L50: */
    }
L60:
    return 0;
} /* dvdrvr_ */

/* ======================================================================= */
/* Subroutine */ int addabs_(op, n, lim, hiend, kpass, nncv, basis, ab, s)
/* Subroutine */ int (*op) ();
integer *n, *lim;
logical *hiend;
integer *kpass, *nncv;
doublereal *basis, *ab, *s;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern doublereal ddot_();
    extern /* Subroutine */ int dscal_();
    static integer iv;
    static doublereal ss;
    static integer ibv, ibstart, idstart, isstart;

/* ======================================================================= */
/*       Called by: DVDRVR, SETUP */

/*       Calculates the new column in the D matrix and the new column */
/*       in the S matrix. The new D column is D(new)=AB(new). S has a */
/*       new row and column, but being symmetric only the new column is */
/*       stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i. */

/*       subroutines called: */
/*       OP, DDOT, DSCAL */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*   on entry */
/*   ------- */
/*   N           The order of the matrix A */
/*   kpass       The current dimension of the expanding sub-basis */
/*   NNCV        Number of new basis vectors. */
/*   Basis       the basis vectors, including the new NNCV ones. */
/*   on exit */
/*   ------- */
/*   AB          The new matrix D=AB. (with new NNCV columns) */
/*   S           The small matrix with NNCV new columns at the last part */
/* ----------------------------------------------------------------------- */

/* The user specified matrix-vector routine is called with the new */
/* basis vector B(*,kpass+1) and the result is assigned to AB(idstart) */

    /* Parameter adjustments */
    --s;
    --ab;
    --basis;

    /* Function Body */
    idstart = *kpass * *n + 1;
    (*op)(n, nncv, &basis[idstart], &ab[idstart]);

/* If highest pairs are sought, use the negative of the matrix */

    if (*hiend) {
	i__1 = *n * *nncv;
	dscal_(&i__1, &c_b18, &ab[idstart], &c__1);
    }

/* The new S is calculated by adding the new last columns */
/* S(new)=B^T D(new). */

    isstart = *kpass * (*kpass + 1) / 2;
    i__1 = *nncv;
    for (iv = 1; iv <= i__1; ++iv) {
	ibstart = 1;
	i__2 = *kpass + iv;
	for (ibv = 1; ibv <= i__2; ++ibv) {
	    ss = ddot_(n, &basis[ibstart], &c__1, &ab[idstart], &c__1);
	    s[isstart + ibv] = ss;
	    ibstart += *n;
/* L10: */
	}
	isstart = isstart + *kpass + iv;
	idstart += *n;
/* L20: */
    }
    return 0;
} /* addabs_ */

/* ======================================================================= */
logical tstsel_(kpass, nume, neig, iselec, svec, eigval, icv, crite, critc, 
	rowlast, ind, oldval, nncv, incv)
integer *kpass, *nume, *neig, *iselec;
doublereal *svec, *eigval;
integer *icv;
doublereal *crite, *critc, *rowlast;
integer *ind;
doublereal *oldval;
integer *nncv, *incv;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;
    logical ret_val;

    /* Local variables */
    static integer nnce;
    static logical done;
    static integer ival, icnt, icur, indx;
    static doublereal temp, tmax;
    static integer i__, l, itemp;
    extern integer idamax_();

/* ======================================================================= */

/*       Called by: DVDRVR */
/*       It first checks if the wanted eigenvalues have reached */
/*       convergence and updates OLDVAL. Second, for each wanted and non */
/*       converged eigenvector, it finds the largest absolute coefficient */
/*       of the NNCV last added vectors (from SVEC) and if not coverged, */
/*       places it in ROWLAST. IND has the corresponding indices. */
/*       Third, it sorts ROWLAST in decreasing order and places the */
/*       corresponding indices in the array INCV. The index array INCV */
/*       and the number of unconverged pairs NNCV, are passed to DVDRVR. */
/*       Later in NEWVEC only the first MBLOCK of NNCV pairs will be */
/*       targeted, since if ROWLAST(i) > ROWLAST(j) */
/*       then approximately RESIDUAL(i) > RESIDUAL(j) */

/*       Subroutines called */
/*       IDAMAX */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*   on entry */
/*   ------- */
/*   KPASS       current dimension of the expanding Basis */
/*   NUME        Largest index of the wanted eigenvalues. */
/*   NEIG        number of wanted eigenvalues of original matrix */
/*   ISELEC      index array of the wanted eigenvalues. */
/*   SVEC        the eigenvectors of the small system */
/*   EIGVAL      The NUME lowest eigenvalues of the small problem */
/*   ICV         Index of converged eigenpairs.ICV(i)=1 iff eigenpair i */
/*               has converged, and ICV(i)=0 if eigenpair i has not. */
/*   CRITE,CRITC Convergence thresholds for eigenvalues and coefficients */
/*   ROWLAST     scratch array, keeping the largest absolute coefficient */
/*               of the NNCV last rows of Svec. */
/*   IND         scratch array, temporary keeping the indices of Rowlast */
/*   OLDVAL      The previous iteration's eigenvalues. */

/*   on exit */
/*   ------- */
/*   NNCV         Number of non converged eigenvectors (to be targeted) */
/*   INCV         Index to these columns in decreasing order of magnitude */
/*   TSTSEL       true if convergence has been reached */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --oldval;
    --icv;
    --eigval;
    --svec;
    --incv;
    --ind;
    --rowlast;
    --iselec;

    /* Function Body */
    done = FALSE_;

/* Test all wanted eigenvalues for convergence under CRITE */

    nnce = 0;
    i__1 = *neig;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ival = iselec[i__];
/* L10: */
	if ((d__1 = oldval[ival] - eigval[ival], abs(d__1)) >= *crite) {
	    ++nnce;
	}
    }
    if (nnce == 0) {
	ret_val = TRUE_;
	return ret_val;
    }

/* Find the maximum element of the last NNCV coefficients of unconverged */
/* eigenvectors. For those unconverged coefficients, put their indices */
/* to IND and find their number NNCV */

    icnt = 0;
    i__1 = *neig;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (icv[iselec[i__]] == 0) {
/*             ..Find coefficient and test for convergence */
	    icur = *kpass * iselec[i__];
	    tmax = (d__1 = svec[icur], abs(d__1));
	    i__2 = *nncv - 1;
	    for (l = 1; l <= i__2; ++l) {
/* L20: */
/* Computing MAX */
		d__2 = tmax, d__3 = (d__1 = svec[icur - l], abs(d__1));
		tmax = max(d__2,d__3);
	    }
	    if (tmax < *critc) {
/*                ..this  coefficient converged */
		icv[iselec[i__]] = 1;
	    } else {
/*                ..Not converged. Add it to the list. */
		++icnt;
		ind[icnt] = iselec[i__];
		rowlast[icnt] = tmax;
	    }
	}
/* L30: */
    }
    *nncv = icnt;
    if (*nncv == 0) {
	done = TRUE_;
    }

/* Sort the ROWLAST elements interchanging their indices as well */

    i__1 = *nncv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nncv - i__ + 1;
	indx = idamax_(&i__2, &rowlast[i__], &c__1);
	incv[i__] = ind[indx + i__ - 1];
	temp = rowlast[indx + i__ - 1];
	rowlast[indx + i__ - 1] = rowlast[i__];
	rowlast[i__] = temp;
	itemp = ind[indx + i__ - 1];
	ind[indx + i__ - 1] = ind[i__];
	ind[i__] = itemp;
/* L40: */
    }
    ret_val = done;
    return ret_val;
} /* tstsel_ */

/* ======================================================================= */
/* Subroutine */ int multbc_(n, k, m, c__, temp, b)
integer *n, *k, *m;
doublereal *c__, *temp, *b;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer irow;
    extern /* Subroutine */ int dgemv_(), dcopy_();

/* ======================================================================= */
/*       called by: DVDRVR */

/*       Multiplies B(N,K)*C(K,M) and stores it in B(N,M) */
/*       Used for collapsing the expanding basis to current estimates, */
/*       when basis becomes too large, or for returning the results back */
/*       Subroutines called */
/*       DINIT, DGEMV, DCOPY */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --b;
    --temp;
    --c__;

    /* Function Body */
    i__1 = *n;
    for (irow = 1; irow <= i__1; ++irow) {
/*              CALL DINIT(M,0.d0,TEMP,1) */
	dgemv_("Transp", k, m, &c_b81, &c__[1], k, &b[irow], n, &c_b38, &temp[
		1], &c__1, (ftnlen)6);
	dcopy_(m, &temp[1], &c__1, &b[irow], n);
/* L10: */
    }
    return 0;
} /* multbc_ */

/* ======================================================================= */
/* Subroutine */ int ovflow_(nume, lim, s, svec, eigval)
integer *nume, *lim;
doublereal *s, *svec, *eigval;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer icur, i__;
    extern /* Subroutine */ int dinit_();
    static integer ind;

/* ======================================================================= */
/*       Called by: DVDRVR */
/*       Called when the upper limit (LIM) has been reached for the basis */
/*       expansion. The new S is computed as S'(i,j)=l(i)delta(i,j) where */
/*       l(i) eigenvalues, and delta of Kronecker, i,j=1,NUME. The new */
/*       eigenvectors of the small matrix are the unit vectors. */

/*       Subroutines called: */
/*       DCOPY, DINIT */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*   on entry */
/*   ------- */
/*   NUME        The largest index of eigenvalues wanted. */
/*   SVEC        the kpass eigenvectors of the smaller system solved */
/*   EIGVAL      the eigenvalues of this small system */
/*   on exit */
/*   ------- */
/*   S           The new small matrix to be solved. */
/* ----------------------------------------------------------------------- */

/* calculation of the new upper S=diag(l1,...,l_NUME) and */
/* its matrix Svec of eigenvectors (e1,...,e_NUME) */

    /* Parameter adjustments */
    --eigval;
    --svec;
    --s;

    /* Function Body */
    i__1 = *nume * (*nume + 1) / 2;
    dinit_(&i__1, &c_b38, &s[1], &c__1);
    i__1 = *nume * *nume;
    dinit_(&i__1, &c_b38, &svec[1], &c__1);
    ind = 0;
    icur = 0;
    i__1 = *nume;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[ind + i__] = eigval[i__];
	svec[icur + i__] = 1.;
	icur += *nume;
/* L10: */
	ind += i__;
    }
    return 0;
} /* ovflow_ */

/* ======================================================================= */
/* Subroutine */ int newvec_(n, nume, lim, mblock, kpass, critr, ortho, nncv, 
	incv, diag, svec, eigval, ab, basis, icv, restart, done)
integer *n, *nume, *lim, *mblock, *kpass;
doublereal *critr, *ortho;
integer *nncv, *incv;
doublereal *diag, *svec, *eigval, *ab, *basis;
integer *icv;
logical *restart, *done;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer icvc, icur, indx, irow;
    extern doublereal dnrm2_();
    static integer newstart, i__;
    extern /* Subroutine */ int dgemv_();
    static doublereal dg;
    static integer nadded, limadd;
    static doublereal ss;
    extern /* Subroutine */ int orthnrm_();

/* ======================================================================= */

/*       Called by: DVDRVR */

/*       It calculates the new expansion vectors of the basis. */
/*       For each one of the vectors in INCV starting with the largest */
/*       megnitude one, calculate its residual Ri= DCi-liBCi and check */
/*       the ||Ri|| for convergence. If it is converged do not add it */
/*       but look for the immediate larger coefficient and its vector. */
/*       The above procedure continues until MBLOCK vectors have been */
/*       added to the basis, or the upper limit has been encountered. */
/*       Thus only  the required MBLOCK residuals are computed. Then, */
/*       calculate the first order correction on the added residuals */
/*       Ri(j) = Ri(j)/(li-Ajj) and orthonormalizes the new vectors */
/*       to the basis and to themselves. */

/*       Subroutines called: */
/*       ORTHNRM, DDOT, DGEMV */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*   on entry */
/*   -------- */
/*   N           The order of the matrix A */
/*   NUME        The largest index of the eigenvalues wanted. */
/*   LIM         The limit on the size of the expanding Basis */
/*   MBLOCK      Maximum number of vectora to enter the basis */
/*   KPASS       the current dimension of the expanding basis */
/*   CRITR       Convergence threshold for residuals */
/*   ORTHO       Orthogonality threshold to be passed to ORTHNRM */
/*   NNCV        Number of Non ConVerged pairs (MBLOCK will be targeted) */
/*   INCV        Index to the corresponding SVEC columns of these pairs. */
/*   DIAG        Array of size N with the diagonal elements of A */
/*   SVEC,EIGVAL Arrays holding the eigenvectors and eigenvalues of S */
/*   AB          Array with the vectors D=AB */
/*   BASIS       the expanding basis having kpass vectors */
/*   ICV         Index of converged eigenpairs (ICV(i)=1 <=>i converged) */
/*   on exit */
/*   ------- */
/*   NNCV        The number of vectors finally added to the basis. */
/*   BASIS       The new basis incorporating the new vectors at the end */
/*   ICV         Index of converged eigenpairs (updated) */
/*   DONE        logical, if covergance has been reached. */
/*   RESTART     logical, if because of extreme loss of orthogonality */
/*               the Basis should be collapsed to current approximations. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --diag;
    --icv;
    --incv;
    --basis;
    --ab;
    --eigval;
    --svec;

    /* Function Body */
    *done = FALSE_;
    newstart = *kpass * *n + 1;
    nadded = 0;
    icvc = 0;
/* Computing MIN */
    i__1 = *lim, i__2 = *mblock + *kpass;
    limadd = min(i__1,i__2);
    icur = newstart;

/* Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors. */

    i__1 = *nncv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	indx = incv[i__];
/*          ..Compute  Newv=BASIS*Svec_indx , then */
/*          ..Compute  Newv=AB*Svec_indx - eigval*Newv and then */
/*          ..compute the norm of the residual of Newv */
	dgemv_("N", n, kpass, &c_b81, &basis[1], n, &svec[(indx - 1) * *kpass 
		+ 1], &c__1, &c_b38, &basis[icur], &c__1, (ftnlen)1);
	d__1 = -eigval[indx];
	dgemv_("N", n, kpass, &c_b81, &ab[1], n, &svec[(indx - 1) * *kpass + 
		1], &c__1, &d__1, &basis[icur], &c__1, (ftnlen)1);
	ss = dnrm2_(n, &basis[icur], &c__1);

/*          ..Check for convergence of this residual */

	if (ss < *critr) {
/*             ..Converged,do not add. Go for next non converged one */
	    ++icvc;
	    icv[indx] = 1;
	    if (icvc < *nncv) {
		goto L10;
	    }
/*             ..All have converged. */
	    *done = TRUE_;
	    return 0;
	} else {
/*             ..Not converged. Add it in the basis */
	    ++nadded;
	    incv[nadded] = indx;
	    if (nadded + *kpass == limadd) {
		goto L20;
	    }
/*             ..More to be added in the block */
	    icur += *n;
	}
L10:
	;
    }
L20:
    *nncv = nadded;

/* Diagonal preconditioning: newvect(i)=newvect(i)/(l-Aii) */
/* If (l-Aii) is very small (or zero) divide by 10.D-6 */

    icur = newstart - 1;
    i__1 = *nncv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (irow = 1; irow <= i__2; ++irow) {
	    dg = eigval[incv[i__]] - diag[irow];
	    if (abs(dg) > 1e-13) {
		basis[icur + irow] /= dg;
	    } else {
		basis[icur + irow] /= 1e-13;
	    }
/* L40: */
	}
	icur += *n;
/* L50: */
    }

/* ORTHONORMALIZATION */

    orthnrm_(n, lim, ortho, kpass, nncv, &ab[newstart], &basis[1], restart);
/* L99: */
    return 0;
} /* newvec_ */

/* ======================================================================= */
/* Subroutine */ int orthnrm_(n, lim, ortho, kpass, nncv, scra1, basis, 
	restart)
integer *n, *lim;
doublereal *ortho;
integer *kpass, *nncv;
doublereal *scra1, *basis;
logical *restart;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern doublereal ddot_();
    static doublereal dcur;
    static integer icur, i__;
    extern /* Subroutine */ int dscal_(), dcopy_();
    static doublereal dprev;
    extern /* Subroutine */ int daxpy_();
    static integer iv, ibstart;

/* ======================================================================= */

/*       It orthogonalizes the new NNCV basis vectors starting from the */
/*       kpass+1, to the previous vectors of the basis and to themselves. */
/*       A Gram-Schmidt method is followed after which the residuals */
/*       should be orthogonal to the BASIS. Because of machine arithmetic */
/*       errors this orthogonality may be lost, and a reorthogonalization */
/*       procedure is adopted whenever orthogonality loss is above a */
/*       ORTHO. If after some reorthogonalizations the procedure does not */
/*       converge to orthogonality, the basis is collapsed to the */
/*       current eigenvector approximations. */

/*       Subroutines called: */
/*       DAXPY, DDOT, DSCAL */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*   on entry */
/*   -------- */
/*   N           The order of the matrix A */
/*   LIM         The limit on the size of the expanding Basis */
/*   ORTHO       The orthogonality threshold */
/*   KPASS       The number of basis vectors already in Basis */
/*   NNCV        The number of new vectors in the basis */
/*   SCRA1       Scratch vector of size N */
/*   BASIS       the expanding basis having kpass vectors */

/*   on exit */
/*   ------- */
/*   BASIS       the new basis orthonormalized */
/*   RESTART     Logical, if true the algoritm will collapse BASIS. */
/* ----------------------------------------------------------------------- */

/* ORTHOGONALIZATION */

    /* Parameter adjustments */
    --scra1;
    --basis;

    /* Function Body */
    *restart = FALSE_;
    icur = *kpass * *n + 1;

/*       .. do iv=1,nncv */
    iv = 1;
L30:
    dprev = 1e7;
L5:
    dcur = 0.;
    ibstart = 1;
    i__1 = *kpass + iv - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	scra1[i__] = ddot_(n, &basis[ibstart], &c__1, &basis[icur], &c__1);
/* Computing MAX */
	d__2 = dcur, d__3 = (d__1 = scra1[i__], abs(d__1));
	dcur = max(d__2,d__3);
	ibstart += *n;
/* L10: */
    }
    ibstart = 1;
    i__1 = *kpass + iv - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = -scra1[i__];
	daxpy_(n, &d__1, &basis[ibstart], &c__1, &basis[icur], &c__1);
	ibstart += *n;
/* L20: */
    }
    if (dcur >= *ortho) {
	if (dcur > dprev) {
	    *restart = TRUE_;
/*                ..Adjust the number of added vectors. */
	    *nncv = iv - 1;
	    return 0;
	} else {
	    dprev = dcur;
	    goto L5;
	}
    }

/* NORMALIZATION */

    scra1[1] = ddot_(n, &basis[icur], &c__1, &basis[icur], &c__1);
    scra1[1] = sqrt(scra1[1]);
    if (scra1[1] < 1e-14) {
	dcopy_(n, &basis[*n * (*nncv - 1) + 1], &c__1, &basis[icur], &c__1);
	--(*nncv);
    } else {
	d__1 = 1 / scra1[1];
	dscal_(n, &d__1, &basis[icur], &c__1);
	icur += *n;
	++iv;
    }
    if (iv <= *nncv) {
	goto L30;
    }
    return 0;
} /* orthnrm_ */

/* ======================================================================= */
/* Subroutine */ int dinit_(n, a, x, incx)
integer *n;
doublereal *a, *x;
integer *incx;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, xaddr;

/* ======================================================================= */
/*       PURPOSE ... INITIALIZES DOUBLE PRECISION VECTOR TO */
/*                   A CONSTANT VALUE 'A' */
/* ======================================================================= */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*incx == 1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[i__] = *a;
/* L100: */
	}
    } else {
	xaddr = 1;
	if (*incx < 0) {
	    xaddr = (-(*n) + 1) * *incx + 1;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[xaddr] = *a;
	    xaddr += *incx;
/* L200: */
	}
    }
    return 0;
} /* dinit_ */
