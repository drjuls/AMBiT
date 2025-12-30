#ifndef INCLUDE_H
#define INCLUDE_H

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <new>
#include <sstream>
#include <fstream>
#include <numeric>
#include <memory>
#include <gsl/gsl_math.h>

#include "Atom/Debug.h"
#include "Atom/OutStreams.h"

#define string_macro(x) _stringify(x)
#define _stringify(x) #x

#ifdef WIN32
    #pragma warning ( disable : 4786 )
    #pragma warning ( disable : 4290 )
    #pragma warning ( disable : 4503 )
    #pragma warning ( disable : 4267 )  // conversion from * to *, possible loss of data

    #define WIN32_LEAN_AND_MEAN     // Exclude rarely-used stuff from Windows headers
    #define mmin std::min
    #define mmax std::max

    #ifdef _DEBUG
        #define PAUSE getchar();
    #endif

#else
    #ifndef UNIX
        #define UNIX
        // Allows us to set file permissions with umask. Only works with UNIX, though
        #include <sys/stat.h>
    #endif

    #ifdef GCC
        #define mmin(i, j) (i <? j)
        #define mmax(i, j) (i >? j)
    #else
        #define mmin(i, j) (((i) < (j))?(i):(j))
        #define mmax(i, j) (((i) > (j))?(i):(j))
    #endif
#endif

#ifndef PAUSE
    #define PAUSE
#endif

/** Default for appending underscore to fortran subroutines. */
#ifndef _FUS
    #define _FUS 1
#endif

namespace Ambit
{
/** I/O streams */
extern std::ostream* outstream;
extern std::ostream* logstream;
extern std::ostream* errstream;
extern FileErrHandler* file_err_handler; 

extern Debug DebugOptions;

/** MPI */
extern int NumProcessors;
extern int ProcessorRank;

extern std::string itoa(int value, unsigned int base = 10);

/** Change endianness when reading. */
//#define fread fread_other_endian
extern size_t fread_other_endian(void* ptr, size_t size, size_t count, FILE * fp);

}
#endif
