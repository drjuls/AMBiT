#include "OutStreams.h"
#include "Include.h"
#include <fstream>
#include <cstdio>
#include <memory>
#include <cstring>
#include <cerrno>
#ifdef UNIX
#include <signal.h>
#endif

namespace Ambit
{
std::ostream* outstream;
std::ostream* logstream;
std::ostream* errstream;
FileErrHandler* file_err_handler; 

void OutStreams::InitialiseStreams()
{
#ifdef WIN32
    outstream = &std::cout;
    errstream = &std::cerr;
    logstream = &std::cerr;
#endif

#ifdef UNIX
    char* jobid = getenv("PBS_JOBID"); // TODO: This should also account for other queueing systems

    std::string logname = "log";
    std::string suffix = ".out";

    if(jobid)
        logname += jobid;
        

    if(NumProcessors == 1)
    {   outstream = &std::cout;
        logstream = new std::ofstream(logname + suffix);
    }
    else
    {   logstream = new std::ofstream(logname + '_' + itoa(ProcessorRank) + suffix);

        if(ProcessorRank == 0)
            outstream = &std::cout;
        else
            outstream = logstream;
    }
    // Send all errors to the output files, since they're more likely to be noticed this way
    errstream = outstream;

    // Also initialise the wrapper around fwrite() and friends
    file_err_handler = new FileErrHandler();
#endif
}

void OutStreams::FinaliseStreams()
{
#ifdef UNIX
    dynamic_cast<std::ofstream*>(logstream)->close();
#endif
}

FILE* FileErrHandler::fopen(const char* path, const char* mode)
{
    
    FILE* ofp = std::fopen(path, mode);
    if(!ofp && errno != ENOENT){ // Don't print a warning if the file doesn't exist - fopen() will automatically create it
        *errstream << "WARNING: error opening file - " << std::strerror(errno) << "\n";
    }
    return(ofp);
}

int FileErrHandler::fclose(FILE* stream)
{

    int ret = std::fclose(stream);
    if(ret){
        *errstream << "WARNING: error writing to file - " << std::strerror(errno) << "\n";
    }
    print_errors = true; // Reset the error counter so we can complain about future files
    return(ret);
}

int FileErrHandler::fwrite(const void* buf, size_t size, size_t nr, FILE* stream)
{
    
    size_t ret = std::fwrite(buf, size, nr, stream);
    if(!ret){
        if(print_errors){
            print_errors = false;
            *errstream << "WARNING: error writing to file - " << std::strerror(errno) << "\n";
        }
    }
    return(ret);
}

int FileErrHandler::fread(void* buf, size_t size, size_t nr, FILE* stream)
{
    size_t ret = std::fread(buf, size, nr, stream);
    
    // fread uses a return value less than nr for both errors and EOF, so we need to explicitly check
    // std::ferror to see which happened.
    // NOTE: we want to alert whenever we get an fread() error, since this can tank the entire
    // calculation.
    if(ret < nr)
    {
        if(std::ferror(stream))
        {
            *errstream << "WARNING: error reading from file - " << std::strerror(errno) << "\n";
            std::clearerr(stream);
        }
        else if(std::feof(stream))
        {
            *errstream << "Error: unexpectedly reached end of data in file\n";
            #ifdef UNIX
                // Send SIGTERM to ourselves so we can emit a stack-trace
                raise(SIGTERM);
            #else
                exit(EXIT_FAILURE);
            #endif
        }
    }
    return(ret);
}
}
