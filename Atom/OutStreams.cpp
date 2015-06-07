#include "OutStreams.h"
#include "Include.h"
#include <fstream>

std::ostream* outstream;
std::ostream* logstream;
std::ostream* errstream;

void OutStreams::InitialiseStreams()
{
#ifdef WIN32
    outstream = &std::cout;
    errstream = &std::cerr;
    logstream = &std::cerr;
#endif

#ifdef UNIX   // Don't output cerr to screen
    char* jobid = getenv("PBS_JOBID");

    std::string errorname = "error";
    std::string logname = "log";
    std::string suffix = ".out";

    if(jobid)
    {   logname += jobid;
        errorname += jobid;
    }

    if(NumProcessors == 1)
    {   outstream = &std::cout;
        errstream = new std::ofstream(errorname + suffix);
        logstream = new std::ofstream(logname + suffix);
    }
    else
    {   errstream = new std::ofstream(errorname + '_' + itoa(ProcessorRank) + suffix);
        logstream = new std::ofstream(logname + '_' + itoa(ProcessorRank) + suffix);

        if(ProcessorRank == 0)
            outstream = &std::cout;
        else
            outstream = logstream;
    }
#endif
}

void OutStreams::FinaliseStreams()
{
#ifdef UNIX
    dynamic_cast<std::ofstream*>(errstream)->close();
    dynamic_cast<std::ofstream*>(logstream)->close();
#endif
}
