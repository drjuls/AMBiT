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
    if(NumProcessors == 1)
    {   outstream = &std::cout;
        errstream = new std::ofstream("error.out");
        logstream = new std::ofstream("log.out");
    }
    else
    {   char buffer[20];
        sprintf(buffer, "error_%d.out", ProcessorRank);
        errstream = new std::ofstream(buffer);

        sprintf(buffer, "log_%d.out", ProcessorRank);
        logstream = new std::ofstream(buffer);

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

//OutStreams& operator<<(OutStreams& out, double in)
//{
//    *outstream << "o " << in;
//    *errstream << "e " << in;
//    return out;
//}
//
//OutStreams& operator<<(OutStreams& out, std::ostream& (*_Pfn)(std::ostream&))
//{
//    *outstream << "o " << *(_Pfn);
//    *errstream << "e " << *(_Pfn);
//    return out;
//}
