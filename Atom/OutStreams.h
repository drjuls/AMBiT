#ifndef OUT_STREAMS_H
#define OUT_STREAMS_H

#include <iostream>

class OutStreams
{
public:
    OutStreams() {}
    ~OutStreams() {}

    static void InitialiseStreams();
    static void FinaliseStreams();
};


#endif
