#ifndef OUT_STREAMS_H
#define OUT_STREAMS_H

#include <iostream>
#include <cstdio>

class OutStreams
{
public:
    OutStreams() {}
    ~OutStreams() {}

    static void InitialiseStreams();
    static void FinaliseStreams();
};

/* Wrapper class around the C buffered IO functions fopen, fwrite and fclose. This class' methods take
 * the same arguments as the standard versions of these functions, but also automatically check for
 * errors and optionally raise an exception if something goes wrong.
 */
class FileErrHandler
{
public:
    FileErrHandler(): print_errors(true) {};
    ~FileErrHandler() = default;

    // Buffered IO functions which do not throw exceptions on failure (but still print error
    // messages)
    FILE* fopen(const char* path, const char* mode);
    int fclose(FILE* stream);
    int fwrite(const void* buf, size_t size, size_t nr, FILE* stream);

protected:
    bool print_errors; // Ensures we don't print too many error messages for the block of writes
};

#endif
