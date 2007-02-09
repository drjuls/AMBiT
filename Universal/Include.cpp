#include "Include.h"

// Use the standard fread inside our own fread.
#undef fread

extern size_t fread_other_endian(void* ptr, size_t size, size_t count, FILE* fp)
{
    size_t num_read;

    unsigned char* lebuffer = new unsigned char[size * count];
    unsigned char* bebuffer = reinterpret_cast<unsigned char*>(ptr);

    num_read = fread(lebuffer, size, count, fp);

    for(size_t i = 0; i < count; i++)
    {
        size_t base = i * size;

        for(size_t j = 0; j < size; j++)
            bebuffer[base + j] = lebuffer[base + size - j - 1];
    }

    return num_read;
}
