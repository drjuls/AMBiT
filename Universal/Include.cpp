#include "Include.h"

// Use the standard fread inside our own fread.
#undef fread

namespace Ambit
{
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

std::string itoa(int value, unsigned int base)
{
	const char digitMap[] = "0123456789abcdef";	
	std::string buf;

	// Guard:
	if (base == 0 || base > 16)
    {	// Error: may add more trace/log output here
		return buf;
	}

	// Take care of negative int:
	std::string sign;	
	int _value = value;

	// Check for case when input is zero:	
	if (_value == 0)
        return "0";

	if (value < 0)
    {	_value = -value;
		sign = "-";
	}

	// Translating number to string with base (max 32 digits):
	for(int i = 32; _value && i; --i)
    {	buf = digitMap[ _value % base ] + buf;
		_value /= base;
	}

	return sign.append(buf);	
}
}
