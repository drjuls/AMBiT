#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <string>
#include <map>

class InputParser
{
    /** Input parsing from a file, similar to Fortran namelists, perhaps a bit more flexible.
        Input file should be of the form
            InputName = "String" or {String} or just String
            InputNum  = 22 or 2.002 or whatever
            List      = 2, 4, 6 or {2, 4, 6} or {{T W O}, "Four", SIX} or whatever.
        Comments marked with ! or //, extra whitespace generally ignored.
        End of input for reading with END or EOF.
     */
public:
    InputParser(const std::string& filename)
    {
        fp = fopen(filename.c_str(), "rt");
    }

    virtual ~InputParser()
    {
        if(fp)
            fclose(fp);
    }

    /** Get a single value from the input data. Return success. */
    template<class T>
    bool GetSingleValue(const std::string& key, T& data) const;

    /** Get a list of values for a single key from the input data.
        PRE: Memory for data[size] is already allocated.
        POST: Return number of data values actually found. 0 <= return <= size.
     */
    template<class T>
    unsigned int GetMultipleValues(const std::string& key, T* data, unsigned int size) const;

    /** Parse the input file and store key/value pairs.
        Return success; errors are generally in the input file. */
    bool Parse();

protected:
    /** Really just value = inputstrings.find(key) with some failsafes. Return success. */
    bool FindKey(const std::string& key, std::string& value) const;

    /** Strip whitespace, quotation marks, brackets, etc from value. */
    bool StripDelimiters(std::string& value) const;

protected:
    FILE* fp;
    std::map<std::string, std::string> inputstrings;

};

template<class T>
bool InputParser::GetSingleValue(const std::string& key, T& data) const
{
    std::string value;
    if(!FindKey(key, value) || !StripDelimiters(value))
        return false;

    std::stringstream sval(value);
    sval >> data;
    return true;
}

template<>
bool InputParser::GetSingleValue<bool>(const std::string& key, bool& data) const;

template<>
bool InputParser::GetSingleValue<std::string>(const std::string& key, std::string& data) const;

template<class T>
unsigned int InputParser::GetMultipleValues(const std::string& key, T* data, unsigned int size) const
{
    unsigned int num_found = 0;
    std::string valuelist;
    if(!FindKey(key, valuelist))
        return 0;

    // Is the list comma-separated?
    bool comma_sep = false;
    unsigned int i;
    for(i = 0; i < valuelist.size(); i++)
        if(valuelist[i] == ',')
        {   comma_sep = true;
            break;
        }

    if(comma_sep == true)
    {   // safe to remove delimiters
        if(!StripDelimiters(valuelist))
            return 0;

        // separate values, strip and convert
        unsigned int last_comma = 0;
        i = 0;
        while(i <= valuelist.size() && num_found < size)
        {
            if((valuelist[i] == ',') || (i == valuelist.size()))
            {   std::string value;
                if(last_comma)
                    value = valuelist.substr(last_comma+1, i - last_comma - 1);
                else
                    value = valuelist.substr(0, i);

                if(!StripDelimiters(value))
                    return 0;

                std::stringstream sdata(value);
                sdata >> data[num_found];
                num_found++;
                last_comma = i;
            }

            i++;
        }
    }

    return num_found;
}

template<>
unsigned int InputParser::GetMultipleValues<std::string>(const std::string& key, std::string* data, unsigned int size) const;

#endif
