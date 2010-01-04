#include "InputParser.h"
#include "Include.h"

bool InputParser::Parse()
{
    inputstrings.clear();
    char buffer[200];
    char key[200], value[200];
    unsigned int i, offset;
    bool atend = false;
    
    std::string skey, svalue;
    unsigned int linenum = 0;

    while(!atend && fgets(buffer, 200, fp))
    {
        linenum++;

        // Get end of buffer
        for(i = 0; i < 200; i++)
            if(iscntrl(buffer[i]) || buffer[i] == '!' ||
               (buffer[i] == '/' && buffer[i+1] == '/'))
            {   buffer[i] = 0;
                break;
            }

        // Get key
        i = 0;
        while(buffer[i] && isspace(buffer[i]))
            i++;
        offset = i;
        while(buffer[i] && !isspace(buffer[i]) && buffer[i] != '=')
        {   key[i-offset] = toupper(buffer[i]);
            i++;
        }
        key[i-offset] = 0;
        skey = key;
        if((skey == "END") || (skey == "&END") || (skey == "$END") || (skey == "EOF"))
            atend = true;

        if(!skey.empty() && !atend)
        {
            // Pass by equals sign
            while(buffer[i] && isspace(buffer[i]))
                i++;
            if(buffer[i++] != '=')
            {   *errstream << "InputParser::Parse(): bad input line " << linenum
                           << ", key = " << skey << std::endl;
                return false;
            }
            while(buffer[i] && isspace(buffer[i]))
                i++;

            // Get value
            offset = i;
            while(value[i-offset] = buffer[i])
                i++;
            // Strip trailing whitespace
            i = i - offset;
            while(isspace(value[--i]))
                value[i] = 0;
            svalue = value;
            
            if(svalue.empty())
            {   *errstream << "InputParser::Parse(): bad input line " << linenum
                           << ", key = " << skey << std::endl;
                return false;
            }

            inputstrings.insert(std::pair<std::string,std::string>(skey, svalue));
        }
    }

    return true;
}

bool InputParser::FindKey(const std::string& key, std::string& value) const
{
    std::map<std::string,std::string>::const_iterator it;
    std::string upperkey = key;
    unsigned int i;
    for(i = 0; i < key.length(); i++)
        upperkey[i] = toupper(upperkey[i]);

    it = inputstrings.find(upperkey);
    if(it == inputstrings.end())
        return false;

    value = it->second;
    return true;
}

bool InputParser::StripDelimiters(std::string& value) const
{
    int first = 0;
    int last  = value.size() - 1;

    // Strip whitespace
    while((first < last) && isspace(value[first]))
        first++;
    while((last > first) && isspace(value[last]))
        last--;

    if(((value[first] == '"') && (value[last] == '"')) ||
       ((value[first] == '{') && (value[last] == '}')) ||
       ((value[first] == '[') && (value[last] == ']')) ||
       ((value[first] == '(') && (value[last] == ')')))
    {   first++;
        last--;
    }
    
    if(first > last)
        return false;
    
    value = value.substr(first, last-first+1);
    return true;
}

template<>
bool InputParser::GetSingleValue<std::string>(const std::string& key, std::string& data) const
{
    std::string value;
    if(!FindKey(key, value) || !StripDelimiters(value))
        return false;

    data = value;
    return true;
}

template<>
bool InputParser::GetSingleValue<bool>(const std::string& key, bool& data) const
{
    std::string value;
    if(!FindKey(key, value) || !StripDelimiters(value))
        return false;

    std::string upperval = value;
    unsigned int i;
    for(i = 0; i < value.length(); i++)
        upperval[i] = toupper(upperval[i]);

    if((upperval == "TRUE") ||
       (upperval == "T"))
    {   data = true;
        return true;
    }
    else if((upperval == "FALSE") ||
            (upperval == "F"))
    {   data = false;
        return true;
    }

    return false;
}

template<>
unsigned int InputParser::GetMultipleValues<std::string>(const std::string& key, std::string* data, unsigned int size) const
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

                data[num_found] = value;
                num_found++;
                last_comma = i;
            }

            i++;
        }
    }

    return num_found;
}
