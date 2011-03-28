#ifndef ENUMS_H
#define ENUMS_H

#include <string>

namespace DisplayOutputType
{
    enum Enum 
    {
        Start,
        Standard = 0,
        Short,
        End
    };
    
    inline std::string Name(Enum from)
    {
        switch(from)
        {
            case Standard:
                return "Standard";
            case Short:
                return "Short";
            default:
                return "";
        }
    }
}

namespace MultipolarityType
{
    enum Enum 
    {
        Start,
        E = 0,
        M,
        End
    };

    inline std::string Name(Enum from)
    {
        switch(from)
        {
            case E:
                return "E";
            case M:
                return "M";
            default:
                return "";
        }
    }
}

namespace TransitionGaugeType
{
    enum Enum 
    {
        Start,
        Length = 0,
        Velocity,
        End
    };

    inline std::string Name(Enum from)
    {
        switch(from)
        {
            case Length:
                return "Length";
            case Velocity:
                return "Velocity";
            default:
                return "";
        }
    }
}

namespace ParityType
{
    enum Enum 
    {
        Start,
        Even = 0,
        Odd,
        End
    };

    inline std::string Name(Enum from)
    {
        switch(from)
        {
            case Even:
                return "Even";
            case Odd:
                return "Odd";
            default:
                return "";
        }
    }
}

#endif
