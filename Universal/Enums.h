#ifndef ENUMS_H
#define ENUMS_H

#include <string>

enum class Parity { start, even = 0, odd, end };

inline std::string Name(Parity from)
{
    switch(from)
    {
        case Parity::even:
            return "Even";
        case Parity::odd:
            return "Odd";
        default:
            return "";
    }
}

inline std::string LowerName(Parity from)
{
    switch(from)
    {
        case Parity::even:
            return "even";
        case Parity::odd:
            return "odd";
        default:
            return "";
    }
}

inline std::string ShortName(Parity from)
{
    switch(from)
    {
        case Parity::even:
            return "+1";
        case Parity::odd:
            return "-1";
        default:
            return "";
    }
}

inline std::string VeryShortName(Parity from)
{
    switch(from)
    {
        case Parity::even:
            return "e";
        case Parity::odd:
            return "o";
        default:
            return "";
    }
}

namespace DisplayOutputType
{
    enum Enum
    {
        Start,
        Standard = 0,
        CommaSeparated,
        SpaceSeparated,
        TabSeparated,
        End
    };

    inline std::string Name(Enum from)
    {
        switch(from)
        {
            case Standard:
                return "Standard";
            case CommaSeparated:
                return "Comma Separated";
            case SpaceSeparated:
                return "Space Separated";
            case TabSeparated:
                return "Tab Separated";
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

namespace TransitionCalculationMethod
{
    enum Enum 
    {
        Start,
        WalterJohnson = 0,
        RobertCowan,
        ZenonasRudzikas,
        End
    };
    
    inline std::string Name(Enum from)
    {
        switch(from)
        {
            case WalterJohnson:
                return "Walter R. Johnson";
            case RobertCowan:
                return "Robert D. Cowan";
            case ZenonasRudzikas:
                return "ZenonasRudzikas";
            default:
                return "";
        }
    }
}

#endif
