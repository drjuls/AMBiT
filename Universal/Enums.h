#ifndef ENUMS_H
#define ENUMS_H

#include <string>

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
    
    inline std::string LowerName(Enum from)
    {
        switch(from)
        {
            case Even:
                return "even";
            case Odd:
                return "odd";
            default:
                return "";
        }
    }

    inline std::string ShortName(Enum from)
    {
        switch(from)
        {
            case Even:
                return "+1";
            case Odd:
                return "-1";
            default:
                return "";
        }
    }
    
    inline std::string VeryShortName(Enum from)
    {
        switch(from)
        {
            case Even:
                return "e";
            case Odd:
                return "o";
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
