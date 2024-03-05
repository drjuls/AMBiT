#ifndef ENUMS_H
#define ENUMS_H

#include <string>

namespace Ambit
{
enum class Parity { even, odd };

enum class ContinuumNormalisation { LandauNu, LandauEnergy, Cowan, Unitary };

enum class MultipolarityType { E, M };

enum class TransitionGauge { Length, Velocity };

inline int Sign(const Parity& p)
{
    if(p == Parity::even)
        return 1;
    else
        return -1;
}

inline Parity operator*(const Parity& first, const Parity& second)
{
    if(first == second)
        return Parity::even;
    else
        return Parity::odd;
}

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

inline std::string LetterName(Parity from)
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

inline std::string ShortName(Parity from)
{
    switch(from)
    {
        case Parity::even:
            return "+";
        case Parity::odd:
            return "-";
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

inline std::string Name(MultipolarityType from)
{
    switch(from)
    {
        case MultipolarityType::E:
            return "E";
        case MultipolarityType::M:
            return "M";
        default:
            return "";
    }
}

inline std::string Name(TransitionGauge from)
{
    switch(from)
    {
        case TransitionGauge::Length:
            return "Length";
        case TransitionGauge::Velocity:
            return "Velocity";
        default:
            return "";
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

}
#endif
