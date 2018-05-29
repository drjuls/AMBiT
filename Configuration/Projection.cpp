#include "Include.h"
#include "Projection.h"
#include "RelativisticConfiguration.h"
#include "HartreeFock/NonRelInfo.h"

namespace Ambit
{
Projection::Projection(const RelativisticConfiguration& relconfig, const std::vector<int>& TwoMs)
{
    config.reserve(relconfig.ParticleNumber());

    int count = 0;
    for(auto& relconfig_it: relconfig)
    {
        const OrbitalInfo& orbital = relconfig_it.first;
        for(int i = 0; i < abs(relconfig_it.second); i++)
            if(relconfig_it.second > 0)
                config.emplace_back(orbital.PQN(), orbital.Kappa(), TwoMs[count++]);
            else
                config.emplace_back(orbital.PQN(), orbital.Kappa(), -TwoMs[count++], true);
    }

    config.shrink_to_fit();
}

ElectronInfo& Projection::operator[](unsigned int i)
{
    return config[i];
}

const ElectronInfo& Projection::operator[](unsigned int i) const
{
    return config[i];
}

Parity Projection::GetParity() const
{
    int sum = 0;
    for_each(config.begin(), config.end(), [&sum](const ElectronInfo& info)
    {   sum += info.L();
    } );

    if(sum%2 == 0)
        return Parity::even;
    else
        return Parity::odd;
}

int Projection::GetTwoM() const
{
    int sum = 0;
    for_each(config.begin(), config.end(), [&sum](const ElectronInfo& info)
    {   sum += (info.IsHole()? -info.TwoM(): info.TwoM());
    } );

    return sum;
}

std::string Projection::Name() const
{
    std::string name = "{";
    if(config.size())
        for(auto& electron: config)
            name.append(" " + electron.Name());
    else    // Vacuum
        name.append("0");
    name.append("}");

    return name;
}
}
