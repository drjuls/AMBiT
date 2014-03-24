#include "Include.h"
#include "Projection.h"
#include "RelativisticConfiguration.h"
#include "HartreeFock/NonRelInfo.h"

Projection::Projection(const RelativisticConfiguration& relconfig, std::vector<int>& TwoMs)
{
    config.reserve(relconfig.ExcitationNumber());

    int i = 0;
    for(auto& relconfig_it: relconfig)
    {
        const OrbitalInfo& orbital = relconfig_it.first;
        config[i] = ElectronInfo(orbital.PQN(), orbital.Kappa(), TwoMs[i], (relconfig_it.second < 0));
    }

    config.shrink_to_fit();
}

Projection::Projection(const Projection& other): config(other.config)
{
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
    for_each(config.begin(), config.end(), [&](const ElectronInfo& info)
    {   sum += info.L();
    } );

    if(sum%2 == 0)
        return even;
    else
        return odd;
}

int Projection::GetTwoM() const
{
    int sum = 0;
    for_each(config.begin(), config.end(), [&](const ElectronInfo& info)
    {   sum += info.TwoM();
    } );

    return sum;
}

//bool Projection::operator<(const Projection& other) const
//{
//    const std::vector<ElectronInfo>& first_config(Config);
//    const std::vector<ElectronInfo>& second_config(other.Config);
//
//    //Sort(first_config);
//    //Sort(second_config);
//
//    std::vector<ElectronInfo>::const_iterator first = first_config.begin();
//    std::vector<ElectronInfo>::const_iterator second = second_config.begin();
//
//    while((first != first_config.end()) && (second != second_config.end()))
//    {
//        // Order by electron info
//        if(*first < *second)
//            return true;
//        else if(*second < *first)
//            return false;
//        
//        first++;
//        second++;
//    }
//
//    if((first == first_config.end()) && (second != second_config.end()))
//        return true;
//    else return false;
//}
//
//bool Projection::operator==(const Projection& other) const
//{
//    const std::vector<ElectronInfo>& first_config(Config);
//    const std::vector<ElectronInfo>& second_config(other.Config);
//
//    std::vector<ElectronInfo>::const_iterator first = first_config.begin();
//    std::vector<ElectronInfo>::const_iterator second = second_config.begin();
//
//    while(first != first_config.end())
//    {
//        if((second == second_config.end()) || (*first != *second))
//            return false;
//        
//        first++;
//        second++;
//    }
//
//    if(second != second_config.end())
//        return false;
//    else 
//        return true;
//}

std::string Projection::Name() const
{
    std::string name;
    for(auto electron : config)
        name.append(" " + electron.Name());

    return name;
}

int Projection::GetProjectionDifferences(const Projection& p1, const Projection& p2, unsigned int* diff)
{
    unsigned int i=0, j=0;
    unsigned int diff1[3], diff2[3];
    unsigned int diff1_size=0, diff2_size=0;

    while((i < p1.size()) && (j < p2.size()) && (diff1_size <= 2) && (diff2_size <= 2))
    {
        if(p1[i] == p2[j])
        {   i++;
            j++;
        }
        else if(p1[i] < p2[j])
        {   diff1[diff1_size++] = i;
            i++;
        }
        else
        {   diff2[diff2_size++] = j;
            j++;
        }
    }
    while((i < p1.size()) && (diff1_size <= 2))
    {   diff1[diff1_size++] = i;
        i++;
    }
    while((j < p2.size()) && (diff2_size <= 2))
    {   diff2[diff2_size++] = j;
        j++;
    }

    if((diff1_size > 2) || (diff2_size > 2))
        return 3;

    if(diff1_size < diff2_size)
        diff1[diff1_size++] = p1.size();
    else if(diff2_size < diff1_size)
        diff2[diff2_size++] = p2.size();

    if(diff1_size != diff2_size)
        return 3;
    else if(diff1_size == 0)
        return 0;

    int sort = abs(int(diff1[0]) - int(diff2[0]));
    
    //Copy first difference
    diff[0] = diff1[0];
    diff[1] = diff2[0];

    if(diff1_size == 2)
    {   sort += abs(int(diff1[1]) - int(diff2[1]));

        // Copy second difference
        diff[2] = diff1[1];
        diff[3] = diff2[1];
    }

    if(sort%2 == 0)
        return int(diff1_size);
    else
        return -int(diff1_size);
}

int Projection::GetProjectionDifferences3(const Projection& p1, const Projection& p2, unsigned int* diff)
{
    unsigned int i=0, j=0;
    unsigned int diff1[4], diff2[4];
    unsigned int diff1_size=0, diff2_size=0;

    while((i < p1.size()) && (j < p2.size()) && (diff1_size <= 3) && (diff2_size <= 3))
    {
        if(p1[i] == p2[j])
        {   i++;
            j++;
        }
        else if(p1[i] < p2[j])
        {   diff1[diff1_size++] = i;
            i++;
        }
        else
        {   diff2[diff2_size++] = j;
            j++;
        }
    }
    while((i < p1.size()) && (diff1_size <= 3))
    {   diff1[diff1_size++] = i;
        i++;
    }
    while((j < p2.size()) && (diff2_size <= 3))
    {   diff2[diff2_size++] = j;
        j++;
    }
    if((diff1_size != diff2_size) || (diff1_size > 3))
        return 4;
    else if(diff1_size == 0)
        return 0;

    int sort = abs(int(diff1[0]) - int(diff2[0]));
    
    //Copy first difference
    diff[0] = diff1[0];
    diff[1] = diff2[0];

    if(diff1_size >= 2)
    {   sort += abs(int(diff1[1]) - int(diff2[1]));

        // Copy second difference
        diff[2] = diff1[1];
        diff[3] = diff2[1];
    }

    if(diff1_size == 3)
    {   sort += abs(int(diff1[2]) - int(diff2[2]));

        // Copy third difference
        diff[4] = diff1[2];
        diff[5] = diff2[2];
    }

    if(sort%2 == 0)
        return int(diff1_size);
    else
        return -int(diff1_size);
}
