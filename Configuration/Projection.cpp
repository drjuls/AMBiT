#include "Include.h"
#include "Projection.h"
#include "HartreeFock/NonRelInfo.h"

void Projection::Add(const ElectronInfo& info)
{
    Config.push_back(info);
}

ElectronInfo& Projection::operator[](unsigned int i)
{
    return Config[i];
}

const ElectronInfo& Projection::operator[](unsigned int i) const
{
    return Config[i];
}

bool Projection::Sort()
{
    return Sort(Config);
}

bool Projection::Sort(std::vector<ElectronInfo>& config) const
{
    unsigned int sort = 0;
    std::vector<ElectronInfo>::iterator m_it = config.begin();
    while(m_it != config.end())
    {
        // Swap sort
        std::vector<ElectronInfo>::iterator r_it = m_it;
        r_it++;
        while(r_it != config.end())
        {
            if(*r_it < *m_it)
            {
                ElectronInfo temp = *r_it;
                *r_it = *m_it;
                *m_it = temp;

                sort++;
            }
            r_it++;
        }
        m_it++;
    }

    if(sort%2)
        return true;
    else
        return false;
}

Parity Projection::GetParity() const
{
    std::vector<ElectronInfo>::const_iterator m_it = Config.begin();
    unsigned int sum = 0;
    while(m_it != Config.end())
    {
        sum += m_it->L();
        m_it++;
    }

    if(sum%2 == 0)
        return even;
    else
        return odd;
}

int Projection::GetTwoM() const
{
    std::vector<ElectronInfo>::const_iterator m_it = Config.begin();
    int sum = 0;
    while(m_it != Config.end())
    {
        sum += m_it->TwoM();
        m_it++;
    }

    return sum;
}

bool Projection::operator<(const Projection& other) const
{
    const std::vector<ElectronInfo>& first_config(Config);
    const std::vector<ElectronInfo>& second_config(other.Config);

    //Sort(first_config);
    //Sort(second_config);

    std::vector<ElectronInfo>::const_iterator first = first_config.begin();
    std::vector<ElectronInfo>::const_iterator second = second_config.begin();

    while((first != first_config.end()) && (second != second_config.end()))
    {
        // Order by electron info
        if(*first < *second)
            return true;
        else if(*second < *first)
            return false;
        
        first++;
        second++;
    }

    if((first == first_config.end()) && (second != second_config.end()))
        return true;
    else return false;
}

bool Projection::operator==(const Projection& other) const
{
    const std::vector<ElectronInfo>& first_config(Config);
    const std::vector<ElectronInfo>& second_config(other.Config);

    std::vector<ElectronInfo>::const_iterator first = first_config.begin();
    std::vector<ElectronInfo>::const_iterator second = second_config.begin();

    while(first != first_config.end())
    {
        if((second == second_config.end()) || (*first != *second))
            return false;
        
        first++;
        second++;
    }

    if(second != second_config.end())
        return false;
    else 
        return true;
}

std::string Projection::Name() const
{
    std::vector<ElectronInfo>::const_iterator m_it = Config.begin();
    std::string name;
    while(m_it != Config.end())
    {
        name.append(" " + m_it->Name());
        m_it++;
    }
    return name;
}

int Projection::GetProjectionDifferences(const Projection& p1, const Projection& p2, unsigned int* diff)
{
    unsigned int i=0, j=0;
    unsigned int diff1[3], diff2[3];
    unsigned int diff1_size=0, diff2_size=0;

    while((i < p1.Size()) && (j < p2.Size()) && (diff1_size <= 2) && (diff2_size <= 2))
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
    while((i < p1.Size()) && (diff1_size <= 2))
    {   diff1[diff1_size++] = i;
        i++;
    }
    while((j < p2.Size()) && (diff2_size <= 2))
    {   diff2[diff2_size++] = j;
        j++;
    }
    if((diff1_size != diff2_size) || (diff1_size > 2))
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

    while((i < p1.Size()) && (j < p2.Size()) && (diff1_size <= 3) && (diff2_size <= 3))
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
    while((i < p1.Size()) && (diff1_size <= 3))
    {   diff1[diff1_size++] = i;
        i++;
    }
    while((j < p2.Size()) && (diff2_size <= 3))
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

Configuration Projection::GetNonRelConfiguration() const
{
    Configuration ret;

    for(unsigned int i = 0; i < Config.size(); i++)
    {
        ret.AddSingleParticle(NonRelInfo(Config[i].PQN(), Config[i].L()));
    }

    return ret;
}
