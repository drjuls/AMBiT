#include "AngularData.h"
#include "Include.h"
#include "RelativisticConfiguration.h"
#include <numeric>

AngularData::AngularData(const RelativisticConfiguration& config, int two_m): CSFs(nullptr)
{
    GenerateProjections(config, two_m);
}

AngularData::AngularData(const RelativisticConfiguration& config, int two_m, int two_j): CSFs(nullptr)
{
    GenerateProjections(config, two_m);
    GenerateCSFs(two_j);
}

AngularData::~AngularData()
{
    if(CSFs)
        delete[] CSFs;
}

bool AngularData::ProjectionCompare(const std::vector<int>& first, const std::vector<int>& second)
{
    auto it_first = first.begin();
    auto it_second = second.begin();
    while(it_first != first.end() && it_second != second.end())
    {
        if(*it_first > *it_second)
            return false;
        else if(*it_first < *it_second)
            return true;

        it_first++;
        it_second++;
    }

    if(first.size() >= second.size())
        return false;
    else
        return true;
}

bool AngularData::GenerateProjections(const RelativisticConfiguration& config, int two_m)
{
    // Start with maximum M.
    // At each step populate all possible projections with M-1.
    // Stop when desired M is reached.

    // Create set of boxes, each with a list of projections to be filled.
    int projection_size = config.ExcitationNumber();
    std::vector<bool> next_in_same_orbital(projection_size);
    std::vector<int> min_projection(projection_size);

    int count = 0;
    for(auto config_it = config.begin(); config_it != config.end(); config_it++)
    {
        for(int i = 0; i < abs(config_it->second)-1; i++)
        {   min_projection[count] = config_it->first.Kappa();
            next_in_same_orbital[count] = true;
            count++;
        }

        min_projection[count] = -config_it->first.TwoJ();
        next_in_same_orbital[count] = false;
        count++;
    }

    // Get maximal projection
    std::vector<int> proj(projection_size);
    proj[0] = -min_projection[0];
    for(int i = 1; i < projection_size; i++)
    {
        if(next_in_same_orbital[i-1])
            proj[i] = proj[i-1] - 2;
        else
            proj[i] = -min_projection[i];
    }
    proj.shrink_to_fit();

    int current_two_m = 0;
    std::accumulate(proj.begin(), proj.end(), current_two_m);

    std::vector< std::list< std::vector<int> > > boxes(projection_size);
    if(current_two_m >= two_m)
        boxes[0].push_back(proj);

    // Iterative method to find all projections with correct m.
    while(current_two_m > two_m)
    {
        std::vector< std::list< std::vector<int> > > next_boxes(projection_size);

        // For each projection in boxes[i], apply M_j-- for i <= j and
        // store the result in next_boxes[j].
        for(int i = 0; i < projection_size; i++)
        {
            for(std::vector<int>& existing_proj: boxes[i])
            {
                for(int j = i; j < projection_size; j++)
                {
                    std::vector<int> new_proj(existing_proj);
                    if(!next_in_same_orbital[j] ||
                       new_proj[j] - new_proj[j+1] > 2)
                    {
                        new_proj[j] -= 2;
                        if(new_proj[j] >= min_projection[j])
                            next_boxes[j].push_back(new_proj);
                    }
                }
            }
        }

        boxes.swap(next_boxes);
        current_two_m -= 2;
    }

    // Sort and merge projections lists
    projections.clear();
    for(auto& list: boxes)
    {   list.sort(ProjectionCompare);
        projections.merge(list);
    }

    projections.unique();

    return (projections.size() != 0);
}

bool AngularData::GenerateCSFs(int two_j)
{
    return true;
}

AngularDataLibrary::AngularDataLibrary(int particle_number, int two_m, int two_j)
{
    filename = itoa(particle_number) + "." + itoa(two_m) + "." + itoa(two_j) + ".angular";
}

pAngularData AngularDataLibrary::operator[](const RelativisticConfiguration& config)
{
    KeyType mykey(GenerateKey(config));
    return library[mykey];
}

AngularDataLibrary::KeyType AngularDataLibrary::GenerateKey(const RelativisticConfiguration& config) const
{
    KeyType key;
    for(auto& config_it: config)
    {   key.push_back(config_it.first.Kappa());
        key.push_back(config_it.second);
    }
    return key;
}
