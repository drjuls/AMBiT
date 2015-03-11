#include "RelativisticConfigList.h"

unsigned int RelativisticConfigList::NumCSFs() const
{
    unsigned int total = 0;
    for(auto& rconfig: m_list)
    {
        total += rconfig.NumCSFs();
    }

    return total;
}

RelativisticConfigList::iterator RelativisticConfigList::operator[](unsigned int i)
{
    iterator it = begin();
    for(unsigned int j = 0; j < i; j++)
        it++;

    return it;
}

RelativisticConfigList::const_iterator RelativisticConfigList::operator[](unsigned int i) const
{
    const_iterator it = begin();
    for(unsigned int j = 0; j < i; j++)
        it++;

    return it;
}

RelativisticConfigList::const_projection_iterator RelativisticConfigList::projection_begin() const
{
    return const_projection_iterator(this);
}

RelativisticConfigList::const_projection_iterator RelativisticConfigList::projection_end() const
{
    const_iterator last = end();

    // Check size = 0 case
    if(begin() == last)
        return const_projection_iterator(this);

    last--;

    // NB: below we use const_projection_iterator(this, last->projection_end(), end(), 0) rather than
    //                  const_projection_iterator(this, last->projection_end(), end(), NumCSFs())
    //     even though the latter makes more sense. But the comparison is only made on
    //     const_projection_iterator.m_base, i.e. last->projection_end(), and so the CSF_index is meaningless
    return const_projection_iterator(this, last.projection_end(), end(), 0);
}

unsigned int RelativisticConfigList::projection_size() const
{
    unsigned int total = 0;
    for(auto& rconfig: m_list)
    {
        total += rconfig.projection_size();
    }

    return total;
}

void RelativisticConfigList::Read(FILE* fp)
{
    clear();
    unsigned int num_configs;
    fread(&num_configs, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < num_configs; i++)
    {
        RelativisticConfiguration config;
        config.Read(fp);
        add(config);
    }
}

void RelativisticConfigList::Write(FILE* fp) const
{
    unsigned int num_configs = size();
    fwrite(&num_configs, sizeof(unsigned int), 1, fp);

    for(const auto& relconfig: *this)
    {
        relconfig.Write(fp);
    }
}
