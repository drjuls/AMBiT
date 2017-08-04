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

unsigned int RelativisticConfigList::NumCSFsSmall() const
{
    unsigned int total = 0;
    for(auto it = begin(); it != small_end(); ++it)
    {
        total += it->NumCSFs();
    }

    return total;
}

RelativisticConfigList::iterator RelativisticConfigList::erase(iterator position)
{
    auto it = m_list.erase(position.base());
    if(it - m_list.begin() < Nsmall)
        Nsmall--;

    return iterator(it);
}

RelativisticConfigList::iterator RelativisticConfigList::erase(const_iterator first, const_iterator last)
{
    int start = first.base() - m_list.begin();
    int num_removed = last.base() - first.base();
    if(start < Nsmall)
        Nsmall -= mmin(Nsmall - start, num_removed);

    return iterator(m_list.erase(first.base(), last.base()));
}

RelativisticConfigList::iterator RelativisticConfigList::operator[](unsigned int i)
{
    return std::next(begin(), i);
}

RelativisticConfigList::const_iterator RelativisticConfigList::operator[](unsigned int i) const
{
    return std::next(begin(), i);
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

void RelativisticConfigList::unique()
{
    auto itsmall = std::next(m_list.begin(), Nsmall);
    itsmall = m_list.erase(std::unique(m_list.begin(), itsmall), itsmall);

    m_list.erase(std::unique(itsmall, m_list.end()), m_list.end());
    Nsmall = itsmall - m_list.begin();
}

void RelativisticConfigList::Read(FILE* fp)
{
    m_list.clear();
    unsigned int num_configs;
    fread(&num_configs, sizeof(unsigned int), 1, fp);
    m_list.reserve(num_configs);

    for(unsigned int i = 0; i < num_configs; i++)
    {
        RelativisticConfiguration config;
        config.Read(fp);
        m_list.push_back(config);
    }
}

void RelativisticConfigList::Write(FILE* fp) const
{
    unsigned int num_configs = m_list.size();
    fwrite(&num_configs, sizeof(unsigned int), 1, fp);

    for(const auto& relconfig: *this)
    {
        relconfig.Write(fp);
    }
}
