#ifndef SORTED_LIST_H
#define SORTED_LIST_H

#include <list>

namespace Ambit
{
template <class Value>
class DefaultComparator
{
public:
    bool operator()(const Value& first, const Value& second) { return (first < second); }
};

/** Template for an extensible, sorted forward_list with a simplified STL interface.
    Compare::operator() is a binary predicate that defines a strict weak ordering between two Value objects.
 */
template <class Value, class Compare =  DefaultComparator<Value> >
class SortedList
{
protected:
    typedef SortedList<Value, Compare> BaseSortedList;
    typedef Compare BaseComparator;

public:
    SortedList() {}
    SortedList(const SortedList<Value, Compare>& other): m_list(other.m_list) {}
    SortedList(SortedList<Value, Compare>&& other): m_list(other.m_list) {}
    SortedList(const Value& val) { m_list.push_back(val); }
    SortedList(Value&& val) { m_list.push_back(val); }
    virtual ~SortedList() {}

    typedef typename std::list<Value>::iterator iterator;
    typedef typename std::list<Value>::const_iterator const_iterator;

    void add(const Value& val) { merge(BaseSortedList(val)); }
    void add(Value&& val)  { merge(BaseSortedList(val)); }

    iterator begin() { return m_list.begin(); }
    const_iterator begin() const { return m_list.begin(); }

    void clear() { m_list.clear(); }

    bool empty() const { return m_list.empty(); }

    iterator end() { return m_list.end(); }
    const_iterator end() const { return m_list.end(); }

    iterator erase(iterator position) { return m_list.erase(position); }
    iterator erase(iterator first, iterator last) { return m_list.erase(first, last); }

    /** PRE: size() > 0 */
    Value& front() { return m_list.front(); }
    const Value& front() const { return m_list.front(); }

    void merge(BaseSortedList& other) { m_list.merge(other.m_list, Compare()); }
    void merge(BaseSortedList&& other) { m_list.merge(other.m_list, Compare()); }

    BaseSortedList& operator=(const BaseSortedList& other) { m_list = other.m_list; return *this; }
    BaseSortedList& operator=(BaseSortedList&& other) { m_list = other.m_list; return *this; }

    unsigned int size() const { return m_list.size(); }

    void unique() { return m_list.unique(); }

protected:
    std::list<Value> m_list;
};


}
#endif
