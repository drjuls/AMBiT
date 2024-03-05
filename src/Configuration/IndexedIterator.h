#ifndef INDEXED_ITERATOR_H
#define INDEXED_ITERATOR_H

#include <boost/iterator_adaptors.hpp>

namespace Ambit
{
/** An iterator class for when one needs an iterator over values or other objects and the position of those values in the sequence.
    This is typical for example when constructing a Hamiltonian matrix or when the position must be correlated on another array.
    The iterator is initialised with a starting index, is incremented or decremented with the iterator, and is accessible via the
    index() function.
 */
template <class Base, class CategoryOrTraversal = boost::use_default>
class indexed_iterator:
    public boost::iterator_adaptor<indexed_iterator<Base, CategoryOrTraversal>, Base, boost::use_default, CategoryOrTraversal>
{
private:
    struct enabler {};

public:
    indexed_iterator(): indexed_iterator::iterator_adaptor_(), m_index(0) {}

    explicit indexed_iterator(Base p, int start_index = 0):
        indexed_iterator::iterator_adaptor_(p), m_index(start_index) {}

    template <class OtherBase>
    indexed_iterator(indexed_iterator<OtherBase, CategoryOrTraversal> const& other,
                     typename boost::enable_if<boost::is_convertible<OtherBase,Base>, enabler>::type = enabler()):
        indexed_iterator::iterator_adaptor_(other.base()), m_index(other.m_index) {}

    int index() const { return m_index; }

private:
    // Needed by iterator_adaptor
    friend class boost::iterator_core_access;
    int m_index;

    void increment()
    {   this->base_reference()++;
        m_index++;
    }

    void decrement()
    {   this->base_reference()--;
        m_index--;
    }

    void advance(typename indexed_iterator::iterator_adaptor_::difference_type n)
    {   this->base_reference() += n;
        m_index += n;
    }
};

}
#endif
