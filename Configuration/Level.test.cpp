#include "Level.h"
#include "gtest/gtest.h"

TEST(LevelTester, LevelIteratorTest)
{
    pLevelMap levels = pLevelMap(new LevelMap());
    Symmetry sym(0, Parity::even);
    int i;
    for(i = 0; i < 3; i++)
        (*levels)[LevelID(sym, i)] = nullptr;
    sym = Symmetry(2, Parity::even);
    for(i = 0; i < 2; i++)
        (*levels)[LevelID(sym, i)] = nullptr;
    sym = Symmetry(2, Parity::odd);
    for(i = 0; i < 4; i++)
        (*levels)[LevelID(sym, i)] = nullptr;

    sym = Symmetry(2, Parity::even);
    int count = 0;
    LevelMap::symmetry_iterator it = levels->begin(sym);
    while(it != levels->end(sym))
    {
        EXPECT_TRUE(it->first.GetSymmetry() == sym);
        count++;
        it++;
    }
    EXPECT_EQ(2, count);

    count = 0;
    sym = Symmetry(2, Parity::odd);
    it = levels->begin(sym);
    while(it.base() != it.end())
    {
        EXPECT_TRUE(it->first.GetSymmetry() == sym);
        count++;
        it++;
    }
    EXPECT_EQ(4, count);

    EXPECT_EQ(9, levels->size());
}
