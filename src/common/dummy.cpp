

#include <common/dummy.hpp>

int getSomeNumber(int x, int y)
{
    int z = x*y;
    // some random special condition (not covered by tests)
    if (z == -13)
    {
        z = 0;
    }
    return x*y;
}