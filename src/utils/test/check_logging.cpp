#include "utils/logging.h"

int foo(int x, int y)
{
  BL_INFO("function foo() was called");
  BL_DEBUG("x is " << x);
  BL_DEBUG("y is " << y);
  if (x == 0)
  {
    BL_ERROR("x is null");
  }

  int z = x * y;
  BL_DEBUG("z is " << z);

  BL_INFO("function foo() is done");
  return z;
}

int main(int argc, char *argv[])
{
  // init the logging
  LOG_INIT();

  BL_INFO("calling foo()");
  int blah = foo(12, 13);
  BL_DEBUG("blah is " << blah);
  USED_BY_LOGGER_ONLY(blah);

  BL_INFO("calling foo() with wrong arguments");
  int bar = foo(0, 0);
  BL_DEBUG("bar is " << bar);
  USED_BY_LOGGER_ONLY(bar);

  return 0;
}
