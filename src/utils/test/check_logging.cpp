#include "utils/logging.h"

int foo(int x, int y)
{
  INFO("function foo() was called");
  DEBUG("x is " << x);
  DEBUG("y is " << y);
  if (x == 0)
  {
    ERROR("x is null");
  }

  int z = x * y;
  DEBUG("z is " << z);

  INFO("function foo() is done");
  return z;
}

int main(int argc, char *argv[])
{
  // init the logging
  LOG_INIT();

  INFO("calling foo()");
  int blah = foo(12, 13);
  DEBUG("blah is " << blah);
  USED_BY_LOGGER_ONLY(blah);

  INFO("calling foo() with wrong arguments");
  int bar = foo(0, 0);
  DEBUG("bar is " << bar);
  USED_BY_LOGGER_ONLY(bar);

  return 0;
}
