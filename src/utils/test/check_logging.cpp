/*
 * Copyright 2015 Georgia Institute of Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
