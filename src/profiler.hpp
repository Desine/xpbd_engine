#pragma once

#define TRACY_ENABLE

#include "Tracy.hpp"

/*
compile on linux for x11
cd tracy/profiler
cmake -DLEGACY=ON -B profiler/build -S profiler
cmake --build profiler/build --config Release --parallel $(nproc)
*/