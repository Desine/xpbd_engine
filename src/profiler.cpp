#include "profiler.hpp"
#include "../tracy/public/TracyClient.cpp"




/*
compile on linux for x11
cd tracy/profiler
cmake -DLEGACY=ON -B profiler/build -S profiler
cmake --build profiler/build --config Release --parallel $(nproc)
*/