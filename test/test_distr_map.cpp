#include <mpi.h>

#include <iostream>

#include <io/CommunicationLayer.hpp>
#include <index/distributed_map.hpp>

#define DEBUG(msg) std::cerr << msg << std::endl;

void test_map(MPI_Comm& comm, int repeat=1000)
{
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);
  CommunicationLayer commLayer(comm, p);
  distributed_counting_map<int, CommunicationLayer> counting_map(commLayer, comm);

  for (int i = 0; i < p; ++i)
  {
    for (int j = 0; j < i*repeat; ++j)
    {
      counting_map.remoteInsert(i);
    }
  }

  counting_map.flush();
}

int main(int argc, char *argv[])
{
  // set up MPI
  MPI_Init(&argc, &argv);

  // get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);

  /* code */
  // TODO

  MPI_Barrier(comm);

  // finalize MPI
  MPI_Finalize();
  return 0;
}

