#include <mpi.h>

#include <iostream>

#include <unistd.h> // for sleep!


#include <io/CommunicationLayer.hpp>
#include <index/distributed_map.hpp>

//#define DEBUG(msg) std::cerr << msg << std::endl;

int glCommSize;

void receiveAnswer(std::pair<int, bliss::index::count_t>& answer)
{
  if (answer.second != 1000u * glCommSize * answer.first)
  {
    std::cerr << "ERROR: distributed count is wrong: received=" << answer.second << ", expected=" << (glCommSize* answer.first*1000) << std::endl;
  }
  else
  {
    std::cerr << "SUCCESS!" << std::endl;
  }

  std::cerr << std::flush;
}


void test_map(MPI_Comm& comm, int repeat=1000)
{
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);
  glCommSize = p;
  bliss::index::distributed_counting_map<int, bliss::io::CommunicationLayer> counting_map(comm, p);
  counting_map.setLookupAnswerCallback(std::function<void(std::pair<int, bliss::index::count_t>&)>(&receiveAnswer));

  sleep(1);

  for (int i = 0; i < p; ++i)
  {
    for (int j = 0; j < i*repeat; ++j)
    {
      counting_map.insert(i);
    }
  }

  counting_map.flush();

  for (int i = 0; i < p; ++i)
  {
    counting_map.asyncLookup(i);
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
  test_map(comm);

  MPI_Barrier(comm);

  // finalize MPI
  MPI_Finalize();
  return 0;
}

