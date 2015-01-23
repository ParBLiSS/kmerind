#include <mpi.h>

#include <iostream>

#include <unistd.h> // for sleep!


#include <io/CommunicationLayer.hpp>
#include <index/distributed_map.hpp>

//#define DEBUG(msg) std::cerr << msg << std::endl;

int glCommSize;
int repeats;


void receiveAnswer(std::pair<int, bliss::index::count_t>& answer)
{
  int key = answer.first;
  int count = answer.second;
  if (count != repeats * (key+1) * glCommSize)
  {
    std::cerr << "ERROR: distributed count is wrong: received=" << count << ", expected=" << ((key+1) * repeats * glCommSize) << std::endl;
  }
  else
  {
    std::cerr << "SUCCESS!" << std::endl;
  }

  std::cerr << std::flush;
}

template<bool ThreadLocal = true>
void test_map(MPI_Comm& comm, int nthreads)
{
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);
  glCommSize = p;

  //printf("INIT COUNTING MAP\n");
  bliss::index::distributed_counting_map<int, bliss::io::CommunicationLayer<ThreadLocal> > counting_map(comm, p, nthreads);
  //printf("REGISTER COUNTING MAP CALLBACK\n");
  counting_map.setLookupAnswerCallback(std::function<void(std::pair<int, bliss::index::count_t>&)>(&receiveAnswer));

  counting_map.init();

  //sleep(1);

  for (int i = 0; i < p; ++i)
  {
    for (int j = 0; j < (i+1)*repeats; ++j)
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
  int nthreads = 1;
  if (argc > 1) {
    nthreads = atoi(argv[1]);
  }

  // set up MPI
  MPI_Init(&argc, &argv);

  repeats = 1000;
  if (argc > 1)
	  repeats = atoi(argv[1]);

  // get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);

  /* code */
//  {
//  test_map<false>(comm, nthreads);
//
//  MPI_Barrier(comm);
//  }

  {
  test_map<true>(comm, nthreads);

  MPI_Barrier(comm);
  }


  // finalize MPI
  MPI_Finalize();
  return 0;
}

