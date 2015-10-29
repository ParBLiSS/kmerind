
#include <mpi.h>
//#include <unistd.h> // for sleep!, for gethosthanme

#include <functional>

#include "config/relacy_config.hpp"
#include "utils/logging.h"
#include "io/communication_layer.hpp"





// template parameter '2' is number of threads
template<typename T, int nThreads, int nelems>
struct commlayerSendOnlyTest : rl::test_suite<commlayerSendOnlyTest<T, nThreads, nelems>, nThreads>
{
    bliss::io::CommunicationLayer<true>* commLayer;
    MPI_Comm comm;
    int p, rank;
    char hostname[256];

    VAR_T(int) msgs_recved;
    VAR_T(int) msgs_error;
    const int FIRST_TAG = 1;

    std::vector<T> msgs;


    T generate_message(const int srcRank, const int dstRank)
    {
      return (srcRank + 1) * 100000 + (dstRank + 1);
    }

    void receivedCallback(uint8_t* msg, const std::size_t count, const int fromRank)
    {
      //DEBUG("Rank " << rank << " received " << count << " message from process: " << fromRank);

      // first: deserialize
      T* msgs = reinterpret_cast<T*>(msg);
      size_t msg_count = count / sizeof(T);

      DEBUGF("R %d <- %d, Callback with %lu elements", rank, fromRank, msg_count);
      T expected = generate_message(fromRank, rank);
      // for all received requests, send the value from the lookup
      for (size_t i = 0; i < msg_count; ++i)
      {
        // check that the message is as expected
        if(msgs[i] != expected)
        {
          //ERRORF("ERROR: message not as expected.  Expected: %d Actual: %d, my rank: %d from rank %d message id = %d, count %lu / %lu",
          //    expected, msgs[i], rank, fromRank, msgs_received.load(), i, msg_count);
          ++VAR(msgs_error);
        }
        else
        {
          // DEBUG("SUCCESS: message received");
          ++VAR(msgs_recved);
        }
      }
    }



    // executed in single thread before main thread function
    void before()
    {
        // get communicator size and my rank
      comm = MPI_COMM_WORLD;
      MPI_Comm_size(comm, &p);
      MPI_Comm_rank(comm, &rank);

//      // get host hame and print out
//      memset(hostname, 0, 256);
//      gethostname(hostname, 256);
//      INFOF("Rank %d hostname [%s]\n", rank, hostname);

      VAR(msgs_recved) = 0;
      VAR(msgs_error) = 0;
      // initialize variables
      for (int i = 0; i < nThreads; ++i) {

      }

      msgs.resize(p * nelems, 0);

      // set up the commlayer,
      commLayer = new bliss::io::CommunicationLayer<true>(comm, p, nThreads);

      using namespace std::placeholders;
      commLayer->addReceiveCallback(FIRST_TAG,
                                   std::bind(&commlayerSendOnlyTest<T, nThreads, nelems>::receivedCallback,
                                             this, _1, _2, _3));

      commLayer->initCommunication();

      MPI_Barrier(comm);
    }

    // main thread function
    void thread(unsigned thread_index)
    {
      size_t idx;
      for (int i = thread_index; i < nelems; i+= nThreads) {
        for (int j = 0; j < p; ++j)
        {
          idx = i * p + j;
          msgs[idx] = generate_message(rank, j);
          commLayer->sendMessage(&(msgs[idx]), 1, j, FIRST_TAG);
//          if (i == 0 || i == els - 1)
//            DEBUGF("W R %d,\tT %d,\tI %d,\tD %d,\tt %d,\ti %d/%d,\tM %d", rank, omp_get_thread_num(), it, j, FIRST_TAG, i, els, msgs[j]);

//          if ((msgs[idx] / 100000 != rank + 1) || (msgs[idx] % 1000 != j + 1)) {
//            ERRORF("ERROR: DEBUG: build not correct: %d -> %d u= %d", rank, j, msgs[idx]);
//          }

        }
      }

    }

    // executed in single thread after main thread function
    void after()
    {

      commLayer->finish(FIRST_TAG);

      assert(VAR(msgs_recved) == nelems * p);

      // R: src rank
      // T: thread id
      // I: iteration
      // D: dest rank
      // t: tag
      // i: message counter
      // M: message
      // L: recv message cont
      //DEBUGF("M R %d,\tT  ,\tD  ,\tt %d,\ti %d,\tM ,\tL%d ", rank, FIRST_TAG, nelems, VAR(msgs_recved));

      commLayer->finishCommunication();

      assert(VAR(msgs_recved) == nelems * p);


      MPI_Barrier(comm);

      delete commLayer;
    }

    // executed in single thread after every 'visible' action in main threads
    // disallowed to modify any state
    void invariant()
    {
    }
};


int main(int argc, char *argv[])
{

  int iterations = 1000;
  if (argc > 1)
  {
    iterations = atoi(argv[1]);
  }

  rl::test_params parm;
  parm.search_type = rl::sched_random;
  parm.iteration_count = iterations;
  parm.execution_depth_limit = 4000;




  // set up MPI
  MPI_Init(&argc, &argv);


    rl::simulate<commlayerSendOnlyTest<int,  1,  1 * 1536> >(parm);
    rl::simulate<commlayerSendOnlyTest<int,  2,  2 * 1536> >(parm);
    rl::simulate<commlayerSendOnlyTest<int,  3,  3 * 1536> >(parm);
    rl::simulate<commlayerSendOnlyTest<int,  4,  4 * 1536> >(parm);
    rl::simulate<commlayerSendOnlyTest<int,  8,  8 * 1536> >(parm);
    rl::simulate<commlayerSendOnlyTest<int, 16, 16 * 1536> >(parm);


  // finalize MPI
  MPI_Finalize();
  return 0;
}

