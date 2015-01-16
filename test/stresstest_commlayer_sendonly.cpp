
#include <mpi.h>
#include <omp.h>
#include <unistd.h> // for sleep!

#include <iostream>
#include <functional>

#include <io/CommunicationLayer.hpp>

//#define DEBUG(msg) std::cerr << msg << std::endl;

int my_rank;
std::atomic<int> msgs_received(0);
std::atomic<int> lookup_received(0);
std::atomic<int> answers_received(0);

template <bool ThreadLocal = false>
struct Tester
{
  const int ANSWER_TAG = 12;
  const int FIRST_TAG = 1;
  const int LOOKUP_TAG = 13;

  int generate_message(int srcRank, int dstRank)
  {
    return (srcRank + 1) * 100000 + (dstRank + 1);
  }

  void receivedCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    //DEBUG("Rank " << my_rank << " received " << count << " message from process: " << fromRank);

    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    // for all received requests, send the value from the lookup
    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != generate_message(fromRank, my_rank))
      {
        ERRORF("ERROR: message not as expected.  Expected: %d Actual: %d", generate_message(fromRank, my_rank), msgs[i]);
        ERRORF("\tmy rank: %d from rank %d message id = %d", my_rank, fromRank, msgs_received.load());

        exit(EXIT_FAILURE);
      }
      else
      {
        // DEBUG("SUCCESS: message received");
        msgs_received.fetch_add(1);
      }
    }
  }

  void lookup_callback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    // for all received requests, send the value from the lookup
    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != (fromRank+1)* (my_rank+1))
      {
        ERRORF("ERROR: LOOKUP message not as expected: %d expected %d", msgs[i], generate_message(fromRank, my_rank));
        exit(EXIT_FAILURE);
      } else {
        // DEBUG("SUCCESS: message received");
        lookup_received.fetch_add(1);

      }
      int msg = msgs[i] + 13;
      commLayer.sendMessage(&msg, sizeof(int), fromRank, ANSWER_TAG);
    }
  }

  void answer_callback(uint8_t* msg, std::size_t count, int fromRank)
  {
    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != (fromRank+1)* (my_rank+1) + 13)
      {
        ERRORF("ERROR: ANSWER message not as expected: %d expected %d", msgs[i], generate_message(my_rank, fromRank) + 1000);
        exit(EXIT_FAILURE);
      }
      answers_received.fetch_add(1);
    }
  }

  void test_comm_layer(int iters, int els)
  {
    DEBUG("Testing Comm Layer");
    DEBUG("Size: " << commLayer.getCommSize());
    DEBUG("Rank: " << commLayer.getCommRank());

    // set global rank
    my_rank = commLayer.getCommRank();
    using namespace std::placeholders;
    commLayer.addReceiveCallback(FIRST_TAG, std::bind(&Tester::receivedCallback, this, _1, _2, _3));
//    commLayer.addReceiveCallback(LOOKUP_TAG, std::bind(&Tester::lookup_callback, this, _1, _2, _3));
//    commLayer.addReceiveCallback(ANSWER_TAG, std::bind(&Tester::answer_callback, this, _1, _2, _3));

    commLayer.initCommunication();

    int nthreads = numThreads;
    int it = 0;
    int i = 0;
    for (; it < iters; ++it) {

      // R: src rank
      // T: thread id
      // I: iteration
      // D: dest rank
      // t: tag
      // i: message counter
      // M: message
      // L: recv message cont

      // start sending one message to each:
#pragma omp parallel for default(none) num_threads(nthreads) shared(els, my_rank, it, stdout)
      for (i = 0; i < els; ++i)
      {
        int msg;
        for (int j = 0; j < commSize; ++j)
        {
          msg = generate_message(my_rank, j);
          if (i == 0 || i == els - 1)
            DEBUGF("W R %d,\tT %d,\tI %d,\tD %d,\tt %d,\ti %d/%d,\tM %d", my_rank, omp_get_thread_num(), it, j, FIRST_TAG, i, els, msg);
          commLayer.sendMessage(&msg, sizeof(int), j, FIRST_TAG);
        }
      }


      DEBUGF("M R %d,\tT  ,\tI %d,\tD  ,\tt %d,\ti %d,\tM ,\tL%d PREFLUSH", my_rank, it, FIRST_TAG, els, msgs_received.load());
      commLayer.flush(FIRST_TAG);
      DEBUGF("M R %d,\tT  ,\tI %d,\tD  ,\tt %d,\ti %d,\tM ,\tL%d POSTFLUSH", my_rank, it, FIRST_TAG, els, msgs_received.load());


//      // check that all messages have been received
//      if (msgs_received.load() != els * commSize)
//      {
//  //      ERROR("ERROR: wrong amount of messages received in phase 1");
//  //      ERROR("received: " << msgs_received.load() << ", should: " << els * commLayer.getCommSize() * iters);
//        ERROR("M R " << my_rank << ",\tT " << " " << ",\tI " << it << ",\tD " << " " << ",\tt " << FIRST_TAG << ",\ti " << " " << ",\tM " << " ,\tL " << msgs_received.load() << "\tFAIL: expected " << els * commSize);
//        exit(EXIT_FAILURE);
//      }
//
//      msgs_received.store(0);
    }
    // call the finish function for this tag  //
    DEBUGF("M R %d,\tT  ,\tI %d,\tD  ,\tt %d,\ti %d,\tM ,\tL%d PREFINISH", my_rank, it, FIRST_TAG, els, msgs_received.load());
    commLayer.finish(FIRST_TAG);
    DEBUGF("M R %d,\tT  ,\tI %d,\tD  ,\tt %d,\ti %d,\tM ,\tL%d POSTFINISH", my_rank, it, FIRST_TAG, els, msgs_received.load());

    // check that all messages have been received
    if (msgs_received.load() != els * commSize * iters)
    {
      ERRORF("M R %d,\tT  ,\tI  ,\tD  ,\tt %d,\ti  ,\tM ,\tL%d, \tFAIL: expected %d", my_rank, FIRST_TAG, msgs_received.load(), els * commSize * iters);
//      exit(EXIT_FAILURE);
    }
    //std::cerr << "INDEX: " << msgs_received << std::endl;


    /* phase 2 communication */

//    // sending one message to each:
//#pragma omp parallel for default(none) num_threads(nthreads) shared(els, my_rank)
//    for (int i = 0; i < els; ++i)
//    {
//
//      for (int j = 0; j < commSize; ++j)
//      {
//        int msg = generate_message(my_rank, j);
//        if (i % 1000 == 0) DEBUG("Rank " << my_rank << "." << omp_get_thread_num() << " " << i << " of " << els << " Querying " << msg);
//        commLayer.sendMessage(&msg, sizeof(int), j, LOOKUP_TAG);
//      }
//    }
//
//    // flush both tags
//    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());
//    commLayer.flush(LOOKUP_TAG);
//    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());
//
//    // flush both tags
//    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());
//    commLayer.flush(LOOKUP_TAG);
//    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " query messages received = " << lookup_received.load());
//
//    //======== call flush twice helps to avoid error below.  HACK
//    // check that all messages have been received correctly
//    if (lookup_received != els * commSize)
//    {
//      std::cerr << "rank " << commRank << " FAIL: wrong amount of lookup messages received in phase 2" << std::endl;
//      std::cerr << "rank " << commRank << " received: " << lookup_received.load() << ", should: " << els * commSize << std::endl;
//      exit(EXIT_FAILURE);
//    }
//
//    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());
//    commLayer.flush(ANSWER_TAG);
//    DEBUG("rank " << commRank << "thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());
//
//    DEBUG("rank " << commRank << " thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());
//    commLayer.flush(ANSWER_TAG);
//    DEBUG("rank " << commRank << "thread " << omp_get_thread_num() << " answer messages received = " << answers_received.load());
//
//    //======== call flush twice helps to avoid error below.  HACK
//    // check that all messages have been received correctly
//    if (answers_received != els * commSize)
//    {
//      std::cerr << "rank " << commRank << " FAIL: wrong amount of answer messages received in phase 2" << std::endl;
//      std::cerr << "rank " << commRank << " received: " << answers_received.load() << ", should: " << els * commSize << std::endl;
//      exit(EXIT_FAILURE);
//    }

    INFOF("M R %d, SEND DONE. ", commRank);



//    commLayer.finish(FIRST_TAG);
//    commLayer.finishTag(LOOKUP_TAG);
//    commLayer.finishTag(ANSWER_TAG);
    commLayer.finishCommunication();

    //std::cerr << "LOOKUP: " << lookup_received << " ANSWERS: " << answers_received << std::endl;

    DEBUGF("This was a triumph.");
//    sleep(1);
    DEBUGF("I'm making a note here: HUGE SUCCESS.");
//    sleep(1);
    DEBUGF("It's hard to overstate my satisfaction.");

    INFOF("M R %d, ALL DONE. ", commRank);
  }

  Tester(MPI_Comm comm, int comm_size, int num_threads) :
    commLayer(comm, comm_size, num_threads), commSize(comm_size), numThreads(num_threads) {
    //commLayer.startThreads();
    MPI_Comm_rank(comm, &commRank);
  }

  bliss::io::CommunicationLayer<ThreadLocal> commLayer;

  int commSize;
  int commRank;
  int numThreads;
};

int main(int argc, char *argv[])
{
  int nthreads = 1;
  if (argc > 1) {
    nthreads = atoi(argv[1]);
  }


  int elems = 1536 * nthreads;
  if (argc > 2) {
    elems = atoi(argv[2]);
  }


  int iters = 10;
  if (argc > 3) {
    iters = atoi(argv[3]);
  }


  // set up MPI
  MPI_Init(&argc, &argv);

  // get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int p, rank;
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);

  /* code */
  {
    msgs_received.store(0);
    lookup_received.store(0);
    answers_received.store(0);

#if defined(THREADLOCAL)
    Tester<true> tester(comm, p, nthreads);
#else
    Tester<false> tester(comm, p, nthreads);
#endif
    tester.test_comm_layer(iters, elems);

    MPI_Barrier(comm);
  }

  // finalize MPI
  MPI_Finalize();
  return 0;
}

