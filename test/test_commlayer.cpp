
#include <mpi.h>
#include <omp.h>
#include <unistd.h> // for sleep!

#include <iostream>
#include <functional>

#include <io/CommunicationLayer.hpp>

//#define DEBUG(msg) std::cerr << msg << std::endl;

int my_rank;
std::atomic<int> msgs_received(0);
std::atomic<int> answers_received(0);
int iters;

struct Tester
{
  const int ANSWER_TAG = 12;
  const int FIRST_TAG = 1;
  const int LOOKUP_TAG = 13;

  int generate_message(int srcRank, int dstRank)
  {
    return (srcRank + 1) * iters * 10 + (dstRank + 1);
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
        ERROR("ERROR: message not as expected.  Expected: " << generate_message(fromRank, my_rank) << " Actual: "<< msgs[i] << "");
        ERROR("my rank: " << my_rank << " from rank " << fromRank);
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
        ERROR("ERROR: LOOKUP message not as expected: " << msgs[i]);
        exit(EXIT_FAILURE);
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
        ERROR("ERROR: ANSWER message not as expected: " << msgs[i]);
        exit(EXIT_FAILURE);
      }
      answers_received.fetch_add(1);
    }
  }

  void test_comm_layer(int repeat_sends=10000, int nthreads = 1)
  {
    DEBUG("Testing Comm Layer");
    DEBUG("Size: " << commLayer.getCommSize());
    DEBUG("Rank: " << commLayer.getCommRank());

    // set global rank
    my_rank = commLayer.getCommRank();
    using namespace std::placeholders;
    commLayer.addReceiveCallback(FIRST_TAG, std::bind(&Tester::receivedCallback, this, _1, _2, _3));
    commLayer.addReceiveCallback(LOOKUP_TAG, std::bind(&Tester::lookup_callback, this, _1, _2, _3));
    commLayer.addReceiveCallback(ANSWER_TAG, std::bind(&Tester::answer_callback, this, _1, _2, _3));

    commLayer.initCommunication();

    int iters = 2;

    for (int it = 0; it < iters; ++it) {

      // start sending one message to each:
#pragma omp parallel for default(none) num_threads(nthreads) shared(repeat_sends, my_rank)
      for (int l = 0; l < repeat_sends; ++l)
      {

        for (int i = 0; i < commLayer.getCommSize(); ++i)
        {
          int msg = generate_message(my_rank, i);
//          if (l % 10000 == 0) DEBUG("Thread " << tid << " Sending " << msg << " to " << i);
          commLayer.sendMessage(&msg, sizeof(int), i, FIRST_TAG);
        }
      }

      if (commLayer.getCommRank() == 0) {
        sleep(1);
      }

      DEBUG("thread " << omp_get_thread_num() << " messages received = " << msgs_received.load());

      // call the flush function for this tag
      commLayer.flush(FIRST_TAG);

      //=== debug messages show that there are no messages waiting.  so where are the missing messages?
      DEBUG("thread " << omp_get_thread_num() << " flushed. messages received = " << msgs_received.load());


    }


    // check that all messages have been received
    if (msgs_received.load() != repeat_sends * commLayer.getCommSize() * iters)
    {
      ERROR("ERROR: wrong amount of messages received in phase 1");
      ERROR("received: " << msgs_received.load() << ", should: " << repeat_sends * commLayer.getCommSize() * iters);
      exit(EXIT_FAILURE);
    }


    /* phase 2 communication */
    msgs_received.store(0);


    // sending one message to each:
#pragma omp parallel for default(none) num_threads(nthreads) shared(repeat_sends, my_rank)
    for (int l = 0; l < repeat_sends; ++l)
    {

      for (int i = 0; i < commLayer.getCommSize(); ++i)
      {
        int msg = (my_rank+1)*(i+1);
//        if (l % 10000 == 0) DEBUG("Thread " << tid << " Querying " << msg << " from " << i);
        commLayer.sendMessage(&msg, sizeof(int), i, LOOKUP_TAG);
      }
    }

    // flush both tags
    DEBUG("thread " << omp_get_thread_num() << " query messages received = " << msgs_received.load());
    commLayer.flush(LOOKUP_TAG);
    DEBUG("thread " << omp_get_thread_num() << " query messages received = " << msgs_received.load());

    DEBUG("thread " << omp_get_thread_num() << " answer messages received = " << msgs_received.load());
    commLayer.flush(ANSWER_TAG);
    DEBUG("thread " << omp_get_thread_num() << " answer messages received = " << msgs_received.load());

    // check that all messages have been received correctly
    if (answers_received != repeat_sends * commLayer.getCommSize())
    {
      ERROR("ERROR: wrong amount of messages received in phase 2");
      ERROR("received: " << answers_received.load() << ", should: " << repeat_sends * commLayer.getCommSize());
      exit(EXIT_FAILURE);
    }

//    commLayer.finish(FIRST_TAG);
//    commLayer.finishTag(LOOKUP_TAG);
//    commLayer.finishTag(ANSWER_TAG);
    commLayer.finishCommunication();

    DEBUG("SENT: " << msgs_received << " ANSWERS: " << answers_received);

    DEBUG("This was a triumph.");
    sleep(1);
    DEBUG("I'm making a note here: HUGE SUCCESS.");
    sleep(1);
    DEBUG("It's hard to overstate my satisfaction.");

    fflush(stdout);
  }

  Tester(MPI_Comm comm, int comm_size) : commLayer(comm, comm_size) {
    //commLayer.startThreads();
  }

  bliss::io::CommunicationLayer commLayer;
};

int main(int argc, char *argv[])
{
  int nthreads = 1;
  if (argc > 1) {
    nthreads = atoi(argv[1]);
  }

  iters = 10 * nthreads;
  if (argc > 2) {
    iters = atoi(argv[2]);
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
  Tester tester(comm, p);
  tester.test_comm_layer(iters, nthreads);

  MPI_Barrier(comm);
  }
  // finalize MPI
  MPI_Finalize();
  return 0;
}

