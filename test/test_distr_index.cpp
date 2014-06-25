
#include <mpi.h>
#include <unistd.h> // for sleep!

#include <iostream>
#include <functional>

#include <io/CommunicationLayer.hpp>

#define DEBUG(msg) std::cerr << msg << std::endl;

int my_rank;
volatile int msgs_received = 0;
volatile int answers_received = 0;

struct Tester
{
  const int ANSWER_TAG = 12;
  const int FIRST_TAG = 1;
  const int LOOKUP_TAG = 13;

  int generate_message(int srcRank, int dstRank)
  {
    return srcRank + dstRank;
  }

  void receivedCallback(uint8_t* msg, std::size_t count, int fromRank)
  {
    DEBUG("Received message from process: " << fromRank);

    // first: deserialize
    int* msgs = reinterpret_cast<int*>(msg);
    int msg_count = count / sizeof(int);

    // for all received requests, send the value from the lookup
    for (int i = 0; i < msg_count; ++i)
    {
      // check that the message is as expected
      if(msgs[i] != generate_message(fromRank, my_rank))
      {
        DEBUG("ERROR: message not as expected: " << msgs[i]);
        exit(EXIT_FAILURE);
      }
      else
      {
        // DEBUG("SUCCESS: message received");
        ++msgs_received;
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
        DEBUG("ERROR: LOOKUP message not as expected: " << msgs[i]);
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
        DEBUG("ERROR: ANSWER message not as expected: " << msgs[i]);
        exit(EXIT_FAILURE);
      }
      ++answers_received;
    }
  }

  void test_comm_layer(int repeat_sends=100)
  {
    DEBUG("Testing Comm Layer");
    DEBUG("Size: " << commLayer.getCommSize());
    DEBUG("Rank: " << commLayer.getCommRank());

    // set global rank
    my_rank = commLayer.getCommRank();
    using namespace std::placeholders;
    commLayer.addReceiveCallback(FIRST_TAG, std::bind(&Tester::receivedCallback, this, _1, _2, _3));

    // start sending one message to each:
    for (int l = 0; l < repeat_sends; ++l)
    {
      for (int i = 0; i < commLayer.getCommSize(); ++i)
      {
        int msg = generate_message(my_rank, i);
        commLayer.sendMessage(&msg, sizeof(int), i, FIRST_TAG);
      }
    }

    // call the flush function for this tag
    commLayer.flush(FIRST_TAG);

    // check that all messages have been received
    if (msgs_received != repeat_sends * commLayer.getCommSize())
    {
      DEBUG("ERROR: wrong amount of messages received in phase 1");
      exit(EXIT_FAILURE);
    }


    /* phase 2 communication */


    commLayer.addReceiveCallback(LOOKUP_TAG, std::bind(&Tester::lookup_callback, this, _1, _2, _3));
    commLayer.addReceiveCallback(ANSWER_TAG, std::bind(&Tester::answer_callback, this, _1, _2, _3));

    // sending one message to each:
    for (int l = 0; l < repeat_sends; ++l)
    {
      for (int i = 0; i < commLayer.getCommSize(); ++i)
      {
        int msg = (my_rank+1)*(i+1);
        commLayer.sendMessage(&msg, sizeof(int), i, FIRST_TAG);
      }
    }

    // flush both tags
    commLayer.flush(LOOKUP_TAG);
    commLayer.flush(ANSWER_TAG);

    // check that all messages have been received correctly
    if (answers_received != repeat_sends * commLayer.getCommSize())
    {
      DEBUG("ERROR: wrong amount of messages received in phase 2");
      exit(EXIT_FAILURE);
    }

    DEBUG("This was a triumph.");
    sleep(1);
    DEBUG("I'm making a note here: HUGE SUCCESS.");
    sleep(1);
    DEBUG("It's hard to overstate my satisfaction.");
  }

  Tester(MPI_Comm comm, int comm_size) : commLayer(comm, comm_size) {}

  CommunicationLayer commLayer;
};

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
  Tester tester(comm, p);
  tester.test_comm_layer(100);

  // finalize MPI
  MPI_Finalize();
  return 0;
}

