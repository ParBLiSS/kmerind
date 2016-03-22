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

/*
 * @file      unix_domain_socket.h
 * @ingroup   io
 * @author  Tony Pan <tpan7@gatech.edu>
 *
 * @brief     helper function to replicate file descriptor across all processes on the same node.
 * @details   based on information from
 *            http://www.microhowto.info/howto/listen_for_and_receive_udp_datagrams_in_c.html
 *            http://stackoverflow.com/questions/27014955/socket-connect-vs-bind
 *            http://beej.us/guide/bgipc/output/html/multipage/unixsock.html
 *            http://www.thomasstover.com/uds.html
 *            http://stackoverflow.com/questions/28003921/sending-file-descriptor-by-linux-socket
 *            http://blog.varunajayasiri.com/passing-file-descriptors-between-processes-using-sendmsg-and-recvmsg
 *
 *
 *  Created on: Mar 20, 2016
 */

#ifndef UNIX_DOMAIN_SOCKET_H_
#define UNIX_DOMAIN_SOCKET_H_


#include <cstdio>
#include <sys/socket.h>
#include <sys/un.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>

#include <sys/stat.h>
#include <fcntl.h>

#include <mxx/comm.hpp>

namespace bliss {

  namespace io {

    namespace util {


      /// serialize a unix domain socket message as a file descriptor
      static int send_fd(int const & socket, int const & fd_to_send)
      {
        struct msghdr socket_message;
        struct iovec io_vector[1];
        struct cmsghdr *control_message = NULL;
        char message_buffer[1];
        /* storage space needed for an ancillary element with a paylod of length is CMSG_SPACE(sizeof(length)) */
        char ancillary_element_buffer[CMSG_SPACE(sizeof(int))];
        int available_ancillary_element_buffer_space;

        /* at least one vector of one byte must be sent */
        message_buffer[0] = 'F';
        io_vector[0].iov_base = message_buffer;
        io_vector[0].iov_len = 1;

        /* initialize socket message */
        memset(&socket_message, 0, sizeof(struct msghdr));

        //  socket_message.msg_name = &target;
        //  socket_message.msg_namelen = sizeof(target.sun_family) + strlen(target.sun_path);

        socket_message.msg_iov = io_vector;
        socket_message.msg_iovlen = 1;

        /* provide space for the ancillary data */
        available_ancillary_element_buffer_space = CMSG_SPACE(sizeof(int));
        memset(ancillary_element_buffer, 0, available_ancillary_element_buffer_space);
        socket_message.msg_control = ancillary_element_buffer;
        socket_message.msg_controllen = available_ancillary_element_buffer_space;

        /* initialize a single ancillary data element for fd passing */
        control_message = CMSG_FIRSTHDR(&socket_message);
        control_message->cmsg_level = SOL_SOCKET;
        control_message->cmsg_type = SCM_RIGHTS;
        control_message->cmsg_len = CMSG_LEN(sizeof(int));
        *((int *) CMSG_DATA(control_message)) = fd_to_send;

        return sendmsg(socket, &socket_message, 0);
      }

      /// deserialize a unix domain socket message as a file descriptor
      static int recv_fd(int const & socket)
      {
        int sent_fd;
        struct msghdr socket_message;
        struct iovec io_vector[1];
        struct cmsghdr *control_message = NULL;
        char message_buffer[1];
        char ancillary_element_buffer[CMSG_SPACE(sizeof(int))];

        /* start clean */
        memset(&socket_message, 0, sizeof(struct msghdr));
        memset(ancillary_element_buffer, 0, CMSG_SPACE(sizeof(int)));

        /* setup a place to fill in message contents */
        //  socket_message.msg_name = &source;
        //  socket_message.msg_namelen = sizeof(source.sun_family) + strlen(source.sun_path);
        io_vector[0].iov_base = message_buffer;
        io_vector[0].iov_len = 1;
        socket_message.msg_iov = io_vector;
        socket_message.msg_iovlen = 1;

        /* provide space for the ancillary data */
        socket_message.msg_control = ancillary_element_buffer;
        socket_message.msg_controllen = CMSG_SPACE(sizeof(int));

        if(recvmsg(socket, &socket_message, MSG_CMSG_CLOEXEC) < 0)
          return -1;

        if(message_buffer[0] != 'F')
        {
          /* this did not originate from the above function */
          return -1;
        }

        if((socket_message.msg_flags & MSG_CTRUNC) == MSG_CTRUNC)
        {
          /* we did not provide enough space for the ancillary element array */
          return -1;
        }

        /* iterate ancillary elements */
        for(control_message = CMSG_FIRSTHDR(&socket_message);
            control_message != NULL;
            control_message = CMSG_NXTHDR(&socket_message, control_message))
        {
          if( (control_message->cmsg_level == SOL_SOCKET) &&
              (control_message->cmsg_type == SCM_RIGHTS) )
          {
            sent_fd = *((int *) CMSG_DATA(control_message));
            return sent_fd;
          }
        }

        return -1;
      }

      /// create an address for unix domain socket
      static struct sockaddr_un make_address(::std::string const & name,  socklen_t &address_length) {
          struct sockaddr_un address;

          memset(&address, 0, sizeof(address));
          address.sun_family = AF_UNIX;
          strcpy(address.sun_path, name.c_str());

          address_length = sizeof(address.sun_family) + strlen(address.sun_path);
          address.sun_path[0] = 0;

          return address;
      }

      ///  create a unix domain socket
      static int create_socket() {
        int socket_fd;

        if((socket_fd = socket(AF_UNIX, SOCK_DGRAM, 0)) < 0)
        {
          perror("server: socket");
          throw ::std::ios_base::failure("ERROR: unable to create socket");
        }

        return socket_fd;
      }

      /// bind the socket for receiving an fd
      static void bind_receiver(int const & socket_fd, struct sockaddr_un const & address, socklen_t const & address_length) {

        unlink(address.sun_path);
        if (bind(socket_fd, (const struct sockaddr *) &address, address_length) < 0)
            {
          close(socket_fd);
          perror("server: bind");
          throw ::std::ios_base::failure("ERROR: unable to bind socket to server address");
            }
      }

      /// close the socket that receives an fd
      static void shutdown_receiver(int const & socket_fd, struct sockaddr_un & address) {
        unlink(address.sun_path);
        close(socket_fd);
      }

      /// send file descriptor from src rank to target_rank
      static void send_file_descriptor(int & fd, int const & id,
    		  int const & nprocs, int const & src_rank, int const & target_rank) {
        assert(src_rank < target_rank);

        if (target_rank >= nprocs) return;  // nothing to do

        // get the file name for the target
        std::stringstream ss;
        ss << "#FD_W" << id << "_H" << target_rank;
        std::string target_name = ss.str();
        socklen_t address_length;


        // create the socket
        int socket_fd = create_socket();

        // create the address
        struct sockaddr_un remote = make_address( target_name, address_length);

        int sent = 0;
        int result1 = -1;
        int result2 = -1;

        // connect and send.
        do {
          // loop until connection is accepted.
          if ((result1 = connect(socket_fd,
                                 (const struct sockaddr *) &remote,
                                 address_length)) < 0) {
            continue;
          }

          // send fd.
          if ( (result2 = send_fd(socket_fd, fd) ) < 0) {
            perror("client: send_fd");
            usleep(1000);
            continue;
          };

          sent = 1;
          //printf("successful send fd %d from %d to %d\n", fd, src_rank, target_rank);

        } while (sent == 0);

        close(socket_fd);
      }

      /// send file descriptor from src rank to target_rank
      static void recv_file_descriptor(int & fd, int const & id,
    		  int const & nprocs, int const & src_rank, int const & target_rank) {
        assert(src_rank < target_rank);

        if (src_rank < 0) return;  // nothing to do

        // get the file name for the target
        std::stringstream ss;
        ss << "#FD_W" << id << "_H" << target_rank;
        std::string target_name = ss.str();
        socklen_t address_length;

        // create the socket
        int socket_fd = create_socket();

        // create the address
        struct sockaddr_un local = make_address( target_name, address_length);

        // start the server
        bind_receiver(socket_fd, local, address_length);

        // wait for receiving fd.
        fd = recv_fd(socket_fd);

        //printf("successful recved fd %d from %d to %d\n", fd, src_rank, target_rank);

        // shutdown.
        shutdown_receiver(socket_fd, local);

      }



      /// broadcast fd from rank 0 to all other processes on the same node, log iterations.
      static void broadcast_file_descriptor(int & fd, int const & id,  int const & nprocs, int const & rank) {

        // first get the smallest power of 2 >= p;
        int step = 1;
        while (step < nprocs) step <<= 1;
        step >>= 1;  // step is half of above

        for (; step > 0; step >>= 1 ) {
          if ((rank % (step * 2)) == 0) {
            // this is source
            send_file_descriptor(fd, id, nprocs, rank, rank + step);
          } else if ((rank % (step * 2)) == step) {
            // this is target.
            recv_file_descriptor(fd, id, nprocs, rank - step, rank);
          }
        }


      }


    } // ns util

  } // ns io

} // ns bliss



#endif /* UNIX_DOMAIN_SOCKET_H_ */
