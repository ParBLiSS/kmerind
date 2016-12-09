#ifndef TCLAP_UTILS_HPP_
#define TCLAP_UTILS_HPP_


#if defined(USE_MPI)
#include <mxx/comm.hpp>
#include <tclap/StdOutput.h>
//#include <tclap/CmdLineInterface.h>
//#include <tclap/ArgException.h>

namespace bliss
{
  namespace utils
  {
    namespace tclap
	{
		class MPIOutput : public TCLAP::StdOutput
		{
			protected:
				mxx::comm comm;

			public:
				MPIOutput(mxx::comm const & _comm) : TCLAP::StdOutput(), comm(_comm.copy()) {};

				virtual void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e)
				{
					if (comm.rank() == 0) {
						TCLAP::StdOutput::failure(c, e);
					}
				}

				virtual void usage(TCLAP::CmdLineInterface& c)
				{
					if (comm.rank() == 0) {
						TCLAP::StdOutput::usage(c);
					}

				}

				virtual void version(TCLAP::CmdLineInterface& c)
				{
					if (comm.rank() == 0) {
						TCLAP::StdOutput::version(c);
					}
				}
		};

	} // ns tclap
  } // ns utils

} // ns bliss

#endif

#endif // TCLAP_UTILS_HPP_
