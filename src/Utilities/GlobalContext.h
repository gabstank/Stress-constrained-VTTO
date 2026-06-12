#ifndef GLOBALCONTEXT_H
#define GLOBALCONTEXT_H

// C++ headers
#include <iostream>

// Deal.II headers
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

using namespace dealii;

/**
 * @class GlobalContext
 * @brief Singleton class to manage global MPI context and utilities.
 *
 * The GlobalContext class provides a singleton instance to manage the MPI
 * communicator, process information, and utilities such as conditional output
 * and timing.
 */

class GlobalContext {
public:
  /**
   * @brief Get the singleton instance of GlobalContext.
   *
   * This function returns the singleton instance of the GlobalContext class.
   *
   * @return GlobalContext& Reference to the singleton instance.
   */
  static GlobalContext &get_instance() {
    static GlobalContext instance;
    return instance;
  }

  /**
   * @var MPI_Comm mpi_communicator
   * @brief The MPI communicator used for parallel communication.
   */

  MPI_Comm mpi_communicator;

  /**
   * @var unsigned int this_mpi_process
   * @brief The rank of the current MPI process.
   */

  unsigned int this_mpi_process;

  /**
   * @var unsigned int n_mpi_processes
   * @brief The total number of MPI processes.
   */

  unsigned int n_mpi_processes;

  /**
   * @var ConditionalOStream pcout
   * @brief Conditional output stream for parallel output.
   *
   * This stream outputs to std::cout only if the current process is the root process.
   */

  ConditionalOStream pcout;

  /**
   * @var TimerOutput compute_timer
   * @brief Timer for measuring computation time.
   *
   * This timer outputs timing information to the conditional output stream.
   */
  TimerOutput compute_timer;

  unsigned int opt_iteration;

private:
  /**
   * @brief Private constructor to initialize the GlobalContext instance.
   *
   * The constructor initializes the MPI communicator, process information,
   * conditional output stream, and computation timer.
   */

  GlobalContext()
      : mpi_communicator(MPI_COMM_WORLD),
        this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)),
        n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)),
        pcout(std::cout, this_mpi_process == 0),
        compute_timer(pcout, TimerOutput::summary, TimerOutput::wall_times), opt_iteration(0) {}

  /**
   * @brief Deleted copy constructor.
   *
   * The copy constructor is deleted to ensure the singleton property.
   */

  GlobalContext(const GlobalContext &) = delete;

  /**
   * @brief Deleted assignment operator.
   *
   * The assignment operator is deleted to ensure the singleton property.
   */
  void operator=(const GlobalContext &) = delete;
};

/**
 * @brief Retrieves the MPI process identifier for the current process.
 *
 * This function returns a reference to the MPI process identifier for the
 * current process. The identifier is managed by the GlobalContext singleton
 * instance.
 *
 * @return A reference to an unsigned int representing the MPI process identifier.
 */
inline MPI_Comm &mpi_communicator() { return GlobalContext::get_instance().mpi_communicator; }

/**
 * @brief Retrieves the number of MPI processes in the communicator.
 *
 * This function returns a reference to the number of MPI processes in the
 * communicator. The number of processes is managed by the GlobalContext
 * singleton instance.
 *
 * @return A reference to an unsigned int representing the number of MPI processes.
 */
inline unsigned int &this_mpi_process() { return GlobalContext::get_instance().this_mpi_process; }

/**
 * @brief Retrieves the number of MPI processes in the communicator.
 *
 * This function returns a reference to the number of MPI processes in the
 * communicator. The number of processes is managed by the GlobalContext
 * singleton instance.
 *
 * @return A reference to an unsigned int representing the number of MPI processes.
 */
inline unsigned int &n_mpi_processes() { return GlobalContext::get_instance().n_mpi_processes; }

/**
 * @brief Retrieves the conditional output stream for parallel output.
 *
 * This function returns a reference to the conditional output stream for parallel
 * output. The stream is managed by the GlobalContext singleton instance.
 *
 * @return A reference to a ConditionalOStream object for parallel output.
 */
inline ConditionalOStream &pcout() { return GlobalContext::get_instance().pcout; }

/**
 * @brief Retrieves the timer for measuring computation time.
 *
 * This function returns a reference to the timer for measuring computation time.
 * The timer is managed by the GlobalContext singleton instance.
 *
 * @return A reference to a TimerOutput object for measuring computation time.
 */
inline TimerOutput &compute_timer() { return GlobalContext::get_instance().compute_timer; }

inline unsigned int &opt_iteration() { return GlobalContext::get_instance().opt_iteration; }

#endif // GLOBALCONTEXT_H