// Written by Abe White and Jonathan Petersen
// Parallel Computing - CS 5500 / 6500

#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <memory>
#include <vector>

#include <mpi.h>

#include "mpi/DatatypeTags.hpp"

namespace mpi
{
class Communicator
{
public:
  const static MPI_Comm COMMUNICATOR;
  const static int DEFAULT_TAG = 0;
  const static int DEFAULT_ROOT = 0;

  /**
   * Communicator constructor. Sets up MPI and builds a communicator.
   * @param argc Argument count. Should be copied directly from main().
   * @param argv Argument list. Should be copied directly from main().
   */
  Communicator(int argc, char* argv[]);

  /**
   * Copy constructor. Deleted to prevent multiple finalizations.
   * @param other The Communicator to copy from.
   */
  Communicator(const Communicator& other) = delete;

  /**
   * Move Constructor. Default implementation.
   * @param other The Communicator to move from.
   */
  Communicator(Communicator&& other) = default;

  /**
   * Communicator destructor. Used to cleanup MPI.
   */
  ~Communicator();

  /**
   * Copy assignment operator. Deleted to prevent multiple finalizations.
   * @param other The Communicator to copy from.
   * @returns a reference to the new Communicator.
   */
  Communicator& operator=(const Communicator& other) = delete;

  /**
   * Move assignment operator. Default implementation.
   * @param other The Communicator to move from.
   * @returns a reference to the new Communicator.
   */
  Communicator& operator=(Communicator&& other) = default;

  /**
   * Sends a message of any type on this communicator.
   * @param DataType    Type of data to send. Should be auto-deduced.
   * @param data        Data to send.
   * @param destination Where to send the data.
   * @param tag         Tag for the data. Defaults to DEFAULT_TAG.
   */
  template <typename DataType>
  void send(DataType data, int destinaiton, int tag = DEFAULT_TAG) const;

  /**
   * Sends a message of vector type on this communicator. Overloaded function
   * for sending collections of data at one time.
   * @param DataType    Type of data to send. Should be auto-deduced.
   * @param data        Vector of data to send.
   * @param destination Where to send the data.
   * @param tag         Tag for the data. Defaults to DEFAULT_TAG.
   */
  template <typename DataType>
  void send(std::vector<DataType> data,
            int destination,
            int tag = DEFAULT_TAG) const;

  /**
   * Receives one item of any type on this communicator.
   * @param DataType Type of data to receive. Should be auto-deduced.
   * @param source   Where to receive data from. Defaults to any source.
   * @param tag      What tag value to receive. Defaults to DEFAULT_TAG.
   * @returns the received message.
   */
  template <typename DataType>
  DataType recvOne(int source = MPI_ANY_SOURCE, int tag = DEFAULT_TAG) const;

  /**
   * Receives an entire buffer of any type on this communicator.
   * @param DataType Type of data to receive. Should be auto-deduced.
   * @param source   Where to receive data from. Defaults to any source.
   * @param tag      What tag value to receive. Defaults to DEFAULT_TAG.
   * @returns the received vector of data.
   */
  template <typename DataType>
  std::vector<DataType> recv(int source = MPI_ANY_SOURCE,
                             int tag = DEFAULT_TAG) const;

  /**
   * Sends an even number of elements from data to each processor.
   * @param DataType Type of data to receive. Should be auto-deduced.
   * @param data     Array of data to send. data.size() should be a divisor of
   *                 the communicator size, and elements are sent in
   *                 contiguous chunks.
   * @param root     The processor scattering data. Defaults to DEFAULT_ROOT.
   * @param tag      Tag for the data. Defaults to DEFAULT_TAG.
   * @returns the slice of data allocated to this processor.
   */
  template <typename DataType>
  std::vector<DataType> scatter(std::vector<DataType> data,
                                int root = DEFAULT_ROOT,
                                int tag = DEFAULT_TAG) const;

  /**
   * Gathers an even number of elements from each processor to a root processor.
   * @param DataType Type of data to receive. Should be auto-deduced.
   * @param data     Data slice to be sent back to the gathering process.
   * @param root     The processor gathering data. Defaults to DEFAULT_ROOT.
   * @param tag      Tag for the data. Defaults to DEFAULT_TAG.
   * @returns a vector of elements gathered from each processor on the root, and
   * implementation-defined values on non-root processors. Elements are sorted
   * by rank in contiguous chunks.
   */
  template <typename DataType>
  std::vector<DataType> gather(std::vector<DataType> data,
                               int root = DEFAULT_ROOT,
                               int tag = DEFAULT_TAG) const;

  /**
   * Gathers an even number of elements from each processor to all processors.
   * @param DataType Type of data to receive. Should be auto-deduced.
   * @param data     Data slice to be sent back to the gathering process.
   * @param tag      Tag for the data. Defaults to DEFAULT_TAG.
   * @returns a vector of elements gathered from each processor. Elements are
   * sorted by rank in contiguous chunks.
   */
  template <typename DataType>
  std::vector<DataType> allGather(std::vector<DataType> data,
                                  int tag = DEFAULT_TAG) const;

  /**
   * Exchanges data with this processor's n-ary hypercube partner.
   * @param DataType Type of data to exchange. Should be auto-deduced.
   * @param level    What dimension to send on.
   * @param data     Data to exchange.
   * @returns the data from the other processor.
   */
  template <typename DataType>
  DataType cube(unsigned int level, DataType data) const;

  // Member variables
  int m_myRank;
  int m_size;

  // Really hacky destinations, but I'm lazy right now
  int m_nextHighest;
  int m_nextLowest;
};

// Template Definitions
template <typename DataType>
void Communicator::send(DataType data, int destinaiton, int tag) const
{
  MPI_Send(std::make_unique<DataType>(data).get(),
           1,
           mpi::datatype_tags::TAG<DataType>,
           destinaiton,
           tag,
           COMMUNICATOR);
}

template <typename DataType>
void Communicator::send(std::vector<DataType> data,
                        int destination,
                        int tag) const
{
  MPI_Send(data.data(),
           data.size(),
           mpi::datatype_tags::TAG<DataType>,
           destination,
           tag,
           COMMUNICATOR);
}

template <typename DataType>
DataType Communicator::recvOne(int source, int tag) const
{
  auto buffer = std::make_unique<DataType>();
  MPI_Recv(buffer.get(),
           1,
           mpi::datatype_tags::TAG<DataType>,
           source,
           tag,
           COMMUNICATOR,
           MPI_STATUS_IGNORE);
  return *buffer;
}

template <typename DataType>
std::vector<DataType> Communicator::recv(int source, int tag) const
{
  // Figure out how large the message is
  MPI_Status status;
  MPI_Probe(source, tag, COMMUNICATOR, &status);

  int size = 0;
  MPI_Get_count(&status, mpi::datatype_tags::TAG<DataType>, &size);

  // Receive the message
  auto buffer = std::vector<DataType>(size);
  MPI_Recv(buffer.data(),
           buffer.size(),
           mpi::datatype_tags::TAG<DataType>,
           source,
           tag,
           COMMUNICATOR,
           MPI_STATUS_IGNORE);
  return buffer;
}

template <typename DataType>
std::vector<DataType> Communicator::scatter(std::vector<DataType> data,
                                            int root,
                                            int tag) const
{
  typename std::vector<DataType>::size_type sliceSize = data.size() / m_size;
  auto slice = std::vector<DataType>(sliceSize);

  MPI_Scatter(data.data(),
              data.size() / m_size,
              mpi::datatype_tags::TAG<DataType>,
              slice.data(),
              sliceSize,
              mpi::datatype_tags::TAG<DataType>,
              root,
              COMMUNICATOR);

  return slice;
}

template <typename DataType>
std::vector<DataType> Communicator::gather(std::vector<DataType> data,
                                           int root,
                                           int tag) const
{
  typename std::vector<DataType>::size_type totalSize = data.size() * m_size;
  auto total = std::vector<DataType>(totalSize);

  MPI_Gather(data.data(),
             data.size(),
             mpi::datatype_tags::TAG<DataType>,
             total.data(),
             data.size(),
             mpi::datatype_tags::TAG<DataType>,
             root,
             COMMUNICATOR);

  return total;
}

template <typename DataType>
std::vector<DataType> Communicator::allGather(std::vector<DataType> data,
                                              int tag) const
{
  typename std::vector<DataType>::size_type totalSize = data.size() * m_size;
  auto total = std::vector<DataType>(totalSize);

  MPI_Allgather(data.data(),
                data.size(),
                mpi::datatype_tags::TAG<DataType>,
                total.data(),
                data.size(),
                mpi::datatype_tags::TAG<DataType>,
                COMMUNICATOR);

  return total;
}

template <typename DataType>
DataType Communicator::cube(unsigned int level, DataType data) const
{
  auto MASK = 0x1 << level;
  auto partner = (m_myRank ^ MASK) % m_size;

  send(data, partner);
  return recv<DataType>(partner);
}

} // namespace mpi

#endif
