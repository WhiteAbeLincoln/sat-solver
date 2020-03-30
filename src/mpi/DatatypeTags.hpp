// Written by Abe White and Jonathan Petersen
// Parallel Computing - CS 5500 / 6500

#ifndef DATATYPE_TAGS_HPP
#define DATATYPE_TAGS_HPP

#include <cstddef>

#include <mpi.h>

namespace mpi::datatype_tags
{
// Base Variable Template
template <typename DataType> constexpr MPI_Datatype TAG = MPI_DATATYPE_NULL;

// Variable Template Specializations
template <> const MPI_Datatype TAG<int8_t> = MPI_INT8_T;
template <> const MPI_Datatype TAG<int16_t> = MPI_INT16_T;
template <> const MPI_Datatype TAG<int32_t> = MPI_INT32_T;
template <> const MPI_Datatype TAG<int64_t> = MPI_INT64_T;
template <> const MPI_Datatype TAG<uint8_t> = MPI_UINT8_T;
template <> const MPI_Datatype TAG<uint16_t> = MPI_UINT16_T;
template <> const MPI_Datatype TAG<uint32_t> = MPI_UINT32_T;
template <> const MPI_Datatype TAG<uint64_t> = MPI_UINT64_T;
template <> const MPI_Datatype TAG<char> = MPI_CHAR;
template <> const MPI_Datatype TAG<long long int> = MPI_LONG_LONG_INT;
template <> const MPI_Datatype TAG<float> = MPI_FLOAT;
template <> const MPI_Datatype TAG<double> = MPI_DOUBLE;
template <> const MPI_Datatype TAG<long double> = MPI_LONG_DOUBLE;
template <> const MPI_Datatype TAG<wchar_t> = MPI_WCHAR;
template <> const MPI_Datatype TAG<bool> = MPI_C_BOOL;
template <> const MPI_Datatype TAG<unsigned long long> = MPI_UNSIGNED_LONG_LONG;

} // namespace mpi::datatype_tags

#endif
