/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tools for indexing arbitrary arrays in a variety of ways.
#pragma once


namespace ATS {
namespace Utils {

enum class Indexer_kind { SCALAR = 0, ELM, ATS };

enum class Domain_kind { SURF, SUBSURF };


template <Indexer_kind ikind, Domain_kind dkind, typename T>
struct Indexer;

//
// Indexer for scalars.
//
template <Domain_kind dkind, typename T>
struct Indexer<Indexer_kind::SCALAR, dkind, T> {
  template <typename D = T>
  typename std::enable_if<dkind == Domain_kind::SUBSURF, D&>::type
  get(T const* const& val, const int col, const int cic) const
  {
    return *val_;
  }

  template <typename D = T>
  typename std::enable_if<dkind == Domain_kind::SURF, D&>::type
  get(T const* const& val, const int i) const
  {
    return *val_;
  }

 private:
  T* val_;
};


//
// Indexer for ELM
//
template <Domain_kind dkind, typename T>
struct Indexer<Indexer_kind::ELM, dkind, T> {
  Indexer(int ncells_per_col = -1) : ncells_per_col_(ncells_per_col) {}

  // for use by subsurface
  template <typename D = T>
  typename std::enable_if<dkind == Domain_kind::SUBSURF, D&>::type
  get(T* const& val, const int col, const int cic) const
  {
    return val[col * ncells_per_col_ + cic];
  }

  // for use by surface
  template <typename D = T>
  typename std::enable_if<dkind == Domain_kind::SURF, D&>::type
  get(T* const& val, const int i) const
  {
    return val[i];
  }

 private:
  int ncells_per_col_;
};


//
// Indexer for ATS
//
template <Domain_kind dkind, typename T>
struct Indexer<Indexer_kind::ATS, dkind, T> {
  Indexer(const Amanzi::AmanziMesh::Mesh& mesh) : mesh_(mesh) {}

  // for use by subsurface
  template <typename D = T>
  typename std::enable_if<dkind == Domain_kind::SUBSURF, T&>::type
  get(T* const& val, const int col, const int cic) const
  {
    return val[mesh_.cells_of_column(col)[cic]];
  }

  // for use by surface
  template <typename D = T>
  typename std::enable_if<dkind == Domain_kind::SURF, D&>::type
  get(T* const& val, const int i) const
  {
    return val[i];
  }

 private:
  const Amanzi::AmanziMesh::Mesh& mesh_;
};


} // namespace Utils
} // namespace ATS
