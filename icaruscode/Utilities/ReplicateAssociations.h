/**
 * @file   icaruscode/Utilities/ReplicateAssociations.h
 * @brief  Utilities to "modify" associations.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   September 2, 2022
 *
 * This is a header-only library
 */

#ifndef ICARUSCODE_UTILITIES_REPLICATEASSOCIATIONS_H
#define ICARUSCODE_UTILITIES_REPLICATEASSOCIATIONS_H


// framework libraries
#include "canvas/Persistency/Common/Assns.h"

// C/C++ standard libraries
#include <utility> // std::declval()
#include <type_traits> // std::is_void_v, ...


//------------------------------------------------------------------------------
namespace util {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Returns a copy of `assns` with the first element replaced.
   * @tparam Assns type of the associations to be processed
   * @tparam LeftMap type of mapping to apply to the left side of the association
   * @param assns the associations to be reprocessed
   * @param leftMap the mapping to apply to the left side of the association
   * @return a copy of `assns` with a modified left side
   * 
   * This function creates a new `art::Assns` association with the left pointer
   * `ptr` of each of the links replaced by its correspondent `leftMap(ptr)`.
   * 
   * In general, `leftMap` should support the mapping of all the left pointers
   * in `assns`. In practice, the `leftMap` object may decide how to deal with
   * any unknown pointers.
   * 
   * Note that the type of associations returned by this function is determined
   * by the type returned by the mapping (if for example the mapping operates
   * on `art::Ptr<A>` and returns `art::Ptr<B>`, an input `art::Assns<A, C, D>`
   * will be turned ito a `art::Assns<B, C, D>`). The mapping must return an
   * _art_ pointer to some object.
   */
  template <typename Assns, typename LeftMap>
  auto replaceAssnsFirst(Assns const& assns, LeftMap const& leftMap);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace util


//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
namespace util::details {
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // something like this is in lardata
  template <typename T> struct AssnsTraits;
  
  template <typename Left, typename Right, typename Data>
  struct AssnsTraits<art::Assns<Left, Right, Data>> {
    using left_t = Left;
    using right_t = Right;
    using data_t = Data;
    static constexpr bool hasMetadata = !std::is_void_v<data_t>;
  };
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  /// Whether `T` is an instance of `art::Ptr`.
  template <typename T> struct isArtPtr: std::false_type {};
  
  template <typename T>
  struct isArtPtr<art::Ptr<T>>: std::true_type {};
  
  template <typename T>
  constexpr bool isArtPtr_v = isArtPtr<T>::value;
  
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
} // util::details


//------------------------------------------------------------------------------
template <typename Assns, typename LeftMap>
auto util::replaceAssnsFirst(Assns const& assns, LeftMap const& leftMap) {
  
  /*
   * Some abundant metaprogramming to detect the type of the type out of the
   * new pointer output by `leftMap`.
   */
  using AssnsTraits_t = details::AssnsTraits<Assns>;
  
  using left_t = typename AssnsTraits_t::left_t;
  using right_t = typename AssnsTraits_t::right_t;
  using data_t = typename AssnsTraits_t::data_t;
  constexpr bool hasMetadata = AssnsTraits_t::hasMetadata;
  
  using destLeftPtr_t
    = std::decay_t<decltype(leftMap(std::declval<art::Ptr<left_t>>()))>;
  static_assert(details::isArtPtr_v<destLeftPtr_t>);
  using destLeft_t = typename destLeftPtr_t::value_type;
  using destAssns_t = art::Assns<destLeft_t, right_t, data_t>;

  // 
  destAssns_t replicaAssns;
  for (auto const& assn: assns) {
    
    destLeftPtr_t destLeft = leftMap(assn.first);
    
    if constexpr(hasMetadata) {
      replicaAssns.addSingle(std::move(destLeft), assn.second, *(assn.data));
    }
    else {
      replicaAssns.addSingle(std::move(destLeft), assn.second);
    }
    
  } // for association link
  
  return replicaAssns;
} // util::replaceAssnsFirst()


//------------------------------------------------------------------------------

#endif // ICARUSCODE_UTILITIES_REPLICATEASSOCIATIONS_H
