/**
 * @file   icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.cxx
 * @brief  Utility to dump a artDAQ fragment on screen.
 * @date   April 5, 2021
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h
 * 
 */

// library header
#include "icaruscode/Decode/DecoderTools/Dumpers/FragmentDumper.h"

// ICARUS/SBN libraries
#include "icarusalg/Utilities/BinaryDumpUtils.h" // icarus::ns::util::hexdump()

// C/C++ standard libraries
#include <ostream>
#include <cstddef> // std::size_t
#include <cassert>


// -----------------------------------------------------------------------------
// forward declarations
namespace { void dumpFragmentImpl(std::ostream&, artdaq::Fragment const&); }


// -----------------------------------------------------------------------------
std::ostream& sbndaq::details::operator<<
  (std::ostream& out, DumpFragWrap const& wrap)
  { ::dumpFragmentImpl(out, wrap.frag); return out; }


// -----------------------------------------------------------------------------
sbndaq::details::DumpFragWrap sbndaq::dumpFragment(artdaq::Fragment const& frag)
  { return details::DumpFragWrap{ frag }; }


// -----------------------------------------------------------------------------
namespace {
  
  void dumpFragmentImpl(std::ostream& out, artdaq::Fragment const& frag) {
    
    out << "Fragment summary: " << frag; // ends with <CR>
    
    artdaq::Fragment::byte_t const* headerData = frag.headerBeginBytes();
    std::size_t const headerSize = frag.headerSizeBytes();
    assert(headerData); // this is the pointer to the begin of everything
    artdaq::Fragment::byte_t const* payloadData = frag.dataBeginBytes();
    std::size_t const payloadSize = frag.dataSizeBytes();
    artdaq::Fragment::byte_t const* metaData
      = frag.metadata<artdaq::Fragment::byte_t>();
    std::size_t const metadataSize = metaData? (payloadData - metaData): 0;
    
    out << " - header (" << headerSize << " bytes):";
    if (headerData)
      out << icarus::ns::util::hexdump(headerData, headerSize);
    else out << " n/a" << "\n";
    out << " - metadata (" << metadataSize << " bytes):";
    if (metaData)
      out << icarus::ns::util::hexdump(metaData, metadataSize);
    else out << " n/a" << "\n";
    out << " - data (" << payloadSize << " bytes):";
    if (payloadData)
      out << icarus::ns::util::hexdump(payloadData, payloadSize);
    else out << " n/a" << "\n";
    
  } // dumpFragmentImpl()
  
} // local namespace


// -----------------------------------------------------------------------------
