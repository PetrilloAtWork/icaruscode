/**
 * @file   icaruscode/Utilities/DumpProvenanceTree_module.cc
 * @brief  Prints in a tree-like structure the provenance of the data products.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July 28, 2022
 * 
 */

// LArSoft libraries
#include "larcorealg/CoreUtils/values.h" // util::values(), util::const_values()
#include "larcorealg/CoreUtils/get_elements.h" // util::get_const_elements()

// framework libraries
#include "art/Framework/Modules/ProvenanceDumper.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/ConfigurationTable.h" // fhiclcpp::WrappedTable
#include "fhiclcpp/types/TableFragment.h"
#include "fhiclcpp/types/Atom.h"

// C/C++ standard libraries
#include <algorithm> // std::sort(), std::lower_bound()
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <string>


// -----------------------------------------------------------------------------
namespace details {
  
  // ---------------------------------------------------------------------------
  template <typename T>
  std::string quantityStr
    (T q, std::string const& unit, std::string const& units = "");
  
  // ---------------------------------------------------------------------------
  
  /// This class follows the interface and prescriptions required by
  /// `art::ProvenanceDumper<>`
  struct DumpProvenanceTreeDetails {
    
    /// Configuration fragment for this algorithm.
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<bool> OmitRedundantProcessNames {
        Name{ "omitRedundantProcessNames"},
        Comment{ "print process names only with ambiguous module labels" },
        true
        };
      
      fhicl::Atom<bool> OmitPreviousModules {
        Name{ "omitPreviousModules"},
        Comment{ "prints the details of modules only the first time" },
        true
        };
      
      fhicl::Atom<bool> PrintOnlyLastGeneration {
        Name{ "printOnlyLastGeneration"},
        Comment{ "prints the tree only for modules with no descendants" },
        true
        };
      
      fhicl::Atom<std::string> LogCategory {
        Name{ "logCategory" },
        Comment{ "message facility stream name for this module output" },
        "DumpProvenanceTree"
        };
      
    }; // Config
    
    
    /// Constructor: gets the configuration parameters from the module.
    DumpProvenanceTreeDetails(fhicl::TableFragment<Config> const& params);
    
    
    // --- BEGIN -- Callbacks --------------------------------------------------
    
    /// Called at the end of each fiscal year.
    void beginJob();
    
    /// Called for each provenance entry in the run.
    void processRunProvenance(art::Provenance const& provenance);
    
    /// Called after all provenance entries in the run are processed
    /// (not called if there were none)
    void postProcessRun();
    
    /// Called for each provenance entry in the subrun.
    void processSubRunProvenance(art::Provenance const& provenance);
    
    /// Called after all provenance entries in the subrun are processed
    /// (not called if there were none)
    void postProcessSubRun();
    
    /// Called for each provenance entry in the event.
    void processEventProvenance(art::Provenance const& provenance);
    
    /// Called after all provenance entries in the event are processed
    /// (not called if there were none)
    void postProcessEvent();
    
    /// Called on the night of full moon.
    void EndJob();
    
    // --- END ---- Callbacks --------------------------------------------------
    
      private:
    
    // --- BEGIN -- Configuration ----------------------------------------------
    
    struct PrintOptions_t {
      std::string indentString = "    ";
      bool omitRedundantProcessNames = true;
      bool omitPreviousModules = true;
      bool printOnlyLastGeneration = true;
    };
    
    PrintOptions_t const fPrintOptions;
    std::string const fLogCategory;
    
    // --- END ---- Configuration ----------------------------------------------
    
    
    struct ParentageInfo {
      
      art::ProductID ID; ///< ID of this product.
      
      std::vector<art::ProductID> parentIDs; ///< ID of all the parents.
      
      std::string className; ///< The class of this data product.
      art::InputTag tag; ///< The input tag of this data product.
      
      std::vector<art::ProductID> descendantIDs; ///< ID of all descendants;
      
      
      bool operator< (ParentageInfo const& other) const
        { return ID < other.ID; }
      
    }; // struct ParentageInfo
    
    /// Summary information of the modules.
    struct ModuleInfo {
      
      art::InputTag tag; ///< The input tag of this data product.
      bool unique = false; ///< If there is no other module with the same label.
      
      std::vector<ParentageInfo const*> products; ///< Products of the module.
      
      /// ID of all parent products.
      std::set<art::ProductID> inputIDs;
      
      /// ID of all descendant products.
      std::set<art::ProductID> descendantIDs;
      
      /// Tag of all modules producing the parent products.
      std::set<art::InputTag> inputModules;
      
      /// Tag of all modules producing the descendant products.
      std::set<art::InputTag> descendantModules;
      
      art::InputTag simpleTag() const
        { return unique? art::InputTag{ tag.label() }: tag; }
      
    }; // ModuleInfo
    
    
    using AllParentageInfo_t = std::map<art::ProductID, ParentageInfo>;
    using AllModuleInfo_t = std::map<art::InputTag, ModuleInfo>;
    
    /// All parent info for current event.
    AllParentageInfo_t fParentInfo;
    
    
    /// Helper class to print a parentage tree.
    template <typename Stream>
    class ModulePrinter {
      
      Stream& fOut; ///< Output stream.
      AllModuleInfo_t const& fParentInfo; ///< All parentage info.
      PrintOptions_t const& fOptions; ///< Printing options.
      
      std::set<art::InputTag> fPrinted; ///< List of printed modules.
      
      /// Returns a pointer to the specified product, `nullptr` if not found.
      ModuleInfo const* findModule(art::InputTag const& tag) const;
      
      /// Starts a new line with its indentation.
      Stream& startLine(unsigned int indentLevel) const;
      
      /// Prints a "friendly" module name.
      void printModuleName
        (art::InputTag const& tag, ModuleInfo const* info) const;
      /// Prints a "friendly" module name (and looks up the info).
      void printModuleName(art::InputTag const& tag) const;
      
      /// Prints a module and its ancestors, and records this printing.
      void printModule
        (art::InputTag const& tag, ModuleInfo const* info, unsigned int level);
      void printModule(art::InputTag const& tag, unsigned int level)
        { printModule(tag, findModule(tag), level); }
      
        public:
      
      ModulePrinter(
        Stream& out, AllModuleInfo_t const& parentInfo,
        PrintOptions_t const& options
        );

      /// Prints the hierarchy of modules starting with `info`.
      void operator() (ModuleInfo const& info, unsigned int level = 0U)
        { printModule(info.tag, &info, level); }
      
      /// Prints the hierarchy of products starting with `ID`.
      void operator() (art::InputTag const& tag, unsigned int level = 0U)
        { printModule(tag, level); }
      
    }; // ModulePrinter<>
    
    
    /// Fills the descendant ID of all the products in `parentInfo`.
    void fillDescendantInformation(AllParentageInfo_t& parentInfo) const;
    
    /// Creates a map of all modules and connections from product information.
    AllModuleInfo_t buildModuleMap(AllParentageInfo_t const& parentInfo) const;

  }; // DumpProvenanceTreeDetails
  
  
  // ---------------------------------------------------------------------------

} // namespace details


// -----------------------------------------------------------------------------
namespace sbn { class DumpProvenanceTree; }
/**
 * @brief Prints on screen the hierarchy of data products.
 * 
 * 
 * 
 * Configuration options
 * ======================
 * 
 * * `OmitRedundantProcessNames` (flag, default: `true`): if set, process names
 *    are printed only when there are different modules with the same label.
 *    If unset, process names are always printed.
 * * `OmitPreviousModules` (flag, default: `true`): if set, modules that have
 *    previously appeared in the printout are not expanded again, instead they
 *    are just mentioned in a list.
 *    If unset, all modules are expanded each time.
 * * `PrintOnlyLastGeneration` (flag, default: `true`): if set, only modules
 *    which have no descendants are printed as root of a tree; the others will
 *    only appear in the expanded trees of their descendants. If unset, each
 *    module is also printed as a tree root.
 * * `LogCategory` (string, default: `DumpProvenanceTree`): name of the
 *   message facility category used to send the module output to screen.
 * 
 * 
 */
class sbn::DumpProvenanceTree
  : public art::ProvenanceDumper<details::DumpProvenanceTreeDetails>
{
  
  using DumperBase_t = art::ProvenanceDumper<details::DumpProvenanceTreeDetails>;
  
    public:
  
  using art::ProvenanceDumper<details::DumpProvenanceTreeDetails>::ProvenanceDumper;
  
  DumpProvenanceTree(DumperBase_t::Parameters const& params)
    : art::ProvenanceDumper<details::DumpProvenanceTreeDetails>{ params }
    {}

  // is there anything to be defined here?
  
}; // sbn::DumpProvenanceTree


// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
template <typename T>
std::string details::quantityStr
  (T q, std::string const& unit, std::string const& units /* = "" */)
{
  std::string s = std::to_string(q);
  if (q == 1) {
    if (!unit.empty()) { s += ' '; s += unit; }
  }
  else { // plural
    if (!units.empty() || !unit.empty())
      { s += ' '; s += (units.empty())? unit + "s": units; }
  }
  return s;
} // quantityStr()


// -----------------------------------------------------------------------------
// ---  details::DumpProvenanceTreeDetails
// -----------------------------------------------------------------------------
template <typename Stream>
auto details::DumpProvenanceTreeDetails::ModulePrinter<Stream>::findModule
  (art::InputTag const& tag) const -> ModuleInfo const*
{
  auto const it = fParentInfo.find(tag);
  return (it == fParentInfo.end())? nullptr: &(it->second);
}


// -----------------------------------------------------------------------------
template <typename Stream>
Stream& details::DumpProvenanceTreeDetails::ModulePrinter<Stream>::startLine
  (unsigned int indentLevel) const
{
  fOut << "\n";
  while (indentLevel-- > 0) fOut << fOptions.indentString;
  return fOut;
}


// -----------------------------------------------------------------------------
template <typename Stream>
void details::DumpProvenanceTreeDetails::ModulePrinter<Stream>::printModuleName
  (art::InputTag const& tag, ModuleInfo const* info) const
{
  fOut << "'" <<
    ((info && fOptions.omitRedundantProcessNames)? info->simpleTag(): tag)
      .encode()
    << "'";
  
}


// -----------------------------------------------------------------------------
template <typename Stream>
void details::DumpProvenanceTreeDetails::ModulePrinter<Stream>::printModuleName
  (art::InputTag const& tag) const
  { printModuleName(tag, findModule(tag)); }


// -----------------------------------------------------------------------------
template <typename Stream>
void details::DumpProvenanceTreeDetails::ModulePrinter<Stream>::printModule
  (art::InputTag const& tag, ModuleInfo const* info, unsigned int level)
{
  static const char* arrow { "`-<-" };
  fPrinted.insert(tag);
  startLine(level);
  if (level > 0) fOut << arrow;
  printModuleName(tag, info);
  if (!info) { fOut << " <no information available>"; return; }
  fOut << " (" << quantityStr(info->products.size(), "product") << ")";
  if (!info->inputModules.empty())
    fOut << " with " << quantityStr(info->inputModules.size(), "parent");
  std::vector<art::InputTag> alreadyPrinted;
  for (art::InputTag const& inputTag: info->inputModules) {
    if (fOptions.omitPreviousModules && fPrinted.count(inputTag))
      alreadyPrinted.push_back(inputTag);
    else printModule(inputTag, level + 1);
  }
  if (!alreadyPrinted.empty()) {
    startLine(level + 1);
    if (alreadyPrinted.size() != info->inputModules.size())
      fOut << "... plus ";
    else
      fOut << arrow;
    auto itParentModule = alreadyPrinted.begin();
    auto const pend = alreadyPrinted.end();
    printModuleName(*itParentModule);
    while (++itParentModule != pend) {
      fOut << ", ";
      printModuleName(*itParentModule);
    }
    fOut << " (see above)";
  } // if some already printed
} // printModule()


// -----------------------------------------------------------------------------
template <typename Stream>
details::DumpProvenanceTreeDetails::ModulePrinter<Stream>::ModulePrinter(
  Stream& out, AllModuleInfo_t const& parentInfo,
  PrintOptions_t const& options
)
  : fOut{ out }
  , fParentInfo{ parentInfo }
  , fOptions{ options }
{}


// -----------------------------------------------------------------------------
details::DumpProvenanceTreeDetails::DumpProvenanceTreeDetails
  (fhicl::TableFragment<Config> const& params)
  : fPrintOptions{
        "    "                                // indentString
      , params().OmitRedundantProcessNames()  // omitRedundantProcessNames
      , params().OmitPreviousModules()        // omitPreviousModules
      , params().PrintOnlyLastGeneration()    // printOnlyLastGeneration
      }
  , fLogCategory{ params().LogCategory() }
{}


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::beginJob() {}


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::processRunProvenance
  (art::Provenance const& provenance)
{
}


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::postProcessRun() {}


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::processSubRunProvenance
  (art::Provenance const& provenance)
{}


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::postProcessSubRun() {}


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::processEventProvenance
  (art::Provenance const& provenance)
{
  
  // collect
  ParentageInfo info;
  info.ID = provenance.productID();
  info.parentIDs = provenance.parents();
  info.className = provenance.friendlyClassName(); // not that friendly...
  info.tag = provenance.inputTag();
  // may also collect flag information: valid, present, produced.
  
  fParentInfo.emplace(provenance.productID(), std::move(info));
  
} // details::DumpProvenanceTreeDetails::processEventProvenance()


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::postProcessEvent() {
  
  /*
   * The plan:
   *  1. build the relations between data products
   *  2. surmise the list of modules involved
   *  3. rebuild the relations between modules from their data products
   *  4. print in a tree structure the modules with no descendants
   */
  
  fillDescendantInformation(fParentInfo);
  
  AllModuleInfo_t const moduleMap = buildModuleMap(fParentInfo);
  
  mf::LogInfo out { fLogCategory };
  out << "Collected information for "
    << quantityStr(fParentInfo.size(), "data product")
    << " in " << quantityStr(moduleMap.size(), "module");
  
  ModulePrinter print{ out, moduleMap, fPrintOptions };
  for (ModuleInfo const& info: util::const_values(moduleMap)) {
    if (fPrintOptions.printOnlyLastGeneration && !info.descendantModules.empty())
      continue;
    
    print(info);
    
  } // for
  
  
} // details::DumpProvenanceTreeDetails::postProcessEvent()


// -----------------------------------------------------------------------------
void details::DumpProvenanceTreeDetails::EndJob() {}


// -----------------------------------------------------------------------------
/// Fills the descendant ID of all the products in `parentInfo`.
void details::DumpProvenanceTreeDetails::fillDescendantInformation
  (AllParentageInfo_t& parentInfo) const
{
  for (auto const& [ ID, info ]: parentInfo) {
    for (art::ProductID const& parentID: info.parentIDs) {
      auto itParent = parentInfo.find(parentID);
      if (itParent != parentInfo.end())
        itParent->second.descendantIDs.push_back(ID);
    } // for all parents
  } // for all products
} // details::DumpProvenanceTreeDetails::fillDescendantInformation()


// -----------------------------------------------------------------------------
/// Creates a map of all modules and connections from product information.
auto details::DumpProvenanceTreeDetails::buildModuleMap
  (AllParentageInfo_t const& parentInfo) const -> AllModuleInfo_t
{
  
  // turns a product label in a module label (i.e. without instance name)
  auto const moduleOf = [](art::InputTag const& tag)
    { return art::InputTag{ tag.label(), "", tag.process() }; };
  
  //
  // first pass: collect all modules
  //
  AllModuleInfo_t moduleMap;
  for (ParentageInfo const& info: util::const_values(parentInfo)) {
    art::InputTag const moduleTag = moduleOf(info.tag);
    ModuleInfo& moduleInfo = moduleMap[moduleTag];
    if (moduleInfo.tag.empty())
      moduleInfo.tag = moduleTag; // never seen before: initialize
    
    moduleInfo.products.push_back(&info);
    
    moduleInfo.inputIDs.insert
      (info.parentIDs.cbegin(), info.parentIDs.cend());
    moduleInfo.descendantIDs.insert
      (info.descendantIDs.cbegin(), info.descendantIDs.cend());
    
  } // for
  
  //
  // label count
  //
  std::unordered_map<std::string, unsigned int> moduleLabelCounts;
  for (art::InputTag const& tag: util::get_const_elements<0U>(moduleMap))
    ++moduleLabelCounts[tag.label()];
  
  //
  // resolve the module names in the relation graph
  //
  auto findProduct
    = [&parentInfo](art::ProductID const& ID) -> ParentageInfo const*
    {
      auto const it = parentInfo.find(ID);
      return (it == parentInfo.end())? nullptr: &(it->second);
    };
  
  for (ModuleInfo& moduleInfo: util::values(moduleMap)) {
    
    for (art::ProductID const& ID: moduleInfo.inputIDs) {
      if (ParentageInfo const* prodInfo = findProduct(ID))
        moduleInfo.inputModules.insert(moduleOf(prodInfo->tag));
    }
    
    for (art::ProductID const& ID: moduleInfo.descendantIDs) {
      if (ParentageInfo const* prodInfo = findProduct(ID))
        moduleInfo.descendantModules.insert(moduleOf(prodInfo->tag));
    }
    
    // also determine whether the label is unique
    moduleInfo.unique = moduleLabelCounts.at(moduleInfo.tag.label()) == 1;
    
  } // for
  return moduleMap;
  
} // details::DumpProvenanceTreeDetails::buildModuleMap()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(sbn::DumpProvenanceTree)


// -----------------------------------------------------------------------------

