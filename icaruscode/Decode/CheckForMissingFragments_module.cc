/**
 * @file   CheckForMissingFragments_module.cc
 * @brief  Checks the presence of raw data fragments from all components.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 26, 2021
 */

// artDAQ
#include "artdaq-core/Data/ContainerFragment.hh"
#include "artdaq-core/Data/Fragment.hh"

// framework libraries
#include "art/Framework/Core/SharedAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h" // art::ModuleLabelSelector
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <memory> // std::make_unique()
#include <unordered_map>
#include <vector>
#include <string>
#include <mutex>
#include <atomic> // std::atomic_uint
#include <functional> // std::reference_wrapper
#include <utility> // std::move()
#include <cassert>


// -----------------------------------------------------------------------------
namespace icarus { class CheckForMissingFragments; }
/**
 * @brief Checks the presence of data fragments from all configured components.
 * 
 * This module only emits warnings on the console when detecting missing
 * data fragments.
 * 
 * 
 * Input
 * ------
 * 
 * All `artdaq::Fragment` collections are considered for input.
 * 
 * 
 * Output
 * -------
 * 
 * Messages are printed on console. No other output is produced.
 * 
 * 
 * Configuration parameters
 * -------------------------
 * 
 * The following configuration parameters are supported:
 * 
 * * **DAQlabel** (input tag, default: `daq`): producer of the data fragments.
 * * **ComponentList** (see below): the list of requested components.
 * * **ThrowAfter** (integer, optional): if specified, throws an exception if
 *     this many events have been found with missing data fragments.
 *     By default, only an error is printed.
 * * **LogCategory** (string, default: `CheckForMissingFragments`):
 *     category tag used for messages to message facility.
 * 
 * The `ComponentList` is currently expressed as a FHiCL table, where the value
 * of each key a sequence of all fragment ID which belong to the component named
 * after the key. For example:
 * ```
 * ComponentList: {
 *   icaruspmt02wwbot: [ 0x2000, 0x2001, 0x2002 ]
 *   icaruspmt02wwtop: [ 0x2003, 0x2004, 0x2005 ]
 * }
 * ```
 * defines two components (`icaruspmt02wwbot` and `icaruspmt02wwtop`), each
 * associated to three fragment IDs.
 * 
 * 
 * Technical notes
 * ----------------
 * 
 * This module may take long to execute because it triggers the reading of
 * _all raw data_ in the event. It is quite a waste, since it's using only
 * the fragment IDs.
 * 
 * Container fragments are currently treated in a special way, in that the
 * ID of their first contained fragment is used when available, instead of the
 * ID of the container fragment itself.
 */
class icarus::CheckForMissingFragments: public art::SharedAnalyzer {
  
    public:
  
  /// Configuration of the module.
  struct Config {
    
    fhicl::DelegatedParameter ComponentList {
      fhicl::Name("ComponentList"),
      fhicl::Comment("list of required components")
      // mandatory
      };
    
    fhicl::Atom<std::string> DAQlabel {
      fhicl::Name("DAQlabel"),
      fhicl::Comment("data fragments producer module label"),
      "daq" // default
      };
    
    fhicl::OptionalAtom<unsigned int> ThrowAfter {
      fhicl::Name("ThrowAfter"),
      fhicl::Comment
        ("raise an exception after this many errors with missing fragments")
      };
    
    fhicl::Atom<std::string> LogCategory {
      fhicl::Name("LogCategory"),
      fhicl::Comment("category tag used for messages to message facility"),
      "PMTconfigurationExtraction" // default
      };
    
  }; // struct Config
  
  using Parameters = art::SharedAnalyzer::Table<Config>;
  
  
  /// Constructor: just reads the configuration.
  CheckForMissingFragments
    (Parameters const& config, art::ProcessingFrame const&);
  
  /// Perform the checks.
  virtual void analyze
    (art::Event const& event, art::ProcessingFrame const&) override;
  
  /// Prints a end-of-job summary.
  virtual void endJob(art::ProcessingFrame const&) override;
  
  
    private:
  
  // --- BEGIN -- Statistics for the components --------------------------------
  class ComponentStatData {
    
      public:
    struct Stat_t {
      std::string name; ///< Component name.
      unsigned int missing { 0U }; ///< Count of missing fragments.
      
      Stat_t& operator++ () { ++missing; return *this; }
    }; // Stat_t
    
    using ComponentStat_t = std::pair<std::string, Stat_t>;
    using ComponentStats_t = std::vector<ComponentStat_t>;
    
    using iterator = ComponentStats_t::iterator;
    using const_iterator = ComponentStats_t::const_iterator;
    using value_type = ComponentStats_t::value_type;
    
    // --- BEGIN -- Query and iteration interface ------------------------------
    bool empty() const { return fStats.empty(); }
    std::size_t size() const { return fStats.size(); }
    const_iterator cbegin() const { return fStats.begin(); }
    const_iterator cend() const { return fStats.end(); }
    const_iterator begin() const { return cbegin(); }
    const_iterator end() const { return cend(); }
    iterator begin() { return fStats.begin(); }
    iterator end() { return fStats.end(); }
    
    Stat_t const* cfind(std::string const& name) const;
    Stat_t const* find(std::string const& name) const { return cfind(name); }
    Stat_t* find(std::string const& name);
    // --- END ---- Query and iteration interface ------------------------------
    
    
    // --- BEGIN -- Initialization interface -----------------------------------
    /// Adds a component with the specified name, if not there already.
    /// @return whether the component was actually added
    bool add(std::string const& name);
    
    // --- END ---- Initialization interface -----------------------------------
    
      public:
    ComponentStats_t fStats;
  }; // ComponentStatData
  
  
  /// Component statistics. Owns data, and exposes it exclusively.
  class ComponentStats {
    
    mutable std::mutex fMutex; ///< Access control.
    ComponentStatData fStats; ///< Statistics.
    
      public:
    
    using LockGuard_t = std::unique_lock<std::mutex>;
    using StatsRef_t = std::reference_wrapper<ComponentStatData>;
    using StatsConstPtr_t = std::reference_wrapper<ComponentStatData const>;
    
    ComponentStats() = default;
    ComponentStats(ComponentStatData data): fStats(std::move(data)) {}
    
    /// Returns a owning lock and a reference to this data.
    [[nodiscard]]
    std::pair<LockGuard_t, ComponentStatData*> exclusiveDataPtr()
      { return { LockGuard_t{ fMutex }, &fStats }; }
    
    /// Returns a constant reference to this data. Not thread-safe.
    ComponentStatData const& data() const
      { return fStats; }
    
    
  }; // class ComponentStats

  /// Type of counters of fragments by component.
  using EventComponents_t = std::unordered_map<std::string, unsigned int>;
  
  // --- END ---- Statistics for the components --------------------------------
  
  
  // --- BEGIN -- Association of components with fragments ---------------------
  
  /// Associates fragment IDs to components.
  class FragmentDB {
    
      public:
    
    using fragment_id_t = artdaq::Fragment::fragment_id_t;
    
    /// Registers the specified IDs for the component `name`.
    void registerComponent
      (std::string const& name, std::vector<fragment_id_t> const& IDs);
    
    /// Returns a ordered list of component names.
    std::vector<std::string> const& componentList() const
      { return fComponents; }
    
    /// Returns the name of the component registered for `id`, empty if none.
    std::string componentOf(fragment_id_t id) const
      {
        if (auto it = fFragments.find(id); it != fFragments.end()) {
          ComponentIndex_t const index = it->second;
          assert(index < fComponents.size());
          return fComponents[index];
        }
        else return {};
      }
    
      private:
    
    using ComponentIndex_t = std::size_t; ///< Internal index for a component.
    
    /// Invalid index mnemonic.
    static constexpr ComponentIndex_t NoIndex
      = std::numeric_limits<ComponentIndex_t>::max();
    
    std::vector<std::string> fComponents;
    std::unordered_map<fragment_id_t, ComponentIndex_t> fFragments;
    
    /// Returns the internal index of the specified component, or `NoIndex`.
    ComponentIndex_t indexOf(std::string const& name) const;
    
    /// Returns the name of component with the specified `index`, empty if none.
    std::string nameOf(ComponentIndex_t index) const;
    
  }; // FragmentDB
  
  // --- END ---- Association of components with fragments ---------------------

  // --- BEGIN -- Configuration parameters -------------------------------------
  
  /// Fragment ID to component name map.
  FragmentDB fFragmentMap;
  
  std::string const fDAQlabel; ///< Label of producer of DAQ fragments.
  
  /// Throw an exception after this many errors.
  unsigned int const fThrowAfterErrors = 0U;
  
  std::string fLogCategory; ///< Category tag for messages.
  
  // --- END ---- Configuration parameters -------------------------------------
  
  std::atomic_uint fNEvents { 0U }; ///< Count of seen events.
  std::atomic_uint fNErrors { 0U }; ///< Count of events with missing fragments.
  
  ComponentStats fComponentStats; ///< Expected components and their usage.
  
  /// Returns the sequence of names of configured components.
  std::vector<std::string> const& componentList() const
    { return fFragmentMap.componentList(); }
  
  /// Returns all the fragment handles from the specified module label.
  std::vector<art::Handle<artdaq::Fragments>> readFragments
    (art::Event const& event) const;
  
  /// Returns the name of the component which produced `fragment`.
  std::string extractComponent(artdaq::Fragment const& fragment) const;
  
  /// Counts the missing `components` out of the ones in the list.
  /// @return the number of missing components
  unsigned int registerEventComponents
    (art::EventID const& id, EventComponents_t const& components);
  
  
  /// Prints the internal counters into the specified output stream.
  template <typename Stream>
  void printCounters(Stream&& out) const;
  
  
  /// Returns a new `FragmentDB` object from the specified configuration.
  static FragmentDB initializeFragmentMap(fhicl::ParameterSet const& config);
  
  /// Returns a `ComponentStatData` object filled with the specified components.
  static ComponentStatData initializeComponentStats
    (std::vector<std::string> const& components);
  
  /// Returns whether `fragment` is a container fragment.
  static bool isContainer(artdaq::Fragment const& fragment);
  
  /// Returns the first element of `fragment` if it is a container,
  /// `fragment` itself otherwise. Null if `fragment` is an empty container.
  static artdaq::FragmentPtr extractFirstFragment(artdaq::Fragment fragment);
  
  
}; // icarus::CheckForMissingFragments



// -----------------------------------------------------------------------------
// ---  implementation
// -----------------------------------------------------------------------------
// ---  icarus::CheckForMissingFragments::ComponentStatData
// -----------------------------------------------------------------------------
auto icarus::CheckForMissingFragments::ComponentStatData::cfind
  (std::string const& name) const -> Stat_t const*
{
  // so far I stick to a linear search
  for (auto& stat: fStats) if (stat.first == name) return &(stat.second);
  return nullptr;
} // icarus::CheckForMissingFragments::ComponentStatData::cfind()


// -----------------------------------------------------------------------------
auto icarus::CheckForMissingFragments::ComponentStatData::find
  (std::string const& name) -> Stat_t*
{
  // so far I stick to a linear search
  for (auto& stat: fStats) if (stat.first == name) return &(stat.second);
  return nullptr;
} // icarus::CheckForMissingFragments::ComponentStatData::find()


// -----------------------------------------------------------------------------
bool icarus::CheckForMissingFragments::ComponentStatData::add
  (std::string const& name)
{
  Stat_t* stat = find(name);
  if (stat) return false;
  fStats.emplace_back(name, Stat_t{ name });
  return true;
} // icarus::CheckForMissingFragments::ComponentStatData::add()


// -----------------------------------------------------------------------------
// ---  icarus::CheckForMissingFragments::FragmentDB
// -----------------------------------------------------------------------------
void icarus::CheckForMissingFragments::FragmentDB::registerComponent
  (std::string const& name, std::vector<fragment_id_t> const& IDs)
{
  ComponentIndex_t index = indexOf(name);
  if (index == NoIndex) { // otherwise we are extending an existing component
    index = fComponents.size();
    fComponents.push_back(name);
  }
  for (fragment_id_t ID: IDs) {
    auto it = fFragments.find(ID);
    if (it == fFragments.end()) {
      // register a new ID
      fFragments[ID] = index;
    }
    else if (it->first != index) { // otherwise it's again the same registration
      throw cet::exception("CheckForMissingFragments")
        << "Attempted to register fragment ID " << std::hex << ID << std::dec
        << " with both '" << nameOf(it->second) << "' and '" << name
        << "' components!\n";
    }
  } // for
    
} // icarus::CheckForMissingFragments::FragmentDB::registerComponent()


// -----------------------------------------------------------------------------
auto icarus::CheckForMissingFragments::FragmentDB::indexOf
  (std::string const& name) const -> ComponentIndex_t
{
  // linear search for now
  for (ComponentIndex_t index = 0; index < fComponents.size(); ++index)
    if (fComponents[index] == name) return index;
  return NoIndex;
} // icarus::CheckForMissingFragments::FragmentDB::indexOf()


// -----------------------------------------------------------------------------
std::string icarus::CheckForMissingFragments::FragmentDB::nameOf
  (ComponentIndex_t index) const
{
  return (index < fComponents.size())? fComponents[index]: std::string{};
} // icarus::CheckForMissingFragments::FragmentDB::nameOf()


// -----------------------------------------------------------------------------
// ---  icarus::CheckForMissingFragments
// -----------------------------------------------------------------------------
icarus::CheckForMissingFragments::CheckForMissingFragments
  (Parameters const& config, art::ProcessingFrame const&)
  : art::SharedAnalyzer{ config }
  // configuration
  , fFragmentMap
    { initializeFragmentMap(config().ComponentList.get<fhicl::ParameterSet>()) }
  , fDAQlabel{ config().DAQlabel() }
  , fThrowAfterErrors{ config().ThrowAfter().value_or(0U) }
  , fLogCategory{ config().LogCategory() }
  // derived quantities
  , fComponentStats{ initializeComponentStats(componentList()) }
{
  
  async<art::InEvent>();
  
  // we don't go `consumesMany()` because we use a selector,
  // so we don't actually consume *all* the fragment collections
  mayConsumeMany<artdaq::Fragments>();
  
  // --- BEGIN -- Configuration dump on screen ---------------------------------
  ComponentStatData const& componentStats = fComponentStats.data();
  mf::LogInfo log { fLogCategory };
  log << "Configuration:"
    "\n * data fragments from tag: '" << fDAQlabel << "'"
    "\n * required components (" << componentStats.size() << "):"
    ;
  for (auto const& nameAndStats: componentStats)
    log << " " << nameAndStats.second.name;
  if (fThrowAfterErrors > 0U) {
    log << "\n * throw an exception after " << fThrowAfterErrors
      << " events with missing events";
  }
  // --- END ---- Configuration dump on screen ---------------------------------
  
} // icarus::CheckForMissingFragments::CheckForMissingFragments()


// -----------------------------------------------------------------------------
void icarus::CheckForMissingFragments::analyze
  (art::Event const& event, art::ProcessingFrame const&)
{
  
  ++fNEvents;
  
  //
  // fill the list of components from all fragment data products:
  //
  EventComponents_t components;
  for (art::Handle<artdaq::Fragments> const& fragsHandle: readFragments(event))
  {
    
    // this will be printed at the end of the iteration:
    mf::LogTrace log{ fLogCategory };
    log << "Fragments '" << fragsHandle.provenance()->inputTag().encode()
      << "' (" << fragsHandle->size() << "):"
      ;
    
    //
    // process each fragment in the handle
    //
    for (artdaq::Fragment const& fragment: *fragsHandle) {
      
      //
      // extract the first actual data fragment if the fragment is a container
      //
      artdaq::FragmentPtr const firstFragment = extractFirstFragment(fragment);
      
      //
      // extract the component of the data fragment and record it
      //
      
      // out of curiosity, we attempt extracting the component also out of
      // empty container fragments:
      artdaq::Fragment const& keyFragment
        = firstFragment? *firstFragment: fragment;
      std::string const componentName = extractComponent(keyFragment);
      if (componentName.empty()) {
        mf::LogTrace log{ fLogCategory };
        log << "Fragment ID=" << std::hex << keyFragment.fragmentID() << std::dec;
        if (isContainer(fragment)) 
          log << " (container ID=" << std::hex << fragment.fragmentID() << std::dec << ")";
        log << " is not registered from any component";
      } // if no known component
      
      log << " " << componentName;
      ++(components[componentName]);
      
    } // for fragments
    
  } // for fragment data products
  
  //
  // check if there are missing components and fill the relative counters
  //
  registerEventComponents(event.id(), components);
  
  if ((fThrowAfterErrors > 0U) && (fNErrors >= fThrowAfterErrors)) {
    
    printCounters(mf::LogError{ fLogCategory });
    
    throw cet::exception("CheckForMissingFragments")
      << "Found at least " << fNErrors << " events with missing fragments!\n";
    
  }
  
} // icarus::CheckForMissingFragments::analyze()


// -----------------------------------------------------------------------------
void icarus::CheckForMissingFragments::endJob(art::ProcessingFrame const&) {
  
  printCounters(mf::LogInfo{ fLogCategory });
  
} // icarus::CheckForMissingFragments::endJob()


// -----------------------------------------------------------------------------
std::vector<art::Handle<artdaq::Fragments>>
icarus::CheckForMissingFragments::readFragments(art::Event const& event) const {
  
  return
    event.getMany<artdaq::Fragments>(art::ModuleLabelSelector{ fDAQlabel });
  
} // icarus::CheckForMissingFragments::readFragments()


// -----------------------------------------------------------------------------
std::string icarus::CheckForMissingFragments::extractComponent
  (artdaq::Fragment const& fragment) const
{
  return fFragmentMap.componentOf(fragment.fragmentID());
} // icarus::CheckForMissingFragments::extractComponent()


// -----------------------------------------------------------------------------
unsigned int icarus::CheckForMissingFragments::registerEventComponents
  (art::EventID const& id, EventComponents_t const& components)
{
  std::vector<std::string> missingComponents;
  
  // locks the statistics
  auto const& [ guard, componentStatsPtr ] = fComponentStats.exclusiveDataPtr();
  for (auto& [ name, stats ]: *componentStatsPtr) {
    
    if (components.count(name) > 0U) continue;
    
    ++stats.missing;
    missingComponents.push_back(name);
    
  } // for registered components
  
  // should we report unregistered components?
  
  if (missingComponents.empty()) return 0;
  
  mf::LogTrace log { fLogCategory };
  log << id << " missing " << missingComponents.size() << " components:";
  for (std::string const& name: missingComponents) log << " " << name;
  ++fNErrors;
  
  return missingComponents.size();
  
} // icarus::CheckForMissingFragments::registerEventComponents()


// -----------------------------------------------------------------------------
template <typename Stream>
void icarus::CheckForMissingFragments::printCounters(Stream&& out) const {
  
  if (fNErrors > 0U) {
    out << "Detected " << fNErrors << "/" << fNEvents
      << " events with missing fragment."
      << "\nEvents with missing fragments from each component:";
    
    unsigned int nMissingComponents = 0U;
    for(auto const& [ key, stats ]: fComponentStats.data()) {
      if (stats.missing == 0U) continue;
      out << "\n - '" << stats.name << ": missed in " << stats.missing
        << " events";
      if (stats.name != key) out << " [component key: '" << key << "']";
      ++nMissingComponents;
    } // for
    out << "\nOverall, " << nMissingComponents << " components were missing.";
  }
  else {
    out << "All " << fNEvents << " events had no missing fragment.";
  }
  
} // icarus::CheckForMissingFragments::printCounters()


// -----------------------------------------------------------------------------
auto icarus::CheckForMissingFragments::initializeFragmentMap
  (fhicl::ParameterSet const& config) -> FragmentDB
{
  FragmentDB fragmentDB;
  
  for (std::string const& name: config.get_names()) {
    
    try {
      fragmentDB.registerComponent
        (name, config.get<std::vector<artdaq::Fragment::fragment_id_t>>(name));
    }
    catch(cet::exception const& e) {
      throw art::Exception{ art::errors::Configuration, "", e }
        << "Component '" << name
        << "' is assigned an invalid list of fragment IDs!\n";
    }
    
  } // for all configured components
  
  return fragmentDB;
} // icarus::CheckForMissingFragments::initializeFragmentMap()


// -----------------------------------------------------------------------------
bool icarus::CheckForMissingFragments::isContainer
  (artdaq::Fragment const& fragment)
  { return fragment.type() == artdaq::Fragment::ContainerFragmentType; }


// -----------------------------------------------------------------------------
auto icarus::CheckForMissingFragments::initializeComponentStats
  (std::vector<std::string> const& components) -> ComponentStatData
{
  ComponentStatData data;
  
  for (std::string const& name: components) data.add(name);
  
  return data;
} // icarus::CheckForMissingFragments::initializeComponentStats()


// -----------------------------------------------------------------------------
artdaq::FragmentPtr icarus::CheckForMissingFragments::extractFirstFragment
  (artdaq::Fragment fragment)
{

  if (!isContainer(fragment))
    return std::make_unique<artdaq::Fragment>(std::move(fragment));
  
  artdaq::ContainerFragment const containerFragment{ std::move(fragment) };
  
  return (containerFragment.block_count() == 0)
    ? artdaq::FragmentPtr{}: containerFragment.at(0);
  
} // icarus::CheckForMissingFragments::extractFirstFragment()


// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::CheckForMissingFragments)


// -----------------------------------------------------------------------------
