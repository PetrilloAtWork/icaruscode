add_subdirectory(tools)
add_subdirectory(Algorithms)

art_make( EXCLUDE AnodeCathodeCrossingClusters_module.cc
          MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nug4_ParticleNavigation
                           lardataobj_Simulation
                           lardata_Utilities
                           larevt_Filters
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardataobj_RawData
                           larevt_CalibrationDBI_Providers
                           nurandom_RandomUtils_NuRandomService_service
			   sbnobj_Common_Analysis
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           cetlib_except
                           ${CLHEP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           icaruscode_RecoUtils
        )

simple_plugin(AnodeCathodeCrossingClusters module
  icaruscode_Analysis_Algorithms
  icaruscode_TPC_SignalProcessing_HitFinder
  larreco_ClusterFinder
  larreco_RecoAlg_ClusterRecoUtil
  lardata_Utilities
  larcore_Geometry_Geometry_service
  lardataalg_DetectorInfo
  larcorealg_Geometry
  lardataobj_RecoBase
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${MF_MESSAGELOGGER}
  ROOT::Core
  )

#install_headers()
install_fhicl()
install_source()
