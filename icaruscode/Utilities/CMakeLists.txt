
art_make(
  EXCLUDE
    "SaveConfigurationIntoTFile_module.cc"
    "DuplicateEventTracker_service.cc"
  LIB_LIBRARIES
    ${ART_ROOT_IO_ROOTDB}
    ${SQLITE3}
    cetlib_except
    cetlib
    fhiclcpp
    hep_concurrency # required by fhicl::ParameterSetRegistry
    ROOT::Tree
    ROOT::Core
    ROOT::RIO
)

simple_plugin(SaveConfigurationIntoTFile "module"
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_CORE}
              canvas
              ${MF_MESSAGELOGGER}
              hep_concurrency
              ${FHICLCPP}
              ${CETLIB}
              ${ROOT_CORE}
)


simple_plugin(DuplicateEventTracker "service"
              icaruscode_Utilities
              ${ART_FRAMEWORK_PRINCIPAL}
              ${MF_MESSAGELOGGER}
)


install_headers()
install_source()
install_fhicl()
