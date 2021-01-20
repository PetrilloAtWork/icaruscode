art_make( LIB_LIBRARIES
  larcorealg_Geometry
  lardataobj_RecoBase
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${MF_MESSAGELOGGER}
  ROOT::Core
  )

install_headers()
install_fhicl()
install_source()
