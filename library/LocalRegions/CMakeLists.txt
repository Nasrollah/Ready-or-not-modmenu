SET(LOCAL_REGIONS_SOURCES
Expansion.cpp
Expansion0D.cpp
Expansion1D.cpp
Expansion2D.cpp
Expansion3D.cpp
QuadExp.cpp
HexExp.cpp
MatrixKey.cpp
NodalTetExp.cpp
NodalTriExp.cpp
PointExp.cpp
PrismExp.cpp
PyrExp.cpp
SegExp.cpp
TetExp.cpp
TriExp.cpp
)

SET(LOCAL_REGIONS_HEADERS
Expansion.h
Expansion0D.h
Expansion1D.h
Expansion2D.h
Expansion3D.h
HexExp.h
LinSys.hpp
LocalRegions.hpp
LocalRegionsDeclspec.h
MatrixKey.h
NodalTetExp.h
NodalTriExp.h
PointExp.h
PrismExp.h
PyrExp.h
QuadExp.h
SegExp.h
TetExp.h
TriExp.h
)

# rpi libraries 
SET(PUMI_INSTALL_DIR
/lore/zaided/qlu/build
#/lore/qlu/scorec_software_dev
)

SET(PUMI_LIBS
${PUMI_INSTALL_DIR}/lib/libcurved_meshadapt_nektarpp.a
${PUMI_INSTALL_DIR}/lib/libadv_loop_base.a
${PUMI_INSTALL_DIR}/lib/libmeshadapt.a
${PUMI_INSTALL_DIR}/lib/libpumi_geom_meshmodel.a
${PUMI_INSTALL_DIR}/lib/libapf_pumi.a
-L${PUMI_INSTALL_DIR}/lib/ -lpumi_mesh
/lore/qlu/packages/lib/libzoltan.a
/lore/qlu/packages/lib/libparmetis.a
/lore/qlu/packages/lib/libmetis.a
/lore/qlu/packages/nektar++-3.4.0/ThirdParty/dist/lib/libmetis.a
${PUMI_INSTALL_DIR}/lib/libpumi_geom_parasolid.a
${PUMI_INSTALL_DIR}/lib/libpumi_geom.a
${PUMI_INSTALL_DIR}/lib/libpumi_util.a
${PUMI_INSTALL_DIR}/lib/libapf.a
${PUMI_INSTALL_DIR}/lib/libpcu.a
/usr/local/parasolid/17.0/shared_object/libpskernel.so
)

INCLUDE_DIRECTORIES( ${PUMI_INSTALL_DIR}/include /lore/qlu/scorec_software_dev/curved_meshadapt_omega3p)


ADD_DEFINITIONS(-DLOCAL_REGIONS_EXPORTS)

ADD_NEKTAR_LIBRARY(LocalRegions lib ${NEKTAR_LIBRARY_TYPE} ${LOCAL_REGIONS_SOURCES} ${LOCAL_REGIONS_HEADERS})
TARGET_LINK_LIBRARIES(LocalRegions
        optimized LibUtilities   debug LibUtilities-g
        optimized SpatialDomains debug SpatialDomains-g
        optimized StdRegions     debug StdRegions-g
        # optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
        ${Boost_THREAD_LIBRARY}
        ${Boost_IOSTREAMS_LIBRARY}
        ${Boost_DATE_TIME_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${ZLIB_LIBRARY}
        ${PUMI_LIBS}
	)
#SET_LAPACK_LINK_LIBRARIES(LocalRegions)

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/LocalRegions COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
