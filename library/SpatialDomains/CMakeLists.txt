SET(SPATIAL_DOMAINS_SOURCES
Conditions.cpp
CrvSegGeom.cpp
CrvTetGeom.cpp
CrvTriGeom.cpp
GeomFactors.cpp
GeomFactorsPUMI.cpp
Geometry0D.cpp
Geometry1D.cpp
Geometry2D.cpp
Geometry3D.cpp
Geometry.cpp
HexGeom.cpp
InterfaceComponent.cpp
MeshComponents.cpp
MeshGraph.cpp
MeshGraph1D.cpp
MeshGraph2D.cpp
MeshGraph3D.cpp
PUMIMeshGraph2D.cpp
PUMIMeshGraph3D.cpp
pchSpatialDomains.cpp
PrismGeom.cpp
PyrGeom.cpp
QuadGeom.cpp
SegGeom.cpp
PointGeom.cpp
TetGeom.cpp
TriGeom.cpp
)

SET(SPATIAL_DOMAINS_HEADERS
Conditions.h
CrvSegGeom.h
CrvTetGeom.h
CrvTriGeom.h
Curve.hpp
GeomFactors.h
GeomFactorsPUMI.h
Geometry0D.h
Geometry1D.h
Geometry2D.h
Geometry3D.h
Geometry.h
HexGeom.h
InterfaceComponent.h
MeshComponents.h
MeshGraph1D.h
MeshGraph2D.h
MeshGraph3D.h
MeshGraph.h
PUMIMeshGraph2D.h
PUMIMeshGraph3D.h
pchSpatialDomains.h
PrismGeom.h
PyrGeom.h
QuadGeom.h
SegGeom.h
PointGeom.h
SpatialDomainsDeclspec.h
TetGeom.h
TriGeom.h
SpatialDomains.hpp
)
# SET(TINYXML_SOURCES
# /lore/zaided/nektar++-4.0.1/ThirdParty/tinyxml-2.6.2/tinyxml.cpp
# /lore/zaided/nektar++-4.0.1/ThirdParty/tinyxml-2.6.2/tinyxmlparser.cpp 
# /lore/zaided/nektar++-4.0.1/ThirdParty/tinyxml-2.6.2/tinyxmlerror.cpp 
# /lore/zaided/nektar++-4.0.1/ThirdParty/tinyxml-2.6.2/tinystr.cpp
# )

# SET(TINYXML_HEADERS
# /lore/zaided/nektar++-4.0.1/ThirdParty/tinyxml-2.6.2/tinystr.h
# /lore/zaided/nektar++-4.0.1/ThirdParty/tinyxml-2.6.2/tinyxml.h
# )


ADD_DEFINITIONS(-DSPATIAL_DOMAINS_EXPORTS -DMA_PARALLEL -DPUMI_PARALLEL -DCURVE)

# rpi libraries 
SET(PUMI_INSTALL_DIR
/lore/zaided/qlu/build
#/lore/qlu/scorec_software_dev
)

SET(PUMI_LIBS
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

INCLUDE_DIRECTORIES( ${PUMI_INSTALL_DIR}/include )

ADD_NEKTAR_LIBRARY(SpatialDomains lib ${NEKTAR_LIBRARY_TYPE} ${SPATIAL_DOMAINS_SOURCES} ${SPATIAL_DOMAINS_HEADERS})
TARGET_LINK_LIBRARIES(SpatialDomains
        optimized LibUtilities debug LibUtilities-g
        optimized StdRegions  debug StdRegions-g
        # optimized ${TINYXML_LIB} debug ${TINYXML_LIB}
        ${Boost_THREAD_LIBRARY}
		${Boost_DATE_TIME_LIBRARY}
		${Boost_IOSTREAMS_LIBRARY}
		${Boost_ZLIB_LIBRARY}
        ${PUMI_LIBS}
	)


INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/SpatialDomains COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# ADD_NEKTAR_LIBRARY(tinyxml lib STATIC ${TINYXML_SOURCES} ${TINYXML_HEADERS})

ADD_DEPENDENCIES(SpatialDomains tinyxml-2.6.2)
# SET_LAPACK_LINK_LIBRARIES(SpatialDomains)
