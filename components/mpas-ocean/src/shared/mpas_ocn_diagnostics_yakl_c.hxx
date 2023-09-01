#ifndef OCN_YAKL_DIAGSOLVE
#define OCN_YAKL_DIAGSOLVE

#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"

namespace diag_solve
{
extern int
    nCells,
    nEdges,
    nVertices,
    nVertLevels
    ;

extern yakl::Stream
     stream1,
     stream2;

extern yakl::Event
    event1;

extern
d_double_1d_t   * ssh
    ;

// temp scratch space
extern
d_double_2d_t   * bTemp,
                * rTemp,
                * cTemp,
                * pressure,
                * density,
                * zMid,
                * vertAleTransportTop,
                * normalizedRelativeVorticityEdge,
                * normalizedPlanetaryVorticityEdge
                ;

extern
d_double_2d_t   * normalVelocity,
                * kineticEnergyCell,
                * thermExpCoeff,
                * salineContractCoeff,
                * layerThickEdgeMean,
                * layerThickEdgeFlux,
                * layerThicknessEdge
                ;

extern
d_double_3d_t   * activeTracers
                ;
extern "C" ocn_yakl_type c_thermExpCoeff;
extern "C" ocn_yakl_type c_salineContractCoeff;

extern "C" ocn_yakl_type c_activeTracers;
extern "C" ocn_yakl_type c_layerThickEdgeFlux;
extern "C" ocn_yakl_type c_kineticEnergyCell;
extern "C" ocn_yakl_type c_normalVelocity;
extern "C" ocn_yakl_type c_vertAleTransportTop;
extern "C" ocn_yakl_type c_pressure;
extern "C" ocn_yakl_type c_diag_density;
extern "C" ocn_yakl_type c_zMid;

};  // namespace diag_solve

#endif
