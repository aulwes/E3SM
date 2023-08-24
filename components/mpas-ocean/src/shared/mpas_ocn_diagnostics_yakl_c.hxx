#ifndef OCN_YAKL_DIAGSOLVE
#define OCN_YAKL_DIAGSOLVE

#include "mpas_ocn_yakl_types.hxx"

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
                * cTemp
                ;

extern
d_double_2d_t   * normalVelocity,
                * kineticEnergyCell,
                * layerThickEdgeMean,
                * layerThicknessEdge
                ;


};  // namespace diag_solve

#endif
