#ifndef OCN_YAKL_DIAGSOLVE
#define OCN_YAKL_DIAGSOLVE

#include "mpas_ocn_yakl_types.hxx"

namespace diag_solve
{
extern int nCells,
    nEdges,
    nVertices,
    nVertLevels
    ;
    
extern
d_double_1d_t   * ssh
    ;

extern
d_double_2d_t   * normalVelocity,
                * layerThicknessEdge
                ;


};  // namespace diag_solve

#endif
