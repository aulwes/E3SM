#ifndef OCN_YAKL_MESH
#define OCN_YAKL_MESH

#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"

namespace mesh
{
extern
d_int_1d_t      * maxLevelCell,
                * minLevelCell,
                * minLevelEdgeTop,
                * maxLevelEdgeTop,
                * minLevelVertexTop,
                * minLevelEdgeBot,
                * maxLevelEdgeBot,
                * maxLevelVertexBot,
                * nAdvCellsForEdge,
                * nEdgesOnEdge,
                * nEdgesOnCell
                ;

extern
d_int_2d_t      * advCellsForEdge,
                * cellsOnEdge,
                * cellsOnCell,
                * cellsOnVertex,
                * edgesOnVertex,
                * kiteIndexOnCell,
                * verticesOnCell,
                * verticesOnEdge,
                * edgesOnEdge,
                * edgesOnCell
                ;

extern
d_double_1d_t   * areaCell,
                * areaTriangle,
                * inverseAreaCell,
                * inverseAreaTriangle,
                * bottomDepth,
                * fVertex,
                * fEdge,
                * dcEdge,
                * dvEdge
                ;
                
extern
d_double_2d_t   * edgeSignOnCell,
                * edgeSignOnVertex,
                * kiteAreasOnVertex,
                * advCoefs,
                * advCoefs3rd,
                * weightsOnEdge,
                * highOrderAdvectionMask
                ;
} // namespace mesh

#endif
