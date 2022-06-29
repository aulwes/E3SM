#include "mpas_ocn_mesh_c.hxx"

struct f2c_type
{
    int shape[10];
    void * c_loc;
};

namespace mesh
{

int nCellsAll,
    nEdgesAll,
    nVertices,
    nVertLevels
    ;
    
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

d_double_1d_t   * areaCell,
                * areaTriangle,
                * invAreaCell,
                * invAreaTriangle,
                * bottomDepth,
                * fVertex,
                * fEdge,
                * dcEdge,
                * dvEdge
                ;
                
d_double_2d_t   * edgeSignOnCell,
                * edgeSignOnVertex,
                * kiteAreasOnVertex,
                * advCoefs,
                * advCoefs3rd,
                * weightsOnEdge,
                * highOrderAdvectionMask
                ;
}

using namespace mesh;

extern "C" int * c_maxLevelCell = nullptr;
extern "C" int * c_minLevelCell = nullptr;
extern "C" int * c_minLevelEdgeTop = nullptr;
extern "C" int * c_maxLevelEdgeTop = nullptr;
extern "C" int * c_minLevelVertexTop = nullptr;
extern "C" int * c_maxLevelEdgeBot = nullptr;
extern "C" int * c_minLevelEdgeBot = nullptr;
extern "C" int * c_maxLevelVertexBot = nullptr;
extern "C" int * c_cellsOnEdge = nullptr;
extern "C" int * c_cellsOnCell = nullptr;
extern "C" int * c_edgesOnCell = nullptr;
extern "C" int * c_edgesOnVertex = nullptr;
extern "C" int * c_nEdgesOnEdge = nullptr;
extern "C" int * c_nEdgesOnCell = nullptr;
extern "C" double * c_fVertex = nullptr;
extern "C" double * c_dvEdge = nullptr;
extern "C" double * c_dcEdge = nullptr;
extern "C" double * c_areaCell = nullptr;
extern "C" double * c_areaTriangle = nullptr;
extern "C" double * c_invAreaCell = nullptr;
extern "C" double * c_invAreaTriangle = nullptr;
extern "C" double * c_edgeSignOnCell = nullptr;
extern "C" double * c_edgeSignOnVertex = nullptr;
extern "C" double * c_highOrderAdvectionMask = nullptr;
extern "C" double * c_bottomDepth = nullptr;

extern "C" ocn_yakl_type c_nAdvCellsForEdge;
extern "C" ocn_yakl_type c_advCellsForEdge;
extern "C" ocn_yakl_type c_advCoefs;
extern "C" ocn_yakl_type c_advCoefs3rd;
extern "C" ocn_yakl_type c_kiteIndexOnCell;
extern "C" ocn_yakl_type c_kiteAreasOnVertex;
extern "C" ocn_yakl_type c_verticesOnCell;
extern "C" ocn_yakl_type c_weightsOnEdge;
extern "C" ocn_yakl_type c_verticesOnEdge;
extern "C" ocn_yakl_type c_edgesOnEdge;
extern "C" ocn_yakl_type c_cellsOnVertex;
extern "C" ocn_yakl_type c_fEdge;

extern "C"
void ocn_tracer_advect_yakl_init()
{
    nAdvCellsForEdge = yakl_wrap_array("nAdvCellsForEdge", 
            static_cast<int *>(c_nAdvCellsForEdge.ptr), c_nAdvCellsForEdge.shape[0]);
    advCellsForEdge = yakl_wrap_array("advCellsForEdge", 
            static_cast<int *>(c_advCellsForEdge.ptr), c_advCellsForEdge.shape[0],
            c_advCellsForEdge.shape[1]);
    advCoefs = yakl_wrap_array("advCoefs", 
            static_cast<double *>(c_advCoefs.ptr), c_advCoefs.shape[0],
            c_advCoefs.shape[1]);
    advCoefs3rd = yakl_wrap_array("advCoefs3rd", 
            static_cast<double *>(c_advCoefs3rd.ptr), c_advCoefs3rd.shape[0],
            c_advCoefs3rd.shape[1]);
}



extern "C"
void ocn_mesh_yakl_init(int nCellsAll, int nEdgesAll, int nVertices, int nVertLevels, int maxNEdges, 
                        int mxConC, int advSize, int esonDim, int vertexDegree)
{
    //std::cerr << " nVertices, nVertLevels, nCells, nEdges = " << nVertices << " " << nVertLevels << " " << nCellsAll << " " << nEdgesAll << std::endl;
    //std::cerr << " esonDim = " << esonDim << std::endl;
    //std::cerr << " maxNEdges, mxConC, advSize = " << maxNEdges << " " << mxConC << " " << advSize << std::endl;
    mesh::nCellsAll = nCellsAll;
    mesh::nEdgesAll = nEdgesAll;
    mesh::nVertices = nVertices;
    mesh::nVertLevels = nVertLevels;
    
    std::cerr << "  c_maxLevelCell = " << c_maxLevelCell << std::endl;
    maxLevelCell = yakl_wrap_array("maxLevelCell", c_maxLevelCell, nCellsAll+1);
    minLevelCell = yakl_wrap_array("minLevelCell", c_minLevelCell, nCellsAll+1);
    bottomDepth = yakl_wrap_array("bottomDepth", c_bottomDepth, nCellsAll+1);
    minLevelEdgeTop = yakl_wrap_array("minLevelEdgeTop", c_minLevelEdgeTop, nEdgesAll+1);
    maxLevelEdgeTop = yakl_wrap_array("maxLevelEdgeTop", c_maxLevelEdgeTop, nEdgesAll+1);
    minLevelVertexTop = yakl_wrap_array("minLevelVertexTop", c_minLevelVertexTop, nVertices+1);
    minLevelEdgeBot = yakl_wrap_array("minLevelEdgeBot", c_minLevelEdgeBot, nEdgesAll+1);
    maxLevelEdgeBot = yakl_wrap_array("maxLevelEdgeBot", c_maxLevelEdgeBot, nEdgesAll+1);
    maxLevelVertexBot = yakl_wrap_array("maxLevelVertexBot", c_maxLevelVertexBot, nVertices+1);
    
    auto iptr = static_cast<int *>(c_kiteIndexOnCell.ptr);
    kiteIndexOnCell = yakl_wrap_array("kiteIndexOnCell", iptr,
                            c_kiteIndexOnCell.shape[0], c_kiteIndexOnCell.shape[1]);

    iptr = static_cast<int *>(c_verticesOnCell.ptr);
    verticesOnCell = yakl_wrap_array("verticesOnCell", iptr,
                            c_verticesOnCell.shape[0], c_verticesOnCell.shape[1]);

    auto dptr = static_cast<double *>(c_kiteAreasOnVertex.ptr);
    kiteAreasOnVertex = yakl_wrap_array("kiteAreasOnVertex", dptr,
                            c_kiteAreasOnVertex.shape[0], c_kiteAreasOnVertex.shape[1]);

    dptr = static_cast<double *>(c_weightsOnEdge.ptr);
    weightsOnEdge = yakl_wrap_array("weightsOnEdge", dptr,
                            c_weightsOnEdge.shape[0], c_weightsOnEdge.shape[1]);

    dptr = static_cast<double *>(c_fEdge.ptr);
    fEdge = yakl_wrap_array("fEdge", dptr, c_fEdge.shape[0]);

    iptr = static_cast<int *>(c_edgesOnEdge.ptr);
    edgesOnEdge = yakl_wrap_array("edgesOnEdge", iptr,
                            c_edgesOnEdge.shape[0], c_edgesOnEdge.shape[1]);

    iptr = static_cast<int *>(c_verticesOnEdge.ptr);
    verticesOnEdge = yakl_wrap_array("verticesOnEdge", iptr,
                            c_verticesOnEdge.shape[0], c_verticesOnEdge.shape[1]);

    iptr = static_cast<int *>(c_cellsOnVertex.ptr);
    cellsOnVertex = yakl_wrap_array("cellsOnVertex", iptr,
                            c_cellsOnVertex.shape[0], c_cellsOnVertex.shape[1]);

    cellsOnEdge = yakl_wrap_array("cellsOnEdge", c_cellsOnEdge, 2, nEdgesAll+1);
    cellsOnCell = yakl_wrap_array("cellsOnCell", c_cellsOnCell, mxConC, nCellsAll+1);
    fVertex = yakl_wrap_array("fVertex", c_fVertex, nVertices+1);
    dvEdge = yakl_wrap_array("dvEdge", c_dvEdge, nEdgesAll+1);
    dcEdge = yakl_wrap_array("dcEdge", c_dcEdge, nEdgesAll+1);
    areaCell = yakl_wrap_array("areaCell", c_areaCell, nCellsAll);
    areaTriangle = yakl_wrap_array("areaTriangle", c_areaTriangle, nVertices);
    edgesOnCell = yakl_wrap_array("edgesOnCell", c_edgesOnCell, maxNEdges, nCellsAll);
    edgesOnVertex = yakl_wrap_array("edgesOnVertex", c_edgesOnVertex, vertexDegree, nVertices+1);
    edgeSignOnCell = yakl_wrap_array("edgeSignOnCell", c_edgeSignOnCell, maxNEdges, nCellsAll+1);
    edgeSignOnVertex = yakl_wrap_array("edgeSignOnVertex", c_edgeSignOnVertex, esonDim, nVertices+1);
    highOrderAdvectionMask = yakl_wrap_array("highOrderAdvectionMask", c_highOrderAdvectionMask, nVertLevels, nEdgesAll+1);
    nEdgesOnCell = yakl_wrap_array("nEdgesOnCell", c_nEdgesOnCell, nCellsAll+1);
    nEdgesOnEdge = yakl_wrap_array("nEdgesOnEdge", c_nEdgesOnEdge, nEdgesAll+1);
    invAreaCell = yakl_wrap_array("invAreaCell", c_invAreaCell, nCellsAll);
    invAreaTriangle = yakl_wrap_array("invAreaTriangle", c_invAreaTriangle, nVertices);
}


extern "C"
void ocn_diag_mesh_yakl_update()
{
    yakl_update_device(mesh::kiteIndexOnCell, static_cast<int *>(c_kiteIndexOnCell.ptr));
    yakl_update_device(mesh::cellsOnVertex, static_cast<int *>(c_cellsOnVertex.ptr));
}

extern "C"
void ocn_mesh_yakl_update()
{
    yakl_update_device(mesh::edgeSignOnCell, c_edgeSignOnCell);
    yakl_update_device(mesh::edgeSignOnVertex, c_edgeSignOnVertex);
    yakl_update_device(mesh::highOrderAdvectionMask, c_highOrderAdvectionMask);
    yakl_update_device(mesh::maxLevelCell, c_maxLevelCell);
    yakl_update_device(mesh::nAdvCellsForEdge,  static_cast<int *>(c_nAdvCellsForEdge.ptr));
    yakl_update_device(mesh::advCellsForEdge,  static_cast<int *>(c_advCellsForEdge.ptr));
    yakl_update_device(mesh::advCoefs, static_cast<double *>(c_advCoefs.ptr));
    yakl_update_device(mesh::advCoefs3rd,  static_cast<double *>(c_advCoefs3rd.ptr));
    yakl_update_device(mesh::fVertex, c_fVertex);
    yakl_update_device(mesh::kiteAreasOnVertex, static_cast<double *>(c_kiteAreasOnVertex.ptr));
    //yakl_update_device(mesh::kiteIndexOnCell, static_cast<int *>(c_kiteIndexOnCell.ptr));
}
