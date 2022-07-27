#include "mpas_ocn_diagnostics_yakl_c.hxx"
#include "mpas_ocn_yakl_c.hxx"
#include "mpas_ocn_diagnostics_variables.hxx"
#include "mpas_ocn_mesh_c.hxx"
#include "mpas_ocn_yakl_types.hxx"

#define GET_DPTR(v) (static_cast<double *>((diag_solve::c_##v).ptr))
#define DEBUG 0

namespace diag_solve
{
int nCells,
    nEdges,
    nVertices,
    nVertLevels
    ;

yakl::yakl_stream_t     stream1 = 0;
yakl::yakl_stream_t     stream2 = 0;
yakl::yakl_event_t     event1 = 0;
yakl::yakl_event_t     event2 = 0;

d_double_1d_t   * ssh,
                * seaIcePressure,
                * surfacePressure,
                * atmosphericPressure
                ; 

d_double_2d_t   * div_hu,
                * div_huTransport,
                * div_huGMBolus;

d_double_2d_t   * normalVelocity,
                * normalTransportVelocity,
                * vertVelocityTop,
                * vertTransportVelocityTop,
                * vertGMBolusVelocityTop,
                * divergence,
                * kineticEnergyCell,
                * normalGMBolusVelocity,
                * BruntVaisalaFreqTop,
                * displacedDensity,
                * density,
                * zMid,
                * zTop,
                * RiTopOfCell,
                * circulation,
                * relativeVorticity,
                * relativeVorticityCell,
                * layerThickness,
                * layerThicknessEdge,
                * tangentialVelocity,
                * edgeAreaFractionOfCell,
                * montgomeryPotential,
                * pressure,
                * normalizedRelativeVorticityEdge,
                * normalizedPlanetaryVorticityEdge,
                * normalizedRelativeVorticityCell,
                * layerThickEdgeFlux,
                * layerThickEdgeMean
                ;

d_double_2d_t   * normalizedPlanetaryVorticityVertex,
                * normalizedRelativeVorticityVertex,
                * vorticityGradientNormalComponent,
                * vorticityGradientTangentialComponent
                ;

extern "C" ocn_yakl_type c_kineticEnergyCell;
extern "C" ocn_yakl_type c_relativeVorticityCell;
extern "C" ocn_yakl_type c_divergence;
extern "C" ocn_yakl_type c_normalGMBolusVelocity;
extern "C" ocn_yakl_type c_normalTransportVelocity;
extern "C" ocn_yakl_type c_vertTransportVelocityTop;
extern "C" ocn_yakl_type c_vertGMBolusVelocityTop;
extern "C" ocn_yakl_type c_vertVelocityTop;
extern "C" ocn_yakl_type c_tangentialVelocity;
extern "C" ocn_yakl_type c_normRelVortEdge;
extern "C" ocn_yakl_type c_normPlanetVortEdge;
extern "C" ocn_yakl_type c_normalizedRelativeVorticityCell;
extern "C" ocn_yakl_type c_displacedDensity;
extern "C" ocn_yakl_type c_density;
extern "C" ocn_yakl_type c_zMid;
extern "C" ocn_yakl_type c_zTop;
extern "C" ocn_yakl_type c_ssh;
extern "C" ocn_yakl_type c_bottomDepth;
extern "C" ocn_yakl_type c_RiTopOfCell;
extern "C" ocn_yakl_type c_BruntVaisalaFreqTop;
extern "C" ocn_yakl_type c_edgeAreaFractionOfCell;
extern "C" ocn_yakl_type c_pressure;
extern "C" ocn_yakl_type c_surfacePressure;
extern "C" ocn_yakl_type c_montgomeryPotential;

};

extern "C"
void ocn_diagnostics_yakl_init(int nCells, int nEdges, int nVertices, int nVertLevels)
{
    diag_solve::nCells = nCells;
    diag_solve::nEdges = nEdges;
    diag_solve::nVertices = nVertices;
    diag_solve::nVertLevels = nVertLevels;
    
    using namespace diag_solve;
    
    yakl::streamCreate(&stream1);
    yakl::streamCreate(&stream2);
    yakl::eventCreate(&event1);
    yakl::eventCreate(&event2);

    if ( DEBUG) {
        std::cerr << " ocn_diagnostics_yakl_init: c_relativeVorticityCell shape = "
            << c_relativeVorticityCell.shape[0] << " " << c_relativeVorticityCell.shape[1]  << std::endl;
        std::cerr << " ocn_diagnostics_yakl_init: nVertLevels, nVertices = "
            << nVertLevels << " " << nVertices << std::endl;
    }
    normalVelocity = yakl_create_real("normalVelocity", nVertLevels, nEdges);
    layerThickness = yakl_create_real("layerThickness", nVertLevels, nCells);
    layerThickEdgeMean = yakl_create_real("layerThickEdgeMean", nVertLevels, nEdges);
    layerThickEdgeFlux = yakl_create_real("layerThickEdgeFlux", nVertLevels, nEdges);
    layerThicknessEdge = yakl_create_real("layerThicknessEdge", nVertLevels, nEdges);
    circulation = yakl_create_real("circulation", nVertLevels, nVertices);
    relativeVorticity = yakl_create_real("relativeVorticity", nVertLevels, nVertices);
    relativeVorticityCell = yakl_create_real("relativeVorticityCell", 
                    c_relativeVorticityCell.shape[0], c_relativeVorticityCell.shape[1]);
    divergence = yakl_create_real("divergence", 
                    c_divergence.shape[0], c_divergence.shape[1]);
    kineticEnergyCell = yakl_create_real("kineticEnergyCell", 
                    c_kineticEnergyCell.shape[0], c_kineticEnergyCell.shape[1]);
    normalGMBolusVelocity = yakl_create_real("normalGMBolusVelocity", 
                    c_normalGMBolusVelocity.shape[0], c_normalGMBolusVelocity.shape[1]);
    normalTransportVelocity = yakl_create_real("normalTransportVelocity", 
                    c_normalTransportVelocity.shape[0], c_normalTransportVelocity.shape[1]);
    vertVelocityTop = yakl_create_real("vertVelocityTop", 
                    c_vertVelocityTop.shape[0], c_vertVelocityTop.shape[1]);
    vertTransportVelocityTop = yakl_create_real("vertTransportVelocityTop", 
                    c_vertTransportVelocityTop.shape[0], c_vertTransportVelocityTop.shape[1]);
    vertGMBolusVelocityTop = yakl_create_real("vertGMBolusVelocityTop", 
                    c_vertGMBolusVelocityTop.shape[0], c_vertGMBolusVelocityTop.shape[1]);
    tangentialVelocity = yakl_create_real("tangentialVelocity", 
                    c_tangentialVelocity.shape[0], c_tangentialVelocity.shape[1]);

    displacedDensity = yakl_create_real("displacedDensity", 
                    c_displacedDensity.shape[0], c_displacedDensity.shape[1]);
    density = yakl_create_real("density", 
                    c_density.shape[0], c_density.shape[1]);
    zMid = yakl_create_real("zMid",
                    c_zMid.shape[0], c_zMid.shape[1]);
    zTop = yakl_create_real("zTop",
                    c_zTop.shape[0], c_zTop.shape[1]);
    RiTopOfCell = yakl_create_real("RiTopOfCell",
                    c_RiTopOfCell.shape[0], c_RiTopOfCell.shape[1]);
    BruntVaisalaFreqTop = yakl_create_real("BruntVaisalaFreqTop", 
                    c_BruntVaisalaFreqTop.shape[0], c_BruntVaisalaFreqTop.shape[1]);

    montgomeryPotential = yakl_create_real("montgomeryPotential", 
                    c_montgomeryPotential.shape[0], c_montgomeryPotential.shape[1]);
    pressure = yakl_create_real("pressure", 
                    c_pressure.shape[0], c_pressure.shape[1]);
    surfacePressure = yakl_create_real("surfacePressure", 
                    c_surfacePressure.shape[0]);

    edgeAreaFractionOfCell = yakl_wrap_array("edgeAreaFractionOfCell",
                            static_cast<double *>(c_edgeAreaFractionOfCell.ptr),
                            c_edgeAreaFractionOfCell.shape[0], c_edgeAreaFractionOfCell.shape[1]);

    normalizedRelativeVorticityEdge = yakl_create_real("normalizedRelativeVorticityEdge", 
                    c_normRelVortEdge.shape[0], c_normRelVortEdge.shape[1]);
    normalizedPlanetaryVorticityEdge = yakl_create_real("normalizedPlanetaryVorticityEdge", 
                    c_normPlanetVortEdge.shape[0], c_normPlanetVortEdge.shape[1]);
    normalizedRelativeVorticityCell = yakl_create_real("normalizedRelativeVorticityCell", 
                    c_normalizedRelativeVorticityCell.shape[0], c_normalizedRelativeVorticityCell.shape[1]);
    
    // scratch arrays
    normalizedPlanetaryVorticityVertex = yakl_create_real("normalizedPlanetaryVorticityVertex", nVertLevels, nVertices);
    normalizedRelativeVorticityVertex = yakl_create_real("normalizedRelativeVorticityVertex", nVertLevels, nVertices);
    vorticityGradientNormalComponent = yakl_create_real("vorticityGradientNormalComponent", nVertLevels, nEdges);
    vorticityGradientTangentialComponent = yakl_create_real("vorticityGradientTangentialComponent", nVertLevels, nEdges);

    ssh = yakl_create_real("ssh", nCells+1);
    seaIcePressure = yakl_create_real("seaIcePressure", nCells);
    atmosphericPressure = yakl_create_real("atmosphericPressure", nCells);

    div_hu = yakl_create_real("div_hu", nVertLevels, nCells);
    div_huTransport = yakl_create_real("div_huTransport", nVertLevels, nCells);
    div_huGMBolus = yakl_create_real("div_huGMBolus", nVertLevels, nCells);
}

extern "C"
void ocn_diagnostics_update_device(double * h_normalVelocity, double * h_layerThickness)
{
    using namespace diag_solve;
    
    //std::cerr << " ocn_diagnostics_update_device..." << std::endl;
    yakl_update_device(normalVelocity, h_normalVelocity, stream1);
    yakl_update_device(layerThickness, h_layerThickness, stream1);
    yakl_record_event(event1, stream1);
    
    yakl_update_device(displacedDensity, static_cast<double *>(c_displacedDensity.ptr), stream1);
    yakl_update_device(density, static_cast<double *>(c_density.ptr), stream1);
    yakl_update_device(BruntVaisalaFreqTop, static_cast<double *>(c_BruntVaisalaFreqTop.ptr), stream1);
    yakl_update_device(zMid, static_cast<double *>(c_zMid.ptr), stream1);
    yakl_update_device(RiTopOfCell, static_cast<double *>(c_RiTopOfCell.ptr), stream1);
    //std::cerr << " ocn_diagnostics_update_device...done" << std::endl;
}

extern "C"
void ocn_diag_solve_circ(int vertexDegree, double * h_normalVelocity, double * h_circulation, 
                        double * h_relativeVorticity)
{
    using diag_solve::stream1;
    using diag_solve::stream2;
    using diag_solve::event1;

    //yakl_update_device(diag_solve::normalVelocity, h_normalVelocity);

    YAKL_LOCAL_NS(diag_solve, circulation);
    YAKL_LOCAL_NS(diag_solve, relativeVorticity);
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(mesh, edgeSignOnVertex);
    YAKL_LOCAL_NS(mesh, minLevelVertexTop);
    YAKL_LOCAL_NS(mesh, maxLevelVertexBot);
    YAKL_LOCAL_NS(mesh, invAreaTriangle);
    //YAKL_LOCAL_NS(mesh, areaTriangle);
    YAKL_LOCAL_NS(mesh, edgesOnVertex);
    YAKL_LOCAL_NS(mesh, dcEdge);
    YAKL_SCOPE(nVertices, diag_solve::nVertices);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    if ( DEBUG ) std::cerr << " ocn_diag_solve_circ. vertexDegree, nVertices, nVertLevels = "
        << vertexDegree << " " << nVertices << " " << nVertLevels << std::endl;

    //std::cerr << "ocn_diag_solve_circ: nVertLevels = " << nVertLevels << std::endl;
    if ( DEBUG) std::cerr << " ocn_diag_solve_circ: loop 1..."  << std::endl;
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nVertices},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iVertex,int k)
    {
         // initialize layerThicknessEdge to avoid divide by zero and NaN problems.
        circulation(k, iVertex) = 0.0;
        relativeVorticity(k, iVertex) = 0.0;
    }, yakl::defaultVectorSize, stream2);

    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diag_solve_circ: loop 2..."  << std::endl;
    }
    yakl_stream_wait(stream2, event1);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nVertices},{1,nVertLevels}) ,
    //yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,1},{1,1}) ,
    YAKL_LAMBDA(int iVertex, int k)
    {
         double invAreaTri1 = invAreaTriangle(iVertex);
         for (int i = 1; i <= vertexDegree; ++i)
         {
            int iEdge = edgesOnVertex(i, iVertex);
            
            if ( (k >= minLevelVertexTop(iVertex)) && (k <= maxLevelVertexBot(iVertex)) )
            {
              double r_tmp = dcEdge(iEdge) * normalVelocity(k, iEdge);
              circulation(k, iVertex) = circulation(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp;
              relativeVorticity(k, iVertex) = relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri1;
            }
         }
    }, yakl::defaultVectorSize, stream2);

    yakl_update_host(diag_solve::circulation, h_circulation, stream2);
    yakl_update_host(diag_solve::relativeVorticity, h_relativeVorticity, stream2);
    
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diag_solve_circ done. " << std::endl;    
    }
}

extern "C"
void ocn_diag_solve_layer_thick_center(int nEdges, double * h_layerThickness,
                                double * h_layerThickEdgeMean)
{
    using diag_solve::stream1;

    if ( DEBUG ) 
    {
        std::cerr << " ocn_diag_solve_layer_thick_center. " << std::endl;    
    }
    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeMean);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        // initialize layerThickEdgeMean to avoid divide by zero and NaN problems.
        layerThickEdgeMean(k, iEdge) = -1.0e34;
    }, yakl::defaultVectorSize, stream1);


    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nEdges}) ,
    YAKL_LAMBDA(int iEdge)
    {
          int cell1 = cellsOnEdge(1,iEdge);
          int cell2 = cellsOnEdge(2,iEdge);
          for (int k = minLevelEdgeBot(iEdge); k <= maxLevelEdgeTop(iEdge); ++k )
          {
            // central differenced
            layerThickEdgeMean(k,iEdge) = 0.5 * (layerThickness(k,cell1) + layerThickness(k,cell2));
          }
    }, yakl::defaultVectorSize, stream1);
    
    yakl_update_host(diag_solve::layerThickEdgeMean, h_layerThickEdgeMean, stream1);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diag_solve_layer_thick_center done. " << std::endl;    
    }
}

extern "C"
void ocn_diag_solve_layer_thick_upwind(int nEdges, double * h_normalVelocity, double * h_layerThickness,
                                double * h_layerThickEdgeFlux)
{
    using diag_solve::stream1;

    if ( DEBUG ) 
    {
        std::cerr << " ocn_diag_solve_layer_thick_upwind. " << std::endl;    
    }
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeFlux);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        // initialize layerThicknessEdge to avoid divide by zero and NaN problems.
        layerThickEdgeFlux(k, iEdge) = -1.0e34;
    }, yakl::defaultVectorSize, stream1);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nEdges}),
    YAKL_LAMBDA(int iEdge)
    {
          int cell1 = cellsOnEdge(1,iEdge);
          int cell2 = cellsOnEdge(2,iEdge);
          for (int k = minLevelEdgeBot(iEdge); k <= maxLevelEdgeTop(iEdge); ++k )
          {
            // upwind
            if (normalVelocity(k, iEdge) > 0.0)
            {
              layerThickEdgeFlux(k,iEdge) = layerThickness(k,cell1);
            }
            else if (normalVelocity(k, iEdge) < 0.0)
            {
              layerThickEdgeFlux(k,iEdge) = layerThickness(k,cell2);
            }
            else
            {
              layerThickEdgeFlux(k,iEdge) = max(layerThickness(k,cell1), layerThickness(k,cell2));
            }
          }
    }, yakl::defaultVectorSize, stream1);

    yakl_update_host(diag_solve::layerThickEdgeFlux, h_layerThickEdgeFlux, stream1);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diag_solve_layer_thick_upwind done. " << std::endl;    
    }
}


extern "C"
void ocn_diagnostic_solve_vort_vel(int nCells, int nCellsH, int nEdges, double * h_normalGMBolusVelocity,
                double * h_normalTransportVelocity,
                double * h_vertVelocityTop, double * h_vertTransportVelocityTop,
                double * h_vertGMBolusVelocityTop, double * h_relativeVorticityCell, 
                double * h_divergence,
                double * h_kineticEnergyCell, double * h_tangentialVelocity)
{
    using diag_solve::stream1;
    using diag_solve::stream2;
    
    if ( DEBUG ) 
    {
        std::cerr << " ocn_diagnostic_solve_vort_vel. " << std::endl;    
    }

    yakl_update_device(diag_solve::normalGMBolusVelocity, h_normalGMBolusVelocity);
    yakl_update_device(diag_solve::normalTransportVelocity, h_normalTransportVelocity);

    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_LOCAL_NS(diag_solve, layerThicknessEdge);
    YAKL_LOCAL_NS(diag_solve, relativeVorticity);
    YAKL_LOCAL_NS(diag_solve, relativeVorticityCell);
    YAKL_LOCAL_NS(diag_solve, divergence);
    YAKL_LOCAL_NS(diag_solve, kineticEnergyCell);
    YAKL_LOCAL_NS(diag_solve, div_hu);
    YAKL_LOCAL_NS(diag_solve, div_huGMBolus);
    YAKL_LOCAL_NS(diag_solve, div_huTransport);
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, normalTransportVelocity);
    YAKL_LOCAL_NS(diag_solve, normalGMBolusVelocity);
    YAKL_LOCAL_NS(diag_solve, vertVelocityTop);
    YAKL_LOCAL_NS(diag_solve, vertTransportVelocityTop);
    YAKL_LOCAL_NS(diag_solve, vertGMBolusVelocityTop);
    YAKL_LOCAL_NS(diag_solve, tangentialVelocity);
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnEdge);
    YAKL_LOCAL_NS(mesh, edgesOnEdge);
    YAKL_LOCAL_NS(mesh, weightsOnEdge);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, invAreaCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, kiteIndexOnCell);
    YAKL_LOCAL_NS(mesh, kiteAreasOnVertex);
    YAKL_LOCAL_NS(mesh, verticesOnCell);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);
    YAKL_LOCAL_NS(mesh, dcEdge);
    YAKL_LOCAL_NS(mesh, dvEdge);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 1. nCells, nVertLevels = "
            << nCells << " " << nVertLevels  << std::endl;
    }
    return;
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
        relativeVorticityCell(k,iCell) = 0.0;
    }, yakl::defaultVectorSize, stream1);
    
    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 2..."  << std::endl;
    }
    int nlvls = tangentialVelocity.extent(0);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nlvls}),
    YAKL_LAMBDA(int iEdge,int k)
    {
        tangentialVelocity(k,iEdge) = 0.0;
    }, yakl::defaultVectorSize, stream2);


    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 3..."  << std::endl;
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nlvls}),
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) )
        {
            double tmp = 0.0;
             for ( int i = 1; i <= nEdgesOnEdge(iEdge); ++i )
             {
                int eoe = edgesOnEdge(i,iEdge);
                tmp += weightsOnEdge(i, iEdge) * normalVelocity(k, eoe);
             }
            tangentialVelocity(k,iEdge) += tmp;
        }
    }, yakl::defaultVectorSize, stream2);

    yakl_update_host(diag_solve::tangentialVelocity, h_tangentialVelocity, stream2);

    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 4..."  << std::endl;
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsH},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
         divergence(k, iCell) = 0.0;
         kineticEnergyCell(k, iCell) = 0.0;
         div_hu(k, iCell) = 0.0;
         div_huTransport(k, iCell) = 0.0;
         div_huGMBolus(k, iCell) = 0.0;
         if ( k < minLevelCell(iCell) )
            vertVelocityTop(k, iCell) = 0.0;
    }, yakl::defaultVectorSize, stream2);
    
    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 5..."  << std::endl;
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell, int k)
    {
        double invAreaCell1 = invAreaCell(iCell);
        for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
          int j = kiteIndexOnCell(i, iCell);
          int iVertex = verticesOnCell(i, iCell);
          if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
          {
            relativeVorticityCell(k, iCell) = relativeVorticityCell(k, iCell) + kiteAreasOnVertex(j, iVertex)
                                            * relativeVorticity(k, iVertex) * invAreaCell1;
          }
        }
    }, yakl::defaultVectorSize, stream1);
    yakl_update_host(diag_solve::relativeVorticityCell, h_relativeVorticityCell, stream1);

    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 6..."  << std::endl;
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsH},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            double invAreaCell1 = invAreaCell(iCell);
            double tdiv = 0.0;
            double tdivhu = 0.0;
            double trke = 0.0;
            double tdivhut = 0.0;
            double tdivhug = 0.0;
            for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
            {
                int iEdge = edgesOnCell(i, iCell);
                double edgeSignOnCell_temp = edgeSignOnCell(i, iCell);
                double dcEdge_temp = dcEdge(iEdge);
                double dvEdge_temp = dvEdge(iEdge);
                double r_tmp = dvEdge_temp * normalVelocity(k, iEdge) * invAreaCell1;
                double h_tmp = layerThicknessEdge(k, iEdge) * edgeSignOnCell_temp
                                * dvEdge_temp * invAreaCell1;

               tdiv += edgeSignOnCell_temp * r_tmp;
               tdivhu += layerThicknessEdge(k, iEdge) * edgeSignOnCell_temp * r_tmp;
               trke += 0.25 * r_tmp * dcEdge_temp * normalVelocity(k,iEdge);
               // Compute vertical velocity from the horizontal total transport
               tdivhut += h_tmp * normalTransportVelocity(k, iEdge);
               // Compute vertical velocity from the horizontal GM Bolus velocity
               tdivhug += h_tmp * normalGMBolusVelocity(k, iEdge);
            }
            
            divergence(k, iCell) = divergence(k, iCell) - tdiv;
            div_hu(k, iCell) = div_hu(k, iCell) - tdivhu;
            kineticEnergyCell(k, iCell) += trke;
            // Compute vertical velocity from the horizontal total transport
            div_huTransport(k, iCell) = div_huTransport(k, iCell) - tdivhut;
            // Compute vertical velocity from the horizontal GM Bolus velocity
            div_huGMBolus(k, iCell) = div_huGMBolus(k, iCell) - tdivhug;
    }
    }, yakl::defaultVectorSize, stream2);
    yakl_update_host(diag_solve::kineticEnergyCell, h_kineticEnergyCell, stream2);
    yakl_update_host(diag_solve::divergence, h_divergence, stream2);

    if ( DEBUG) {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel: loop 7..."  << std::endl;
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nCellsH}),
    YAKL_LAMBDA(int iCell)
    {
        // Vertical velocity at bottom (maxLevelCell(iCell)+1) is zero, initialized above.
        vertVelocityTop(maxLevelCell(iCell)+1, iCell) = 0.0;
        vertTransportVelocityTop(maxLevelCell(iCell)+1, iCell) = 0.0;
        vertGMBolusVelocityTop(maxLevelCell(iCell)+1, iCell) = 0.0;
        // TO DO: this is essentially a reverse inclusive scan.  need to find
        // existing GPU algorithm.
        for (int k = maxLevelCell(iCell); k >= 1; --k )
        {
            vertVelocityTop(k,iCell) = vertVelocityTop(k+1,iCell) - div_hu(k,iCell);
            vertTransportVelocityTop(k,iCell) = vertTransportVelocityTop(k+1,iCell) - div_huTransport(k,iCell);
            vertGMBolusVelocityTop(k,iCell) = vertGMBolusVelocityTop(k+1,iCell) - div_huGMBolus(k,iCell);
        }
    }, yakl::defaultVectorSize, stream2);
        
    yakl_update_host(diag_solve::vertVelocityTop, h_vertVelocityTop, stream2);
    yakl_update_host(diag_solve::vertTransportVelocityTop, h_vertTransportVelocityTop, stream2);
    yakl_update_host(diag_solve::vertGMBolusVelocityTop, h_vertGMBolusVelocityTop, stream2);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort_vel done. " << std::endl;    
    }
}

extern "C"
void ocn_diagnostic_solve_vort(int nVerticesH, int nEdgesH, int nCellsH, int vertexDegree,
                double * h_normalizedRelativeVorticityEdge,
                double * h_normalizedPlanetaryVorticityEdge,
                double * h_normalizedRelativeVorticityCell)
{
    using diag_solve::stream1;
    using diag_solve::stream2;

    if ( DEBUG ) 
    {
        std::cerr << " ocn_diagnostic_solve_vort. " << std::endl;    
    }
    YAKL_LOCAL_NS(mesh, invAreaTriangle);
    YAKL_LOCAL_NS(mesh, minLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, maxLevelVertexBot);
    YAKL_LOCAL_NS(mesh, kiteAreasOnVertex);
    YAKL_LOCAL_NS(mesh, cellsOnVertex);
    YAKL_LOCAL_NS(mesh, verticesOnEdge);
    YAKL_LOCAL_NS(mesh, fVertex);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, kiteIndexOnCell);
    YAKL_LOCAL_NS(mesh, verticesOnCell);
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, invAreaCell);
    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_LOCAL_NS(diag_solve, relativeVorticity);
    YAKL_LOCAL_NS(diag_solve, normalizedRelativeVorticityEdge);
    YAKL_LOCAL_NS(diag_solve, normalizedPlanetaryVorticityEdge);
    YAKL_LOCAL_NS(diag_solve, normalizedRelativeVorticityVertex);
    YAKL_LOCAL_NS(diag_solve, normalizedPlanetaryVorticityVertex);
    YAKL_LOCAL_NS(diag_solve, normalizedRelativeVorticityCell);
    YAKL_SCOPE(nCells, diag_solve::nCells);
    YAKL_SCOPE(nEdges, diag_solve::nEdges);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);
    YAKL_SCOPE(nVertices, diag_solve::nVertices);

    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort loop 1. " << std::endl;    
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}),
    YAKL_LAMBDA(int iEdge,int k)
    {
        normalizedRelativeVorticityEdge(k, iEdge) = 0.0;
        normalizedPlanetaryVorticityEdge(k, iEdge) = 0.0;
    });

    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort loop 2. " << std::endl;    
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
        normalizedRelativeVorticityCell(k, iCell) = 0.0;
    });

    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort loop 3. nVerticesH, nVertLevels, vertexDegree = "
            <<  nVerticesH << " " <<  nVertLevels << " " << vertexDegree << std::endl;    
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nVerticesH},{1,nVertLevels}),
    YAKL_LAMBDA(int iVertex,int k)
    {
        double invAreaTri1 = invAreaTriangle(iVertex);
        if ( k <= maxLevelVertexBot(iVertex) )
        {
            double layerThicknessVertex = 0.0;
            for ( int i = 1; i <= vertexDegree; ++i )
               layerThicknessVertex = layerThicknessVertex + layerThickness(k,cellsOnVertex(i,iVertex))
                                    * kiteAreasOnVertex(i,iVertex);

            layerThicknessVertex = layerThicknessVertex * invAreaTri1;
            if (layerThicknessVertex != 0)
            {
                layerThicknessVertex = 1.0 / layerThicknessVertex;
                normalizedRelativeVorticityVertex(k,iVertex) = relativeVorticity(k,iVertex) * layerThicknessVertex;
                normalizedPlanetaryVorticityVertex(k,iVertex) = fVertex(iVertex) * layerThicknessVertex;
            }
        }
    }, yakl::defaultVectorSize, stream1);

    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort loop 4. " << std::endl;    
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesH},{1,nVertLevels}),
    YAKL_LAMBDA(int iEdge,int k)
    {
        int vertex1 = verticesOnEdge(1, iEdge);
        int vertex2 = verticesOnEdge(2, iEdge);
        if ( (k >= minLevelEdgeTop(iEdge)) && (k <= maxLevelEdgeBot(iEdge)) )
        {
          normalizedRelativeVorticityEdge(k, iEdge) = 0.5 * (normalizedRelativeVorticityVertex(k, vertex1)
                                                    + normalizedRelativeVorticityVertex(k, vertex2));
          normalizedPlanetaryVorticityEdge(k, iEdge) = 0.5 * (normalizedPlanetaryVorticityVertex(k, vertex1)
                                                     + normalizedPlanetaryVorticityVertex(k, vertex2));
        }
    }, yakl::defaultVectorSize, stream1);

    yakl_update_host(diag_solve::normalizedRelativeVorticityEdge, h_normalizedRelativeVorticityEdge);
    yakl_update_host(diag_solve::normalizedPlanetaryVorticityEdge, h_normalizedPlanetaryVorticityEdge);

    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort loop 5. " << std::endl;    
    }
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsH},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
        double invAreaCell1 = invAreaCell(iCell);

        for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
          int j = kiteIndexOnCell(i, iCell);
          int iVertex = verticesOnCell(i, iCell);
          if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
          {
            normalizedRelativeVorticityCell(k, iCell) = normalizedRelativeVorticityCell(k, iCell)
              + kiteAreasOnVertex(j, iVertex) * normalizedRelativeVorticityVertex(k, iVertex) * invAreaCell1;
          }
        }
    }, yakl::defaultVectorSize, stream1);
    
    yakl_update_host(diag_solve::normalizedRelativeVorticityCell, h_normalizedRelativeVorticityCell);

    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_vort done. " << std::endl;    
    }
}

extern "C"
void ocn_diagnostic_apvm(int nEdgesH,
                double * h_normalizedRelativeVorticityEdge)
{

}

extern "C"
void ocn_diagnostic_solve_zcoord(double * h_ssh)
{
    if ( DEBUG ) 
    {
        std::cerr << " ocn_diagnostic_solve_zcoord. " << std::endl;    
    }
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, bottomDepth);
    YAKL_LOCAL_NS(diag_solve, zMid);
    YAKL_LOCAL_NS(diag_solve, zTop);
    YAKL_LOCAL_NS(diag_solve, ssh);
    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_SCOPE(nCells, diag_solve::nCells);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    using diag_solve::stream1;

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
        if ( k >= maxLevelCell(iCell) )
        {
            double fixed = layerThickness(maxLevelCell(iCell),iCell);
            zMid(k,iCell) = -bottomDepth(iCell) + 0.5*fixed;
            zTop(k,iCell) = -bottomDepth(iCell) + fixed;
        }
    }, yakl::defaultVectorSize, stream1);

    // Compute zMid, the z-coordinate of the middle of the layer.
    // Compute zTop, the z-coordinate of the top of the layer.
    // Note the negative sign, since bottomDepth is positive
    // and z-coordinates are negative below the surface.
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nCells}),
    YAKL_LAMBDA(int iCell)
    {
        for ( int k = maxLevelCell(iCell)-1; k >= minLevelCell(iCell); --k )
        {
            zMid(k,iCell) = zMid(k+1,iCell)
                            + 0.5*(  layerThickness(k+1,iCell)
                            + layerThickness(k  ,iCell));
            zTop(k,iCell) = zTop(k+1,iCell)
              + layerThickness(k  ,iCell);
        }
        // copy zTop(1,iCell) into sea-surface height array
        ssh(iCell) = zTop(minLevelCell(iCell),iCell);
    }, yakl::defaultVectorSize, stream1);
    yakl_update_host(diag_solve::zMid, GET_DPTR(zMid), stream1);
    yakl_update_host(diag_solve::zTop, GET_DPTR(zTop), stream1);
    yakl_update_host(diag_solve::ssh, h_ssh, stream1);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_zcoord done. " << std::endl;    
    }
}


extern "C"
void ocn_diagnostic_solve_rich(double coef, int nCellsH,
                double * h_BruntVaisalaFreqTop,
                double * h_RiTopOfCell)
{
    if ( DEBUG ) 
    {
        std::cerr << " ocn_diagnostic_solve_rich. " << std::endl;    
    }
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(diag_solve, edgeAreaFractionOfCell);
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, displacedDensity);
    YAKL_LOCAL_NS(diag_solve, density);
    YAKL_LOCAL_NS(diag_solve, zMid);
    YAKL_LOCAL_NS(diag_solve, BruntVaisalaFreqTop);
    YAKL_LOCAL_NS(diag_solve, RiTopOfCell);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);

    YAKL_SCOPE(nCells, diag_solve::nCells);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    using diag_solve::stream1;

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
         RiTopOfCell(k,iCell) = 100.0;
    }, yakl::defaultVectorSize, stream1);
    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
         if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
         {
            if ( k == minLevelCell(iCell) )
            {
                BruntVaisalaFreqTop(minLevelCell(iCell),iCell) = 0.0;
            }
            else
            {
                BruntVaisalaFreqTop(k,iCell) = coef * (displacedDensity(k-1,iCell) - density(k,iCell))
                    / (zMid(k-1,iCell) - zMid(k,iCell));
            }
         }
    }, yakl::defaultVectorSize, stream1);

    yakl_update_host(diag_solve::BruntVaisalaFreqTop, 
                     h_BruntVaisalaFreqTop, stream1);


    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsH},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell,int k)
    {
         if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
         {
            if ( k == minLevelCell(iCell) )
            {
             double shearSquared = 0.0;
             int kp = minLevelCell(iCell) + 1;
             for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
             {
                int iEdge = edgesOnCell(i, iCell);
                double delU2 = (normalVelocity(kp-1,iEdge) - normalVelocity(kp,iEdge));
                delU2 *= delU2;
                shearSquared = shearSquared + edgeAreaFractionOfCell(i,iCell) * delU2;
             }
             //  Note that the factor of two is from averaging dot product to cell center on a C-grid
             double shearMean = sqrt(2.0*shearSquared );
             shearMean = shearMean / (zMid(kp-1,iCell) - zMid(kp,iCell));
             RiTopOfCell(k,iCell) = BruntVaisalaFreqTop(kp,iCell) / (shearMean*shearMean + 1.0e-10);            
            }
            else
            {
            
             double shearSquared = 0.0;
             for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
             {
                int iEdge = edgesOnCell(i, iCell);
                double delU2 = (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge));
                delU2 *= delU2;
                shearSquared = shearSquared + edgeAreaFractionOfCell(i,iCell) * delU2;
             }
             //  Note that the factor of two is from averaging dot product to cell center on a C-grid
             double shearMean = sqrt(2.0*shearSquared );
             shearMean = shearMean / (zMid(k-1,iCell) - zMid(k,iCell));
             RiTopOfCell(k,iCell) = BruntVaisalaFreqTop(k,iCell) / (shearMean*shearMean + 1.0e-10);
            }
         }
    }, yakl::defaultVectorSize, stream1);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nCellsH}),
    YAKL_LAMBDA(int iCell)
    {
        int k = minLevelCell(iCell);
          RiTopOfCell(k,iCell) = RiTopOfCell(k+1,iCell);
    }, yakl::defaultVectorSize, stream1);
    
    yakl_update_host(diag_solve::RiTopOfCell, h_RiTopOfCell, stream1);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_solve_rich done. " << std::endl;    
    }
}

extern "C"
void ocn_diagnostic_montgomery()
{
    if ( DEBUG ) 
    {
        std::cerr << " ocn_diagnostic_montgomery. " << std::endl;    
    }
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(diag_solve, montgomeryPotential);
    YAKL_LOCAL_NS(diag_solve, RiTopOfCell);
    YAKL_LOCAL_NS(diag_solve, pressure);

    YAKL_SCOPE(nCells, diag_solve::nCells);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    using diag_solve::stream1;

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell, int k)
    {
        double ptop = 0.0;
        for ( int i = 2; i <= k; ++i )
        int k = minLevelCell(iCell);
          RiTopOfCell(k,iCell) = RiTopOfCell(k+1,iCell);
    }, yakl::defaultVectorSize, stream1);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diagnostic_montgomery done. " << std::endl;    
    }
}

extern "C"
void ocn_diag_update_density()
{
    using diag_solve::stream1;
    yakl_update_device(diag_solve::density, GET_DPTR(density), stream1);
}

extern "C"
void ocn_diag_pressure(double gravity)
{
    if ( DEBUG ) 
    {
        std::cerr << " ocn_diag_pressure. " << std::endl;    
    }
    yakl_update_device(diag_solve::surfacePressure, GET_DPTR(surfacePressure));

    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(diag_solve, density);
    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_LOCAL_NS(diag_solve, pressure);
    YAKL_LOCAL_NS(diag_solve, surfacePressure);

    YAKL_SCOPE(nCells, diag_solve::nCells);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);

    using diag_solve::stream1;

    double coef = gravity*0.5;
    //yakl::timer_start("ocn_diag_pressure");
    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell, int k)
    {
          if ( k == minLevelCell(iCell) )
          {
            pressure(k,iCell) = surfacePressure(iCell)
                + density(k,iCell)*coef*layerThickness(k,iCell);
          }
          else if ( k <= maxLevelCell(iCell) )
          {
            pressure(k,iCell) = surfacePressure(iCell);
            double tmp = 0.0;
            for ( int j = minLevelCell(iCell); j <= k-1; ++j )
                tmp = tmp + density(j,iCell)*layerThickness(j,iCell);
            pressure(k,iCell) = pressure(k,iCell) +
                    coef * (2*tmp + density(k  ,iCell)*layerThickness(k  ,iCell));
          }
    }, yakl::defaultVectorSize, stream1);

    /*
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nCells}),
    YAKL_LAMBDA(int iCell)
    {
       pressure(minLevelCell(iCell),iCell) = surfacePressure(iCell)
         + density(minLevelCell(iCell),iCell)*coef*layerThickness(minLevelCell(iCell),iCell);

       for (int k = minLevelCell(iCell)+1; k<= maxLevelCell(iCell); ++k )
          pressure(k,iCell) = pressure(k-1,iCell)
            + coef*(  density(k-1,iCell)*layerThickness(k-1,iCell)
                                 + density(k  ,iCell)*layerThickness(k  ,iCell));
    });
    */
    //yakl::timer_stop("ocn_diag_pressure");

    
    yakl_update_host(diag_solve::pressure, GET_DPTR(pressure), stream1);
    if ( DEBUG ) 
    {
        yakl::fence();
        std::cerr << " ocn_diag_pressure done. " << std::endl;    
    }
}

