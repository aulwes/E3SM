#include "mpas_ocn_diagnostics_variables.hxx"
#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_mesh_c.hxx"
#include "mpas_ocn_yakl_c.hxx"

/***
 mpas_ocn_mesh_c.hxx:
   advCoefs, advCoefs3rd,nAdvCellsForEdge,advCellsForEdge
**/


namespace
{

// This should match values in mpas_ocn_tracer_advection_mono.F
enum {
      vertOrder2 = 2,    //< 2nd order scheme
      vertOrder3 = 3,    //< 3rd order scheme
      vertOrder4 = 4     //< 4th order scheme
} vertOrder;

}

namespace tracer_mono
{

yakl::yakl_stream_t     stream = 0;

double          coef3rdOrder;
double          eps;

d_double_2d_t   * hProvInv, 
                * flxIn, 
                * flxOut,
                * highOrderFlx,
                * hNewInv,
                * hProv,
                * lowOrderFlx,
                * tracerCur,
                * tracerMax,
                * tracerMin,
                * workTend,
                * layerThickness,
                * normalThicknessFlux,
                * wgtTmp,
                * sgnTmp,
                * w
                ;
                
d_double_3d_t   * tend
                ;

};

extern "C"
void ocn_advect_mono_copyin(int * tendSize, double * h_normalThicknessFlux, 
                                double * h_w, double * h_tend, double * h_layerThickness)
{
    using namespace tracer_mono;
    
    //std::cerr << " tend size = " << tendSize[0] << " " << tendSize[1] << " " << tendSize[2] << std::endl;
    tend = yakl_wrap_array("tend", h_tend, tendSize[0], tendSize[1], tendSize[2]);
    yakl_update_device(normalThicknessFlux, h_normalThicknessFlux, stream);
    yakl_update_device(w, h_w, stream);
    yakl_update_device(layerThickness, h_layerThickness, stream);
}

extern "C"
void ocn_advect_mono_end(double * h_tend)
{
    yakl_update_host(tracer_mono::tend, h_tend, tracer_mono::stream);
    yakl::fence();

    delete tracer_mono::tend;
} 

extern "C"
void ocn_advect_mono_cell_init(int nCells)
{
    YAKL_LOCAL_NS(tracer_mono,highOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,workTend);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels+1}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        highOrderFlx(k, iCell) = 0.0;
        if ( k <= nVertLevels ) 
               workTend(k, iCell) = 0.0;
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
}

extern "C"
void ocn_advect_mono_check(double * h_advCoefs, 
                            double * h_dvEdge, 
                            double * h_advCoefs3rd)
{
    yakl_update_host(mesh::dvEdge, h_dvEdge);
    yakl_update_host(mesh::advCoefs, h_advCoefs);
    yakl_update_host(mesh::advCoefs3rd, h_advCoefs3rd);
}

extern "C"
void ocn_advect_mono_check2(int * h_nAdvCellsForEdge, 
                            int * h_advCellsForEdge)
{
    yakl_update_host(mesh::nAdvCellsForEdge, h_nAdvCellsForEdge);
    yakl_update_host(mesh::advCellsForEdge, h_advCellsForEdge);
}


extern "C"
void ocn_advect_mono_pre(double dt, double * h_hProv, double * h_hProvInv, double * h_hNewInv)
{
    YAKL_LOCAL_NS(tracer_mono,hProv);
    YAKL_LOCAL_NS(tracer_mono,hProvInv);
    YAKL_LOCAL_NS(tracer_mono,hNewInv);
    YAKL_LOCAL_NS(tracer_mono,layerThickness);
    YAKL_LOCAL_NS(tracer_mono,normalThicknessFlux);
    YAKL_LOCAL_NS(tracer_mono,w);

    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(mesh, dvEdge);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);
    YAKL_LOCAL_NS(mesh, invAreaCell);

    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    int nCells = invAreaCell.extent(0);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        double invAreaCell1 = dt * invAreaCell(iCell);
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            hProv(k, iCell) = layerThickness(k, iCell);
            for (int i = 1; i <= nEdgesOnCell(iCell); ++i )
            {
              int iEdge = edgesOnCell(i,iCell);
              double signedFactor = invAreaCell1*dvEdge(iEdge)*edgeSignOnCell(i,iCell);
              // Provisional layer thickness is after horizontal thickness flux only          
              hProv(k, iCell) = hProv(k, iCell)
                                + signedFactor*normalThicknessFlux(k,iEdge);
            }
            hProvInv(k,iCell) = 1.0 / hProv(k,iCell);
            hNewInv (k,iCell) = 1.0 /(hProv(k,iCell) -
                                         dt*w(k,iCell) + dt*w(k+1, iCell));
        }
        /*
        for (int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
          int iEdge = edgesOnCell(i,iCell);
          double signedFactor = invAreaCell1*dvEdge(iEdge)*edgeSignOnCell(i,iCell);
          // Provisional layer thickness is after horizontal thickness flux only          
          for (int k = 1; k <= kmax; ++k ) 
          {
            hProv(k, iCell) = hProv(k, iCell)
                            + signedFactor*normalThicknessFlux(k,iEdge);
          }
        }
        
        // New layer thickness is after horizontal and vertical thickness flux
        for (int k = 1; k <= kmax; ++k ) 
        {
          hProvInv(k,iCell) = 1.0 / hProv(k,iCell);
          hNewInv (k,iCell) = 1.0 /(hProv(k,iCell) -
                                         dt*w(k,iCell) + dt*w(k+1, iCell));
        }
        */
        
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
}

extern "C"
//void ocn_advect_mono_set_minmax(int nCells, double * h_tracerCur, double * h_tracerMin,
void ocn_advect_mono_set_flux(int nCells, int nEdges, double * h_tracerCur, double * h_tracerMin,
                                double * h_tracerMax, double * h_lowOrderFlx, double * h_highOrderFlx)
{
    YAKL_LOCAL_NS(tracer_mono,tracerCur);
    YAKL_LOCAL_NS(tracer_mono,tracerMin);
    YAKL_LOCAL_NS(tracer_mono,tracerMax);
    YAKL_LOCAL_NS(tracer_mono,lowOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,highOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,normalThicknessFlux);
    YAKL_LOCAL_NS(tracer_mono,wgtTmp);
    YAKL_LOCAL_NS(tracer_mono,sgnTmp);
    YAKL_SCOPE(coef3rdOrder, tracer_mono::coef3rdOrder);

    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, cellsOnCell);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, highOrderAdvectionMask);
    YAKL_LOCAL_NS(mesh, advCellsForEdge);
    YAKL_LOCAL_NS(mesh, nAdvCellsForEdge);
    YAKL_LOCAL_NS(mesh, advCoefs);
    YAKL_LOCAL_NS(mesh, advCoefs3rd);
    YAKL_LOCAL_NS(mesh, dvEdge);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    yakl_update_device(tracer_mono::tracerCur, h_tracerCur, tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        if ( (k >= minLevelCell(iCell)) &&  (k <= maxLevelCell(iCell)) )
        {
            tracerMin(k,iCell) = tracerCur(k,iCell);
            tracerMax(k,iCell) = tracerCur(k,iCell);
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        for (int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
            int cell2 = cellsOnCell(i,iCell);
            int kmin  = max(minLevelCell(iCell),
                        minLevelCell(cell2));
            int kmax  = std::min(maxLevelCell(iCell), maxLevelCell(cell2));
            if ( (k >= kmin) && (k <= kmax) )
            {
              tracerMax(k,iCell) = std::max(tracerMax(k,iCell), tracerCur(k,cell2));
              tracerMin(k,iCell) = std::min(tracerMin(k,iCell), tracerCur(k,cell2));
            }
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        // compute some common intermediate factors
        wgtTmp(k,iEdge) = normalThicknessFlux(k,iEdge)*highOrderAdvectionMask(k,iEdge);
        sgnTmp(k,iEdge) = SIGN(1.0, normalThicknessFlux(k,iEdge));
        highOrderFlx(k,iEdge) = 0.0;
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        for ( int i = 1; i <= nAdvCellsForEdge(iEdge); ++i )
        {
            int iCell = advCellsForEdge(i,iEdge);
            double coef1 = advCoefs       (i,iEdge);
            double coef3 = advCoefs3rd    (i,iEdge)*coef3rdOrder;
            if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
            {
               highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)*
                          wgtTmp(k,iEdge)*(coef1 + coef3*sgnTmp(k,iEdge));
            }
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) &&  (k<= maxLevelEdgeTop(iEdge)) )
        {
            int cell1 = cellsOnEdge(1, iEdge);
            int cell2 = cellsOnEdge(2, iEdge);

            // Compute 2nd order fluxes where needed.
            // Also compute low order upwind horizontal flux (monotonic and diffused)
            // Remove low order flux from the high order flux
            // Store left over high order flux in highOrderFlx array
            double tracerWeight = ((int(highOrderAdvectionMask(k, iEdge))+1) & 1)
                         * (dvEdge(iEdge) * 0.5)
                         * normalThicknessFlux(k, iEdge);

            lowOrderFlx(k,iEdge) = dvEdge(iEdge) * 
               (std::max(0.0,normalThicknessFlux(k,iEdge))*tracerCur(k,cell1)
              + std::min(0.0,normalThicknessFlux(k,iEdge))*tracerCur(k,cell2));

            highOrderFlx(k, iEdge) = highOrderFlx(k, iEdge)
                                   + tracerWeight * (tracerCur(k, cell1)
                                                   + tracerCur(k, cell2));

            highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge)
                                  -  lowOrderFlx(k,iEdge);
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
}

extern "C"
void ocn_advect_mono_set_flux_inout(int initNCells, int nCells, int nEdges, double dt,
                                double * h_highOrderFlx, double * h_workTend, double * h_flxIn, double * h_flxOut)
{
    YAKL_LOCAL_NS(tracer_mono,workTend);
    YAKL_LOCAL_NS(tracer_mono,flxIn);
    YAKL_LOCAL_NS(tracer_mono,flxOut);
    YAKL_LOCAL_NS(tracer_mono,lowOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,highOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,tracerCur);
    YAKL_LOCAL_NS(tracer_mono,hProvInv);
    YAKL_LOCAL_NS(tracer_mono,layerThickness);
    YAKL_LOCAL_NS(tracer_mono,tracerMax);
    YAKL_LOCAL_NS(tracer_mono,tracerMin);

    YAKL_LOCAL_NS(mesh, invAreaCell);
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, cellsOnCell);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);
    YAKL_LOCAL_NS(mesh, areaCell);
    YAKL_SCOPE(eps, tracer_mono::eps);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,initNCells+1},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        workTend(k, iCell) = 0.0;
        flxIn   (k, iCell) = 0.0;
        flxOut  (k, iCell) = 0.0;
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        double invAreaCell1 = invAreaCell(iCell);

        // Finish computing the low order horizontal fluxes
        // Upwind fluxes are accumulated in workTend
        for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
            int iEdge = edgesOnCell(i, iCell);
            if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) )
            {
              int cell1 = cellsOnEdge(1,iEdge);
              int cell2 = cellsOnEdge(2,iEdge);
              double signedFactor = edgeSignOnCell(i, iCell) * invAreaCell1;
            
              // Here workTend is the advection tendency due to the
              // upwind (low order) fluxes.
              workTend(k,iCell) = workTend(k,iCell)
                                + signedFactor*lowOrderFlx(k,iEdge);

              // Accumulate remaining high order fluxes
              flxOut(k,iCell) = flxOut(k,iCell) + min(0.0,
                                signedFactor*highOrderFlx(k,iEdge));
              flxIn (k,iCell) = flxIn (k,iCell) + max(0.0,
                                signedFactor*highOrderFlx(k,iEdge));
            }
       }
    });

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
       // Build the factors for the FCT
       // Computed using the bounds that were computed previously,
       // and the bounds on the newly updated value
       // Factors are placed in the flxIn and flxOut arrays
       if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
       {
            // Here workTend is the upwind tendency
            double tracerUpwindNew = (tracerCur(k,iCell)*layerThickness(k,iCell)
                             + dt*workTend(k,iCell)) * hProvInv(k,iCell);
            double tracerMinNew = tracerUpwindNew
                         + dt*flxOut(k,iCell)*hProvInv(k,iCell);
            double tracerMaxNew = tracerUpwindNew
                         + dt*flxIn (k,iCell)*hProvInv(k,iCell);

            double scaleFactor = (tracerMax(k,iCell) - tracerUpwindNew)/
                          (tracerMaxNew - tracerUpwindNew + eps);
            flxIn (k,iCell) = std::min(1.0, std::max(0.0, scaleFactor));

            scaleFactor = (tracerUpwindNew - tracerMin(k,iCell))/
                          (tracerUpwindNew - tracerMinNew + eps);
            flxOut(k,iCell) = std::min(1.0, std::max(0.0, scaleFactor));
       }
    });

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) )
        {
            int cell1 = cellsOnEdge(1, iEdge);
            int cell2 = cellsOnEdge(2, iEdge);
            highOrderFlx(k,iEdge) = max(0.0,highOrderFlx(k,iEdge))*
                                    min(flxOut(k,cell1), flxIn (k,cell2))
                                  + min(0.0,highOrderFlx(k,iEdge))*
                                    min(flxIn (k,cell1), flxOut(k,cell2));
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
}


extern "C"
void ocn_advect_mono_compute_tend(int iTracer, int nCells, double dt, 
                                double * h_tend)
{
    using yakl::COLON;
    
    YAKL_LOCAL_NS(tracer_mono,workTend);
    YAKL_LOCAL_NS(tracer_mono,tend);
    YAKL_LOCAL_NS(tracer_mono,highOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,tracerCur);
    YAKL_LOCAL_NS(tracer_mono,hProvInv);
    YAKL_LOCAL_NS(tracer_mono,layerThickness);

    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);
    YAKL_LOCAL_NS(mesh, invAreaCell);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    auto tslice = tend.slice<2>(COLON, COLON, iTracer);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell, int k)
    {
        double invAreaCell1 = invAreaCell(iCell);
        
        // Accumulate the scaled high order horizontal tendencies
        for ( int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
            int iEdge = edgesOnCell(i, iCell);
            double signedFactor = invAreaCell1 * edgeSignOnCell(i, iCell);
            if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) )
            {
              // workTend on the RHS is the upwind tendency
              // workTend on the LHS is the total horizontal advection tendency
              workTend(k,iCell) = workTend(k,iCell)
                                + signedFactor * highOrderFlx(k, iEdge);
            }
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
        
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell, int k)
    {
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            // workTend on the RHS is the total horizontal advection tendency
            // tracerCur on LHS is the  provisional tracer after horizontal fluxes only.
            tracerCur(k,iCell) = (tracerCur(k,iCell)*layerThickness(k,iCell)
                                  + dt*workTend(k,iCell))*hProvInv(k,iCell);
            tslice(k,iCell) = tslice(k,iCell) + workTend(k,iCell);
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
}

extern "C"
void ocn_advect_mono_vert_flux(int nCells, int mxlvl, int vertOrder,
                                double * h_flxIn, double * h_flxOut,
                                double * h_tracerCur, double * h_workTend)
{
    YAKL_LOCAL_NS(tracer_mono,highOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,lowOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,flxIn);
    YAKL_LOCAL_NS(tracer_mono,flxOut);
    YAKL_LOCAL_NS(tracer_mono,workTend);
    YAKL_LOCAL_NS(tracer_mono,tracerCur);
    YAKL_LOCAL_NS(tracer_mono,tracerMin);
    YAKL_LOCAL_NS(tracer_mono,tracerMax);
    YAKL_LOCAL_NS(tracer_mono,hProv);
    YAKL_LOCAL_NS(tracer_mono,w);

    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);
    YAKL_LOCAL_NS(mesh, areaCell);
    YAKL_SCOPE(coef3rdOrder, tracer_mono::coef3rdOrder);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCells) ,
    YAKL_LAMBDA(int iCell)
    {
        int kmax = maxLevelCell(iCell);

        // take care of top cell
        tracerMax(1,iCell) = std::max(tracerCur(1,iCell),
                                   tracerCur(2,iCell));
        tracerMin(1,iCell) = std::min(tracerCur(1,iCell),
                                   tracerCur(2,iCell));
        for ( int k=2; k <= kmax-1; ++k )
        {
            tracerMax(k,iCell) = std::max(std::max(tracerCur(k-1,iCell),
                                     tracerCur(k  ,iCell)),
                                     tracerCur(k+1,iCell));
            tracerMin(k,iCell) = std::min(std::min(tracerCur(k-1,iCell),
                                     tracerCur(k  ,iCell)),
                                     tracerCur(k+1,iCell));
        }
        // finish with bottom cell
        tracerMax(kmax,iCell) = std::max(tracerCur(kmax  ,iCell),
                                      tracerCur(kmax-1,iCell));
        tracerMin(kmax,iCell) = std::min(tracerCur(kmax  ,iCell),
                                      tracerCur(kmax-1,iCell));
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    double inv12 = 1.0 / 12.0;

    switch ( vertOrder )
    {
    case vertOrder4:
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,mxlvl}),
        YAKL_LAMBDA(int iCell, int k)
        {
            if ( (k >= minLevelCell(iCell)+2) && (k < maxLevelCell(iCell)) )
                highOrderFlx(k, iCell) = w(k,iCell)*(
                     7.0*(tracerCur(k  ,iCell) + tracerCur(k-1,iCell)) -
                               (tracerCur(k+1,iCell) + tracerCur(k-2,iCell))) * inv12;
        }, yakl::LaunchConfig<>(), tracer_mono::stream);
        break;
    case vertOrder3:
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,mxlvl}),
        YAKL_LAMBDA(int iCell, int k)
        {
            if ( (k >= minLevelCell(iCell)+2) && (k < maxLevelCell(iCell)) )
                highOrderFlx(k, iCell) = (w(k,iCell)*
                      (7.0*(tracerCur(k  ,iCell) + tracerCur(k-1,iCell)) -
                                 (tracerCur(k+1,iCell) + tracerCur(k-2,iCell))) -
                        coef3rdOrder*abs(w(k,iCell))*
                                ((tracerCur(k+1,iCell) - tracerCur(k-2,iCell)) -
                       3.0*(tracerCur(k  ,iCell) - tracerCur(k-1,iCell)))) * inv12;
        }, yakl::LaunchConfig<>(), tracer_mono::stream);
        break;
    case vertOrder2:
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{3,mxlvl}) ,
        YAKL_LAMBDA(int iCell, int k)
        {
            if ( (k >= minLevelCell(iCell)+2) && (k < maxLevelCell(iCell)) )
            {
              double verticalWeightK   = hProv(k-1,iCell) /
                                 (hProv(k  ,iCell) + hProv(k-1,iCell));
              double verticalWeightKm1 = hProv(k  ,iCell) /
                                 (hProv(k  ,iCell) + hProv(k-1,iCell));
              highOrderFlx(k,iCell) = w(k,iCell)*
                                 (verticalWeightK  *tracerCur(k  ,iCell) +
                                  verticalWeightKm1*tracerCur(k-1,iCell));
            }
        }, yakl::LaunchConfig<>(), tracer_mono::stream);
        break;
    }

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        int kmin = minLevelCell(iCell);
        int kmax = maxLevelCell(iCell);
        // Next-to-top cell in column is second-order
        if ( k == kmin )
        {
            highOrderFlx(k,iCell) = 0.0;
            lowOrderFlx(k,iCell) = 0.0;
        }

        if ( (kmax > kmin) && (k == kmin+1) )
        {
            double verticalWeightK   = hProv(k-1,iCell) /
                               (hProv(k,iCell) + hProv(k-1, iCell));
            double verticalWeightKm1 = hProv(k,iCell) /
                               (hProv(k,iCell) + hProv(k-1, iCell));
            highOrderFlx(k,iCell) = w(k,iCell)*
                               (verticalWeightK  *tracerCur(k,iCell) +
                                verticalWeightKm1*tracerCur(k-1,iCell));
        }

        // Deepest vertical cell in column is second order
        int km = max(kmin+1,kmax);
        if ( k == km )
        {
            double verticalWeightK   = hProv(k-1,iCell) /
                                 (hProv(k  ,iCell) + hProv(k-1,iCell));
            double verticalWeightKm1 = hProv(k  ,iCell) /
                                 (hProv(k  ,iCell) + hProv(k-1,iCell));
            highOrderFlx(k,iCell) = w(k,iCell)*
                                 (verticalWeightK  *tracerCur(k  ,iCell) +
                                  verticalWeightKm1*tracerCur(k-1,iCell));
        }
        else if ( k == km + 1 )
            highOrderFlx(k+1,iCell) = 0.0;

        if ( k == kmax)
            lowOrderFlx(k+1,iCell) = 0.0;

        // Compute low order (upwind) flux and remove from high order
        if ( (k >= kmin+1) && (k <= kmax) )
        {
            lowOrderFlx(k,iCell) =
                std::min(0.0,w(k,iCell))*tracerCur(k-1,iCell) +
                std::max(0.0,w(k,iCell))*tracerCur(k  ,iCell);
            highOrderFlx(k,iCell) = highOrderFlx(k,iCell) -
                                     lowOrderFlx(k,iCell);
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        // Upwind fluxes are accumulated in workTend
        // flxIn contains the total remaining high order flux into iCell
        //          it is positive.
        // flxOut contains the total remaining high order flux out of iCell
        //           it is negative
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            workTend(k,iCell) = lowOrderFlx(k+1,iCell)
                              - lowOrderFlx(k  ,iCell);
            flxIn (k, iCell) = std::max(0.0, highOrderFlx(k+1, iCell))
                             - std::min(0.0, highOrderFlx(k  , iCell));
            flxOut(k, iCell) = std::min(0.0, highOrderFlx(k+1, iCell))
                             - std::max(0.0, highOrderFlx(k  , iCell));
        }

    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    /*
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nCells}) ,
    YAKL_LAMBDA(int iCell)
    {
        int kmin = minLevelCell(iCell);
        int kmax = maxLevelCell(iCell);
        // Next-to-top cell in column is second-order
        highOrderFlx(kmin,iCell) = 0.0;
        if (kmax > kmin)
        {
            int k = kmin+1;
            double verticalWeightK   = hProv(k-1,iCell) /
                               (hProv(k,iCell) + hProv(k-1, iCell));
            double verticalWeightKm1 = hProv(k,iCell) /
                               (hProv(k,iCell) + hProv(k-1, iCell));
            highOrderFlx(k,iCell) = w(k,iCell)*
                               (verticalWeightK  *tracerCur(k,iCell) +
                                verticalWeightKm1*tracerCur(k-1,iCell));
        }
        // Deepest vertical cell in column is second order
        int k = max(kmin+1,kmax);
        double verticalWeightK   = hProv(k-1,iCell) /
                             (hProv(k  ,iCell) + hProv(k-1,iCell));
        double verticalWeightKm1 = hProv(k  ,iCell) /
                             (hProv(k  ,iCell) + hProv(k-1,iCell));
        highOrderFlx(k,iCell) = w(k,iCell)*
                             (verticalWeightK  *tracerCur(k  ,iCell) +
                              verticalWeightKm1*tracerCur(k-1,iCell));
        highOrderFlx(k+1,iCell) = 0.0;

        // Compute low order (upwind) flux and remove from high order
        lowOrderFlx(kmin,iCell) = 0.0;
        for ( k = kmin+1; k <= kmax; ++k )
        {
            lowOrderFlx(k,iCell) =
                std::min(0.0,w(k,iCell))*tracerCur(k-1,iCell) +
                std::max(0.0,w(k,iCell))*tracerCur(k  ,iCell);
            highOrderFlx(k,iCell) = highOrderFlx(k,iCell) -
                                     lowOrderFlx(k,iCell);
        }
        lowOrderFlx(kmax+1,iCell) = 0.0;

        // Upwind fluxes are accumulated in workTend
        // flxIn contains the total remaining high order flux into iCell
        //          it is positive.
        // flxOut contains the total remaining high order flux out of iCell
        //           it is negative
        for ( k = kmin; k <= kmax; ++k )
        {
            workTend(k,iCell) = lowOrderFlx(k+1,iCell)
                              - lowOrderFlx(k  ,iCell);
            flxIn (k, iCell) = std::max(0.0, highOrderFlx(k+1, iCell))
                             - std::min(0.0, highOrderFlx(k  , iCell));
            flxOut(k, iCell) = std::min(0.0, highOrderFlx(k+1, iCell))
                             - std::max(0.0, highOrderFlx(k  , iCell));
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
    */
}

extern "C"
void ocn_advect_mono_update_tend(double dt, int iTracer, int nCells, int mxlvl)
{
    YAKL_LOCAL_NS(tracer_mono,highOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,lowOrderFlx);
    YAKL_LOCAL_NS(tracer_mono,flxIn);
    YAKL_LOCAL_NS(tracer_mono,flxOut);
    YAKL_LOCAL_NS(tracer_mono,workTend);
    YAKL_LOCAL_NS(tracer_mono,tracerCur);
    YAKL_LOCAL_NS(tracer_mono,tracerMin);
    YAKL_LOCAL_NS(tracer_mono,tracerMax);
    YAKL_LOCAL_NS(tracer_mono,hProv);
    YAKL_LOCAL_NS(tracer_mono,hNewInv);
    YAKL_LOCAL_NS(tracer_mono,tend);
    YAKL_LOCAL_NS(mesh,minLevelCell);
    YAKL_LOCAL_NS(mesh,maxLevelCell);
    YAKL_SCOPE(eps, tracer_mono::eps);

    // Build the scale factors to limit flux for FCT
    // Computed using the bounds that were computed previously,
    // and the bounds on the newly updated value
    // Factors are placed in the flxIn and flxOut arrays

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,mxlvl}),
    YAKL_LAMBDA(int iCell, int k)
    {
        // workTend on the RHS is the upwind tendency
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            double tracerMinNew = (tracerCur(k,iCell)*hProv(k,iCell)
                         + dt*(workTend(k,iCell)+flxOut(k,iCell)))
                         * hNewInv(k,iCell);
            double tracerMaxNew = (tracerCur(k,iCell)*hProv(k,iCell)
                         + dt*(workTend(k,iCell)+flxIn(k,iCell)))
                         * hNewInv(k,iCell);
            double tracerUpwindNew = (tracerCur(k,iCell)*hProv(k,iCell)
                            + dt*workTend(k,iCell)) * hNewInv(k,iCell);

            double scaleFactor = (tracerMax(k,iCell)-tracerUpwindNew)/
                          (tracerMaxNew-tracerUpwindNew+eps);
            flxIn (k,iCell) = std::min(1.0, std::max(0.0, scaleFactor));

            scaleFactor = (tracerUpwindNew-tracerMin(k,iCell))/
                          (tracerUpwindNew-tracerMinNew+eps);
            flxOut(k,iCell) = std::min(1.0, std::max(0.0, scaleFactor));
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    // Accumulate the scaled high order vertical tendencies
    // and the upwind tendencies
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,mxlvl}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        // rescale the high order vertical flux
        if ( (k > minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            double flux =  highOrderFlx(k,iCell);
            highOrderFlx(k,iCell) = std::max(0.0,flux)*
                                    std::min(flxOut(k  ,iCell),
                                        flxIn (k-1,iCell))
                                  + std::min(0.0,flux)*
                                    std::min(flxOut(k-1,iCell),
                                        flxIn (k  ,iCell));
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,mxlvl}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            // workTend on the RHS is the upwind tendency
            // workTend on the LHS is the total vertical advection tendency
            workTend(k, iCell) = workTend(k, iCell)
                               + (highOrderFlx(k+1,iCell)
                                - highOrderFlx(k  ,iCell));
            tend(k,iCell,iTracer) = tend(k,iCell,iTracer) + workTend(k,iCell);
        }
    }, yakl::LaunchConfig<>(), tracer_mono::stream);
}

extern "C"
void ocn_advect_mono_init(double f_coef3rdOrder, double f_eps)
{
    using namespace tracer_mono;
    using namespace mesh;
    
    yakl::streamCreate(&stream);
    
    coef3rdOrder = f_coef3rdOrder;
    eps = f_eps;
    hProvInv = yakl_create_real("hProvInv", nVertLevels, nCellsAll);
    flxIn    = yakl_create_real("flxIn", nVertLevels, nCellsAll+1);
    flxOut   = yakl_create_real("flxOut", nVertLevels, nCellsAll+1);
    highOrderFlx = yakl_create_real("highOrderFlx", nVertLevels+1,std::max(nCellsAll,nEdgesAll)+1);
    hNewInv  = yakl_create_real("hNewInv", nVertLevels, nCellsAll);
    hProv    = yakl_create_real("hProv", nVertLevels, nCellsAll);
    lowOrderFlx = yakl_create_real("lowOrderFlx", nVertLevels+1,std::max(nCellsAll,nEdgesAll)+1);
    tracerCur = yakl_create_real("tracerCur", nVertLevels, nCellsAll+1);
    tracerMax = yakl_create_real("tracerMax", nVertLevels, nCellsAll);
    tracerMin = yakl_create_real("tracerMin", nVertLevels, nCellsAll);
    workTend = yakl_create_real("workTend", nVertLevels, nCellsAll+1);
    wgtTmp = yakl_create_real("wgtTmp", nVertLevels, nEdgesAll);
    sgnTmp = yakl_create_real("sgnTmp", nVertLevels, nEdgesAll);
    normalThicknessFlux = yakl_create_real("normalThicknessFlux", nVertLevels, nEdgesAll+1);
    w = yakl_create_real("w", nVertLevels+1, nCellsAll+1);
    layerThickness = yakl_create_real("layerThickness", nVertLevels, nCellsAll);
}
