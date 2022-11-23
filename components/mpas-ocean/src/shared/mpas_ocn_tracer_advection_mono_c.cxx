#include "mpas_ocn_diagnostics_variables.hxx"
#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_mesh_c.hxx"
#include "mpas_ocn_yakl_c.hxx"
#include <cmath>
/***
 mpas_ocn_mesh_c.hxx:
   advCoefs, advCoefs3rd,nAdvCellsForEdge,advCellsForEdge
**/

#define DMIN(a,b) std::min<double>(a,b)
#define DMAX(a,b) std::max<double>(a,b)
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

yakl::Stream     stream;

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

extern "C" ocn_yakl_type c_transTend;

extern "C"
void ocn_advect_mono_copyin(double * h_normalThicknessFlux, 
                                double * h_w, double * h_layerThickness)
{
    using namespace tracer_mono;
    
    //std::cerr << " tend size = " << tendSize[0] << " " << tendSize[1] << " " << tendSize[2] << std::endl;
    tend = yakl_wrap_array("tend", static_cast<double *>(c_transTend.ptr),
            c_transTend.shape[0], c_transTend.shape[1], c_transTend.shape[2]);
    yakl_update_device(normalThicknessFlux, h_normalThicknessFlux);
    yakl_update_device(w, h_w);
    yakl_update_device(layerThickness, h_layerThickness);
}

extern "C"
void ocn_advect_mono_end()
{
    //yakl_update_host(tracer_mono::tend, h_tend, tracer_mono::stream);
    yakl_update_host(tracer_mono::tend, static_cast<double *>(c_transTend.ptr));
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
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());
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

    YAKL_LOCAL_NS(mesh, areaCell);
    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);
    YAKL_LOCAL_NS(mesh, nEdgesOnCell);
    YAKL_LOCAL_NS(mesh, edgesOnCell);
    YAKL_LOCAL_NS(mesh, dvEdge);
    YAKL_LOCAL_NS(mesh, edgeSignOnCell);
    YAKL_LOCAL_NS(mesh, invAreaCell);

    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nCells, mesh::nCellsAll);

//#if 1
#ifdef MPAS_DEBUG

    // original Fortran loop structure
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCells) ,
    YAKL_LAMBDA(int iCell) {
         double invAreaCell1 = dt/areaCell(iCell);
         int kmin = minLevelCell(iCell);
         int kmax = maxLevelCell(iCell);
         for ( int k = kmin; k <= kmax; ++k ) {
            hProv(k, iCell) = layerThickness(k, iCell);
         }

         for ( int i = 1; i <= nEdgesOnCell(iCell); ++i ) {
            int iEdge = edgesOnCell(i,iCell);
            double signedFactor = invAreaCell1*dvEdge(iEdge)*
                           edgeSignOnCell(i,iCell);
            // Provisional layer thickness is after horizontal
            // thickness flux only
            for ( int k = kmin; k <= kmax; ++k ) {
               hProv(k,iCell) = hProv(k,iCell)
                              + signedFactor*normalThicknessFlux(k,iEdge);
            }
         }

         // New layer thickness is after horizontal and vertical
         // thickness flux
         for ( int k = kmin; k <= kmax; ++k ) {
            hProvInv(k,iCell) = 1.0/ hProv(k,iCell);
            hNewInv (k,iCell) = 1.0/(hProv(k,iCell) -
                                dt*w(k,iCell) + dt*w(k+1, iCell));
         }
    }, yakl::LaunchConfig<>());

#else

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        // Using invAreaCell gives diffs
        //double invAreaCell1 = dt * invAreaCell(iCell);
         double invAreaCell1 = dt/areaCell(iCell);
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
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());
#endif

//    yakl_update_host(tracer_mono::hProvInv, h_hProvInv);
//    yakl_update_host(tracer_mono::hNewInv, h_hNewInv);
//    yakl_update_host(tracer_mono::hProv, h_hProv);
//    yakl::fence();

}

extern "C"
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
    YAKL_LOCAL_NS(mesh, advMaskHighOrder);
    YAKL_LOCAL_NS(mesh, advCellsForEdge);
    YAKL_LOCAL_NS(mesh, nAdvCellsForEdge);
    YAKL_LOCAL_NS(mesh, advCoefs);
    YAKL_LOCAL_NS(mesh, advCoefs3rd);
    YAKL_LOCAL_NS(mesh, dvEdge);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);

    //yakl_update_device(tracer_mono::tracerCur, h_tracerCur, tracer_mono::stream);
    yakl_update_device(tracer_mono::tracerCur, h_tracerCur);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        if ( (k >= minLevelCell(iCell)) &&  (k <= maxLevelCell(iCell)) )
        {
            tracerMin(k,iCell) = tracerCur(k,iCell);
            tracerMax(k,iCell) = tracerCur(k,iCell);
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        for (int i = 1; i <= nEdgesOnCell(iCell); ++i )
        {
            int cell2 = cellsOnCell(i,iCell);
            int kmin  = DMAX(minLevelCell(iCell),
                        minLevelCell(cell2));
            int kmax  = DMIN(maxLevelCell(iCell), maxLevelCell(cell2));
            if ( (k >= kmin) && (k <= kmax) )
            {
              tracerMax(k,iCell) = DMAX(tracerMax(k,iCell), tracerCur(k,cell2));
              tracerMin(k,iCell) = DMIN(tracerMin(k,iCell), tracerCur(k,cell2));
            }
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        // compute some common intermediate factors
        wgtTmp(k,iEdge) = normalThicknessFlux(k,iEdge)*advMaskHighOrder(k,iEdge);
        sgnTmp(k,iEdge) = SIGN(1.0, normalThicknessFlux(k,iEdge));
        highOrderFlx(k,iEdge) = 0.0;
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

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
               //highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + 
                 //   normalThicknessFlux(k,iEdge);
               highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge) + tracerCur(k,iCell)*
                          wgtTmp(k,iEdge)*(coef1 + coef3*sgnTmp(k,iEdge));
            }
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

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
            double tracerWeight = (1.0 - advMaskHighOrder(k,iEdge))
                         * (dvEdge(iEdge) * 0.5)
                         * normalThicknessFlux(k, iEdge);

            lowOrderFlx(k,iEdge) = dvEdge(iEdge) * 
               (DMAX(0.0,normalThicknessFlux(k,iEdge))*tracerCur(k,cell1)
              + DMIN(0.0,normalThicknessFlux(k,iEdge))*tracerCur(k,cell2));

            highOrderFlx(k, iEdge) = highOrderFlx(k, iEdge)
                                   + tracerWeight * (tracerCur(k, cell1)
                                                   + tracerCur(k, cell2));

            highOrderFlx(k,iEdge) = highOrderFlx(k,iEdge)
                                  -  lowOrderFlx(k,iEdge);
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

//    yakl_update_host(tracer_mono::highOrderFlx, h_highOrderFlx);
//    yakl_update_host(tracer_mono::lowOrderFlx, h_lowOrderFlx);
//    yakl_update_host(tracer_mono::tracerMin, h_tracerMin);
//    yakl_update_host(tracer_mono::tracerMax, h_tracerMax);
//    yakl::fence();
}

extern "C"
void ocn_advect_mono_set_flux_inout(int mpas_myrank, int initNCells, int nCellsHalo, int nEdges, double dt,
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
    YAKL_SCOPE(nCells, mesh::nCellsAll);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells+1},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        workTend(k, iCell) = 0.0;
        flxIn   (k, iCell) = 0.0;
        flxOut  (k, iCell) = 0.0;
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());


//#if 1
#ifdef MPAS_DEBUG
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCellsHalo) ,
    YAKL_LAMBDA(int iCell)
    {
        //double invAreaCell1 = invAreaCell(iCell);
        double invAreaCell1 = 1.0 / areaCell(iCell);

        for ( int i = 1; i <= nEdgesOnCell(iCell); ++i ) {
            int iEdge = edgesOnCell(i, iCell);
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            double signedFactor = edgeSignOnCell(i, iCell) * invAreaCell1;

            for ( int k = minLevelEdgeBot(iEdge); k <= maxLevelEdgeTop(iEdge); ++k ) {
              workTend(k,iCell) = workTend(k,iCell)
                                + signedFactor*lowOrderFlx(k,iEdge);

              // Accumulate remaining high order fluxes
              flxOut(k,iCell) = flxOut(k,iCell) + DMIN(0.0,
                                signedFactor*highOrderFlx(k,iEdge));
              flxIn (k,iCell) = flxIn (k,iCell) + DMAX(0.0,
                                signedFactor*highOrderFlx(k,iEdge));
            }
        }
        
        for ( int k = minLevelCell(iCell); k <= maxLevelCell(iCell); ++k ) {
            // Here workTend is the upwind tendency
            double tracerUpwindNew = (tracerCur(k,iCell)*layerThickness(k,iCell)
                             + dt*workTend(k,iCell)) * hProvInv(k,iCell);
            double tracerMinNew = tracerUpwindNew
                         + dt*flxOut(k,iCell)*hProvInv(k,iCell);
            double tracerMaxNew = tracerUpwindNew
                         + dt*flxIn (k,iCell)*hProvInv(k,iCell);

            double scaleFactor = (tracerMax(k,iCell) - tracerUpwindNew)/
                          (tracerMaxNew - tracerUpwindNew + eps);
            flxIn (k,iCell) = DMIN(1.0, DMAX(0.0, scaleFactor));

            scaleFactor = (tracerUpwindNew - tracerMin(k,iCell))/
                          (tracerUpwindNew - tracerMinNew + eps);
            flxOut(k,iCell) = DMIN(1.0, DMAX(0.0, scaleFactor));
        }        
    }, yakl::LaunchConfig<>());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdges) ,
    YAKL_LAMBDA(int iEdge) {
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        for ( int k = minLevelEdgeBot(iEdge); k <= maxLevelEdgeTop(iEdge); ++k ) {
              highOrderFlx(k,iEdge) = DMAX(0.0,highOrderFlx(k,iEdge))*
                                      DMIN(flxOut(k,cell1), flxIn (k,cell2))
                                    + DMIN(0.0,highOrderFlx(k,iEdge))*
                                      DMIN(flxIn (k,cell1), flxOut(k,cell2));        
        }
    }, yakl::LaunchConfig<>());

#else

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
            flxIn (k,iCell) = DMIN(1.0, DMAX(0.0, scaleFactor));

            scaleFactor = (tracerUpwindNew - tracerMin(k,iCell))/
                          (tracerUpwindNew - tracerMinNew + eps);
            flxOut(k,iCell) = DMIN(1.0, DMAX(0.0, scaleFactor));
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
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());
#endif
    
//    yakl_update_host(tracer_mono::highOrderFlx, h_highOrderFlx);
//    yakl_update_host(tracer_mono::flxIn, h_flxIn);
//    yakl_update_host(tracer_mono::flxOut, h_flxOut);
//    yakl_update_host(tracer_mono::workTend, h_workTend);
//    yakl::fence();
}


extern "C"
void ocn_advect_mono_compute_tend(int iTracer, int nCells, double dt 
                                )
                                  //double * h_tend, double * h_tracerCur)
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
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());
        
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
    YAKL_LAMBDA(int iCell, int k)
    {
        if ( (k >= minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            // workTend on the RHS is the total horizontal advection tendency
            // tracerCur on LHS is the  provisional tracer after horizontal fluxes only.
            //tracerCur(k,iCell) = (layerThickness(k,iCell)
              //                    + dt*workTend(k,iCell))*hProvInv(k,iCell);
            tracerCur(k,iCell) = (tracerCur(k,iCell)*layerThickness(k,iCell)
                                  + dt*workTend(k,iCell))*hProvInv(k,iCell);
            tend(k,iCell,iTracer) += workTend(k,iCell);
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());


    //yakl_update_host(tracer_mono::tracerCur, h_tracerCur);
    //yakl_update_host(&tslice, h_tend);
    //yakl::fence();
}

extern "C"
void ocn_advect_mono_vert_flux(int rank, int nCells, int vertOrder,
                                //double * h_flxIn, double * h_flxOut,
                                double * h_lowOrderFlx, double * h_highOrderFlx,
                                double * h_tracerMin, double * h_tracerMax)
                                //double * h_workTend)
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
    YAKL_SCOPE(nCellsAll, mesh::nCellsAll);
    
    
    //std::cerr << " vertOrder = " << vertOrder << std::endl;
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCells) ,
    YAKL_LAMBDA(int iCell)
    {
        int kmin = minLevelCell(iCell);
        int kmax = maxLevelCell(iCell);

        // take care of top cell
        tracerMax(kmin,iCell) = DMAX(tracerCur(1,iCell),
                                   tracerCur(2,iCell));
        tracerMin(kmin,iCell) = DMIN(tracerCur(1,iCell),
                                   tracerCur(2,iCell));
        for ( int k=kmin+1; k <= kmax-1; ++k )
        {
            tracerMax(k,iCell) = DMAX(DMAX(tracerCur(k-1,iCell),
                                     tracerCur(k  ,iCell)),
                                     tracerCur(k+1,iCell));
            tracerMin(k,iCell) = DMIN(DMIN(tracerCur(k-1,iCell),
                                     tracerCur(k  ,iCell)),
                                     tracerCur(k+1,iCell));
        }
        // finish with bottom cell
        tracerMax(kmax,iCell) = DMAX(tracerCur(kmax  ,iCell),
                                      tracerCur(kmax-1,iCell));
        tracerMin(kmax,iCell) = DMIN(tracerCur(kmax  ,iCell),
                                      tracerCur(kmax-1,iCell));
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

    double inv12 = 1.0 / 12.0;

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsAll},{1,nVertLevels+1}) ,
    YAKL_LAMBDA(int iCell, int k) {
         highOrderFlx(k,iCell) = 0.0;
    }, yakl::LaunchConfig<>());


    switch ( vertOrder )
    {
    case vertOrder4:
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsAll},{1,nVertLevels}),
        YAKL_LAMBDA(int iCell, int k)
        {
            if ( (k >= minLevelCell(iCell)+2) && (k < maxLevelCell(iCell)) )
                highOrderFlx(k, iCell) = w(k,iCell)*(
                     7.0*(tracerCur(k  ,iCell) + tracerCur(k-1,iCell)) -
                               (tracerCur(k+1,iCell) + tracerCur(k-2,iCell))) * inv12;
        //}, yakl::LaunchConfig<>(), tracer_mono::stream);
        }, yakl::LaunchConfig<>());
        break;
    case vertOrder3:
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsAll},{1,nVertLevels}),
        YAKL_LAMBDA(int iCell, int k)
        {
            if ( (k >= minLevelCell(iCell)+2) && (k < maxLevelCell(iCell)) )
                highOrderFlx(k, iCell) = (w(k,iCell)*
                      (7.0*(tracerCur(k  ,iCell) + tracerCur(k-1,iCell)) -
                                 (tracerCur(k+1,iCell) + tracerCur(k-2,iCell))) -
                        coef3rdOrder*std::abs<double>(w(k,iCell))*
                                ((tracerCur(k+1,iCell) - tracerCur(k-2,iCell)) -
                       3.0*(tracerCur(k  ,iCell) - tracerCur(k-1,iCell)))) / 12.0;
        //}, yakl::LaunchConfig<>(), tracer_mono::stream);
        }, yakl::LaunchConfig<>());
        break;
    case vertOrder2:
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsAll},{3,nVertLevels}) ,
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
        //}, yakl::LaunchConfig<>(), tracer_mono::stream);
        }, yakl::LaunchConfig<>());
        break;
    }

//#if 1
#ifdef MPAS_DEBUG

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCellsAll) ,
    YAKL_LAMBDA(int iCell) {
        int kmin = minLevelCell(iCell);
        int kmax = maxLevelCell(iCell);
        int k;
        // at top, flux is zero (already initialized)
        // at next-to-top (kmin+1), reduce to 2nd order
        //   but avoid case where 0 or 1 active layer (kmax <= kmin)
         if (kmax > kmin) {
            k = kmin+1;
            double verticalWeightK   = hProv(k-1,iCell) /
                               (hProv(k,iCell) + hProv(k-1, iCell));
            double verticalWeightKm1 = hProv(k,iCell) /
                               (hProv(k,iCell) + hProv(k-1, iCell));
            highOrderFlx(k,iCell) = w(k,iCell)*
                               (verticalWeightK  *tracerCur(k,iCell) +
                                verticalWeightKm1*tracerCur(k-1,iCell));
            // Deepest active level also at 2nd order
            k = kmax;
            verticalWeightK   = hProv(k-1,iCell) /
                               (hProv(k  ,iCell) +
                                hProv(k-1,iCell));
            verticalWeightKm1 = hProv(k  ,iCell) /
                               (hProv(k  ,iCell) +
                                hProv(k-1,iCell));
            highOrderFlx(k,iCell) = w(k,iCell)*
                          (verticalWeightK  *tracerCur(k  ,iCell) +
                           verticalWeightKm1*tracerCur(k-1,iCell));
         }
    }, yakl::LaunchConfig<>());


    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCells) ,
    YAKL_LAMBDA(int iCell) {
        int kmin = minLevelCell(iCell);
        int kmax = maxLevelCell(iCell);
        lowOrderFlx(kmin,iCell) = 0.0;
        for ( int k = kmin+1; k <= kmax; ++k ) {
            lowOrderFlx(k,iCell) =
                    DMIN(0.0,w(k,iCell))*tracerCur(k-1,iCell) +
                    DMAX(0.0,w(k,iCell))*tracerCur(k  ,iCell);
            highOrderFlx(k,iCell) = highOrderFlx(k,iCell)
                                    -  lowOrderFlx(k,iCell);
        }
        lowOrderFlx(kmax+1,iCell) = 0.0;

        // Upwind fluxes are accumulated in workTend
        // flxIn  contains total remaining high order flux into iCell
        //          it is positive.
        // flxOut contains total remaining high order flux out of iCell
        //          it is negative
        for ( int k = kmin; k <= kmax; ++k ) {
            workTend(k,iCell) = lowOrderFlx(k+1,iCell)
                                - lowOrderFlx(k  ,iCell);
            flxIn (k,iCell) = DMAX(0.0, highOrderFlx(k+1,iCell))
                              - DMIN(0.0, highOrderFlx(k  ,iCell));
            flxOut(k,iCell) = DMIN(0.0, highOrderFlx(k+1,iCell))
                              - DMAX(0.0, highOrderFlx(k  ,iCell));
        }
    }, yakl::LaunchConfig<>());

#else

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
        int km = std::max(kmin+1,kmax);
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
                DMIN(0.0,w(k,iCell))*tracerCur(k-1,iCell) +
                DMAX(0.0,w(k,iCell))*tracerCur(k  ,iCell);
            highOrderFlx(k,iCell) = highOrderFlx(k,iCell) -
                                     lowOrderFlx(k,iCell);
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());
    
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
            flxIn (k, iCell) = DMAX(0.0, highOrderFlx(k+1, iCell))
                             - DMIN(0.0, highOrderFlx(k  , iCell));
            flxOut(k, iCell) = DMIN(0.0, highOrderFlx(k+1, iCell))
                             - DMAX(0.0, highOrderFlx(k  , iCell));
        }

    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());
#endif // MPAS_DEBUG


    //yakl_update_host(tracer_mono::flxIn, h_flxIn);
    //yakl_update_host(tracer_mono::flxOut, h_flxOut);
    //yakl_update_host(tracer_mono::workTend, h_workTend);
    //yakl_update_host(tracer_mono::tracerMin, h_tracerMin);
    //yakl_update_host(tracer_mono::tracerMax, h_tracerMax);
    //yakl_update_host(tracer_mono::lowOrderFlx, h_lowOrderFlx);
    //yakl_update_host(tracer_mono::highOrderFlx, h_highOrderFlx);
    //yakl::fence();

}

extern "C"
void ocn_advect_mono_update_tend(double dt, int iTracer, int nCells, double * h_workTend)
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
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(eps, tracer_mono::eps);

    // Build the scale factors to limit flux for FCT
    // Computed using the bounds that were computed previously,
    // and the bounds on the newly updated value
    // Factors are placed in the flxIn and flxOut arrays

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}),
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
            flxIn (k,iCell) = DMIN(1.0, DMAX(0.0, scaleFactor));

            scaleFactor = (tracerUpwindNew-tracerMin(k,iCell))/
                          (tracerUpwindNew-tracerMinNew+eps);
            flxOut(k,iCell) = DMIN(1.0, DMAX(0.0, scaleFactor));
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

    // Accumulate the scaled high order vertical tendencies
    // and the upwind tendencies
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell, int k)
    {
        // rescale the high order vertical flux
        if ( (k > minLevelCell(iCell)) && (k <= maxLevelCell(iCell)) )
        {
            double flux =  highOrderFlx(k,iCell);
            highOrderFlx(k,iCell) = DMAX(0.0,flux)*
                                    DMIN(flxOut(k  ,iCell),
                                        flxIn (k-1,iCell))
                                  + DMIN(0.0,flux)*
                                    DMIN(flxOut(k-1,iCell),
                                        flxIn (k  ,iCell));
        }
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
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
    //}, yakl::LaunchConfig<>(), tracer_mono::stream);
    }, yakl::LaunchConfig<>());

    //yakl_update_host(tracer_mono::workTend, h_workTend);
    //yakl::fence();
}

extern "C"
void ocn_advect_mono_update_for_budgets(double * h_lowOrderFlx, double * h_highOrderFlx, double * h_workTend)
{
    yakl_update_host(tracer_mono::lowOrderFlx, h_lowOrderFlx);
    yakl_update_host(tracer_mono::highOrderFlx, h_highOrderFlx);
    yakl_update_host(tracer_mono::workTend, h_workTend);

    yakl::fence();
}

extern "C"
void ocn_advect_mono_update_for_flxchk(double * h_tracerCur, double * h_hProv, 
                                       double * h_hNewInv, double * h_workTend)
{
    yakl_update_host(tracer_mono::tracerCur, h_tracerCur);
    yakl_update_host(tracer_mono::hProv, h_hProv);
    yakl_update_host(tracer_mono::hNewInv, h_hNewInv);
    yakl_update_host(tracer_mono::workTend, h_workTend);

    yakl::fence();   
}

extern "C"
void ocn_advect_mono_update_for_trchk(double * h_tracerCur, double * h_tracerMin, double * h_tracerMax)
{
    yakl_update_host(tracer_mono::tracerCur, h_tracerCur);
    yakl_update_host(tracer_mono::tracerMin, h_tracerMin);
    yakl_update_host(tracer_mono::tracerMax, h_tracerMax);

    yakl::fence();   
}

extern "C"
void ocn_advect_mono_init(double f_coef3rdOrder, double f_eps)
{
    using namespace tracer_mono;
    using namespace mesh;
    
    stream = yakl::create_stream();
    
    coef3rdOrder = f_coef3rdOrder;
    eps = f_eps;
    hProvInv = yakl_create_real("hProvInv", nVertLevels, nCellsAll);
    flxIn    = yakl_create_real("flxIn", nVertLevels, nCellsAll+1);
    flxOut   = yakl_create_real("flxOut", nVertLevels, nCellsAll+1);
    highOrderFlx = yakl_create_real("highOrderFlx", nVertLevels+1,DMAX(nCellsAll,nEdgesAll)+1);
    hNewInv  = yakl_create_real("hNewInv", nVertLevels, nCellsAll);
    hProv    = yakl_create_real("hProv", nVertLevels, nCellsAll);
    lowOrderFlx = yakl_create_real("lowOrderFlx", nVertLevels+1,DMAX(nCellsAll,nEdgesAll)+1);
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
