#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"
#include "mpas_ocn_mesh_c.hxx"
#include "mpas_ocn_diagnostics_yakl_c.hxx"

#define GET_DPTR(v) (static_cast<double *>((vmix::c_##v).ptr))
#define GET_DDPTR(v) (static_cast<double *>((diag_solve::c_##v).ptr))

namespace vmix {
d_double_1d_t   * bottomDrag
                ;
                
d_double_2d_t   * vertDiffTopOfCell,
                * vertViscTopOfCell,
                * vertViscTopOfEdge,
                * tracerGroupSurfaceFlux = nullptr
                ;

d_double_3d_t   * tracersGroup = nullptr,
                * vertNonLocalFlux
                ;
                
extern "C" ocn_yakl_type c_vertViscTopOfCell;
extern "C" ocn_yakl_type c_vertDiffTopOfCell;
extern "C" ocn_yakl_type c_vertViscTopOfEdge;
extern "C" ocn_yakl_type c_vertNonLocalFlux;
extern "C" ocn_yakl_type c_tracerGroupSurfaceFlux;
extern "C" ocn_yakl_type c_bottomDrag;
extern "C" ocn_yakl_type c_tracersGroup;

}

extern "C"
void ocn_vmix_yakl_init() {
    using namespace vmix;
    YAKL_SCOPE(nEdgesAll, mesh::nEdgesAll);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    
    vertDiffTopOfCell = yakl_create_real("vertDiffTopOfCell",c_vertDiffTopOfCell.shape[0],c_vertDiffTopOfCell.shape[1]);
    vertViscTopOfCell = yakl_create_real("vertViscTopOfCell",c_vertViscTopOfCell.shape[0],c_vertViscTopOfCell.shape[1]);
    vertViscTopOfEdge = yakl_create_real("vertViscTopOfEdge",c_vertViscTopOfEdge.shape[0],c_vertViscTopOfEdge.shape[1]);
    vertNonLocalFlux = yakl_create_real("vertNonLocalFlux",c_vertNonLocalFlux.shape[0],c_vertNonLocalFlux.shape[1],c_vertNonLocalFlux.shape[1]);
    bottomDrag = yakl_create_real("bottomDrag",c_bottomDrag.shape[0]);
}

extern "C"
void ocn_vmix_upd_device() {
    yakl_update_device(vmix::bottomDrag, GET_DPTR(bottomDrag));
    yakl_update_device(vmix::vertDiffTopOfCell, GET_DPTR(vertDiffTopOfCell));
    yakl_update_device(vmix::vertViscTopOfEdge, GET_DPTR(vertViscTopOfEdge));
    yakl_update_device(vmix::vertViscTopOfCell, GET_DPTR(vertViscTopOfCell));
    yakl_update_device(vmix::vertNonLocalFlux, GET_DPTR(vertNonLocalFlux));
    yakl_update_device(diag_solve::normalVelocity, GET_DDPTR(normalVelocity));
    yakl_update_device(diag_solve::layerThickEdgeMean, GET_DDPTR(layerThickEdgeMean));
    yakl_update_device(diag_solve::kineticEnergyCell, GET_DDPTR(kineticEnergyCell));
}

extern "C"
void ocn_vmix_wrap_sfcflux() {
    using namespace vmix;

    if ( nullptr == tracerGroupSurfaceFlux ) {
        tracerGroupSurfaceFlux = yakl_wrap_array("tracerGroupSurfaceFlux", GET_DPTR(tracerGroupSurfaceFlux),
            c_tracerGroupSurfaceFlux.shape[0], c_tracerGroupSurfaceFlux.shape[1]);
    }
    else if ( tracerGroupSurfaceFlux->extent(0) != c_tracerGroupSurfaceFlux.shape[0] ) {
        yakl_delete(tracerGroupSurfaceFlux);
        tracerGroupSurfaceFlux = yakl_wrap_array("tracerGroupSurfaceFlux", GET_DPTR(tracerGroupSurfaceFlux),
            c_tracerGroupSurfaceFlux.shape[0], c_tracerGroupSurfaceFlux.shape[1]);
    }
    else {
        yakl_update_device(vmix::tracerGroupSurfaceFlux, GET_DPTR(tracerGroupSurfaceFlux));
    }

    if ( nullptr == tracersGroup ) {
        tracersGroup = yakl_wrap_array("tracersGroup", GET_DPTR(tracersGroup),
            c_tracersGroup.shape[0], c_tracersGroup.shape[1], c_tracersGroup.shape[2]);
    }
    else if ( tracersGroup->extent(0) != c_tracersGroup.shape[0] ) {
        yakl_delete(tracersGroup);
        tracersGroup = yakl_wrap_array("tracersGroup", GET_DPTR(tracersGroup),
            c_tracersGroup.shape[0], c_tracersGroup.shape[1], c_tracersGroup.shape[2]);
    }
    else {
        yakl_update_device(vmix::tracersGroup, GET_DPTR(tracersGroup));
    }
}


extern "C"
void ocn_tracer_vmix_tend_implicit(double dt, bool config_cvmix_kpp_nonlocal_with_implicit_mix) {
    YAKL_LOCAL_NS(vmix, vertDiffTopOfCell);
    YAKL_LOCAL_NS(vmix, vertNonLocalFlux);
    YAKL_LOCAL_NS(vmix, tracerGroupSurfaceFlux);
    YAKL_LOCAL_NS(vmix, tracersGroup);

    YAKL_LOCAL_NS(diag_solve, layerThickness);
    YAKL_LOCAL_NS(diag_solve, bTemp);
    YAKL_LOCAL_NS(diag_solve, rTemp);
    YAKL_LOCAL_NS(diag_solve, cTemp);

    YAKL_LOCAL_NS(mesh, minLevelCell);
    YAKL_LOCAL_NS(mesh, maxLevelCell);

    YAKL_SCOPE(nCellsOwned, mesh::nCellsOwned);

    int num_tracers = tracersGroup.extent(0);
    
    if ( config_cvmix_kpp_nonlocal_with_implicit_mix ) {
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCellsOwned},{1,num_tracers}) ,
        YAKL_LAMBDA(int iCell,int iTracer) {
            int Nsurf = minLevelCell(iCell);
            int N = maxLevelCell(iCell);
            rTemp(iTracer,Nsurf) = tracersGroup(iTracer,Nsurf,iCell) + dt*tracerGroupSurfaceFlux(iTracer,iCell)
                    * (-vertNonLocalFlux(1, Nsurf+1,iCell) )/ layerThickness(Nsurf,iCell);
        }, yakl::DefaultLaunchConfig());
    }
    else
    {
    }
}


extern "C"
void ocn_vmix_vel_tend_impl(double dt, double implicitBottomDragCoef) {

    YAKL_LOCAL_NS(vmix, vertViscTopOfEdge);

    YAKL_LOCAL_NS(diag_solve, bTemp);
    YAKL_LOCAL_NS(diag_solve, rTemp);
    YAKL_LOCAL_NS(diag_solve, cTemp);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, layerThicknessEdge);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeFlux);
    YAKL_LOCAL_NS(diag_solve, kineticEnergyCell);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeMean);

    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k < minLevelEdgeBot(iEdge)) || (k > maxLevelEdgeTop(iEdge)) )
            normalVelocity(k,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned),
    YAKL_LAMBDA(int iEdge)
    {
        int N = maxLevelEdgeTop(iEdge);
        int Nsurf = minLevelEdgeBot(iEdge);
        if (N > 0) {
            if (N == Nsurf) {
                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                normalVelocity(N,iEdge) = normalVelocity(N,iEdge)
                    / (1.0 + dt*implicitBottomDragCoef
                       * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) );
            }
            else {
                // tridiagonal matrix algorithm
                cTemp(Nsurf,iEdge) = -2.0 * dt * vertViscTopOfEdge(Nsurf+1,iEdge)
                              / (layerThickEdgeMean(Nsurf,iEdge) + layerThickEdgeMean(Nsurf+1,iEdge))
                              / layerThickEdgeMean(Nsurf,iEdge);
                bTemp(Nsurf,iEdge) = 1.0 - cTemp(Nsurf,iEdge);
                rTemp(Nsurf,iEdge) = normalVelocity(Nsurf,iEdge);
            }
        }
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        int N = maxLevelEdgeTop(iEdge);
        int Nsurf = minLevelEdgeBot(iEdge);
        if ( (k > Nsurf) && (k < N) && (Nsurf < N) ) {
            double A        = -2.0*dt*vertViscTopOfEdge(k,iEdge)
                 / (layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge))
                 / layerThickEdgeMean(k,iEdge);
            double m        = A/bTemp(k-1, iEdge);
            cTemp(k,iEdge)     = -2.0*dt*vertViscTopOfEdge(k+1,iEdge)
                 / (layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k+1,iEdge))
                 / layerThickEdgeMean(k,iEdge);
            bTemp(k,iEdge) = 1.0 - A - cTemp(k,iEdge) - m*cTemp(k-1,iEdge);
            rTemp(k,iEdge) = normalVelocity(k,iEdge) - m*rTemp(k-1,iEdge);
        }
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned),
    YAKL_LAMBDA(int iEdge)
    {
        int N = maxLevelEdgeTop(iEdge);
        int Nsurf = minLevelEdgeBot(iEdge);
        if ( (N > 0) && (Nsurf < N) ) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);

            double A = -2.0*dt*vertViscTopOfEdge(N,iEdge)
                 / (layerThickEdgeMean(N-1,iEdge) + layerThickEdgeMean(N,iEdge))
                 / layerThickEdgeMean(N,iEdge);
            double m = A/bTemp(N-1, iEdge);

           // x(N) = rTemp(N) / bTemp(N)
           // Apply bottom drag boundary condition on the viscous term
           // using sqrt(2.0*kineticEnergyEdge(k,iEdge))
           normalVelocity(N,iEdge) = (normalVelocity(N,iEdge) - m*rTemp(N-1, iEdge))
               / (1.0 - A + dt*implicitBottomDragCoef
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge)
               - m*cTemp(N-1,iEdge));

           // second pass: back substitution
           for ( int k = N-1; k >= Nsurf; --k ) {
              normalVelocity(k,iEdge) = (rTemp(k,iEdge) - cTemp(k,iEdge)*normalVelocity(k+1,iEdge)) / bTemp(k,iEdge);
           }
        }
    }, yakl::DefaultLaunchConfig());

}


extern "C"
void ocn_vmix_vel_tend_impl_spatial(double dt, double implicitBottomDragCoef) {


    YAKL_LOCAL_NS(vmix, vertViscTopOfEdge);
    YAKL_LOCAL_NS(vmix, bottomDrag);

    YAKL_LOCAL_NS(diag_solve, bTemp);
    YAKL_LOCAL_NS(diag_solve, rTemp);
    YAKL_LOCAL_NS(diag_solve, cTemp);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, layerThicknessEdge);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeFlux);
    YAKL_LOCAL_NS(diag_solve, kineticEnergyCell);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeMean);

    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k < minLevelEdgeBot(iEdge)) || (k > maxLevelEdgeTop(iEdge)) )
            normalVelocity(k,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned),
    YAKL_LAMBDA(int iEdge)
    {
        int N = maxLevelEdgeTop(iEdge);
        int Nsurf = minLevelEdgeBot(iEdge);
        if (N > 0) {
            if (N == Nsurf) {
                int cell1 = cellsOnEdge(1,iEdge);
                int cell2 = cellsOnEdge(2,iEdge);
                //  average cell-based implicit bottom drag to edges
                double implicitCD = 0.5 * (bottomDrag(cell1) + bottomDrag(cell2));
                normalVelocity(N,iEdge) = normalVelocity(N,iEdge)
                    / (1.0 + dt*implicitCD
                       * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge) );
            }
            else {
                // tridiagonal matrix algorithm
                cTemp(Nsurf,iEdge) = -2.0 * dt * vertViscTopOfEdge(Nsurf+1,iEdge)
                              / (layerThickEdgeMean(Nsurf,iEdge) + layerThickEdgeMean(Nsurf+1,iEdge))
                              / layerThickEdgeMean(Nsurf,iEdge);
                bTemp(Nsurf,iEdge) = 1.0 - cTemp(Nsurf,iEdge);
                rTemp(Nsurf,iEdge) = normalVelocity(Nsurf,iEdge);
            }
        }
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        int N = maxLevelEdgeTop(iEdge);
        int Nsurf = minLevelEdgeBot(iEdge);
        if ( (k > Nsurf) && (k < N) && (Nsurf < N) ) {
            double A        = -2.0*dt*vertViscTopOfEdge(k,iEdge)
                 / (layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge))
                 / layerThickEdgeMean(k,iEdge);
            double m        = A/bTemp(k-1, iEdge);
            cTemp(k,iEdge)     = -2.0*dt*vertViscTopOfEdge(k+1,iEdge)
                 / (layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k+1,iEdge))
                 / layerThickEdgeMean(k,iEdge);
            bTemp(k,iEdge) = 1.0 - A - cTemp(k,iEdge) - m*cTemp(k-1,iEdge);
            rTemp(k,iEdge) = normalVelocity(k,iEdge) - m*rTemp(k-1,iEdge);
        }
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned),
    YAKL_LAMBDA(int iEdge)
    {
        int N = maxLevelEdgeTop(iEdge);
        int Nsurf = minLevelEdgeBot(iEdge);
        if ( (N > 0) && (Nsurf < N) ) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);

            double implicitCD = 0.5*(bottomDrag(cell1) + bottomDrag(cell2));

            double A = -2.0*dt*vertViscTopOfEdge(N,iEdge)
                 / (layerThickEdgeMean(N-1,iEdge) + layerThickEdgeMean(N,iEdge))
                 / layerThickEdgeMean(N,iEdge);
            double m = A/bTemp(N-1, iEdge);

           // x(N) = rTemp(N) / bTemp(N)
           // Apply bottom drag boundary condition on the viscous term
           // using sqrt(2.0*kineticEnergyEdge(k,iEdge))
           normalVelocity(N,iEdge) = (normalVelocity(N,iEdge) - m*rTemp(N-1, iEdge))
               / (1.0 - A + dt*implicitCD
               * sqrt(kineticEnergyCell(N,cell1) + kineticEnergyCell(N,cell2)) / layerThickEdgeMean(N,iEdge)
               - m*cTemp(N-1,iEdge));

           // second pass: back substitution
           for ( int k = N-1; k >= Nsurf; --k ) {
              normalVelocity(k,iEdge) = (rTemp(k,iEdge) - cTemp(k,iEdge)*normalVelocity(k+1,iEdge)) / bTemp(k,iEdge);
           }
        }
    }, yakl::DefaultLaunchConfig());

}


extern "C"
void ocn_vmix_avg() {
    YAKL_LOCAL_NS(vmix, vertViscTopOfEdge);
    YAKL_LOCAL_NS(vmix, vertViscTopOfCell);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k < minLevelEdgeBot(iEdge)) || (k > maxLevelEdgeTop(iEdge)) )
            vertViscTopOfEdge(k,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            vertViscTopOfEdge(k,iEdge) = 0.5*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1));
        }
    }, yakl::DefaultLaunchConfig());
}

extern "C"
void ocn_vmix_gotm_avg() {
    YAKL_LOCAL_NS(vmix, vertViscTopOfEdge);
    YAKL_LOCAL_NS(vmix, vertViscTopOfCell);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k < minLevelEdgeBot(iEdge)) || (k > maxLevelEdgeTop(iEdge)) )
            vertViscTopOfEdge(k,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            vertViscTopOfEdge(k,iEdge) = 0.5*(vertViscTopOfCell(k,cell2)+vertViscTopOfCell(k,cell1));
        }
    }, yakl::DefaultLaunchConfig());
}
