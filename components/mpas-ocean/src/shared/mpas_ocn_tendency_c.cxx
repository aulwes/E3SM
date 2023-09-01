#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"
#include "mpas_ocn_mesh_c.hxx"
#include "mpas_ocn_diagnostics_yakl_c.hxx"

#define GET_DPTR(v) (static_cast<double *>((diag_solve::c_##v).ptr))
#define GET_TPTR(v) (static_cast<double *>((tendency::c_##v).ptr))

namespace tendency {
d_double_2d_t   * tendVel,
                * w_dudzTopEdge,
                * JacobianTz,
                * JacobianSz
                ;

extern "C" ocn_yakl_type c_tendVel;
}

extern "C"
void ocn_tend_init() {
    using namespace tendency;
    YAKL_SCOPE(nEdgesAll, mesh::nEdgesAll);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    
    tendVel = yakl_create_real("tendVel",c_tendVel.shape[0],c_tendVel.shape[1]);
    w_dudzTopEdge = yakl_create_real("w_dudzTopEdge",nVertLevels+1,nEdgesAll);
    JacobianTz = yakl_create_real("JacobianTz",nVertLevels,nEdgesAll);
    JacobianSz = yakl_create_real("JacobianSz",nVertLevels,nEdgesAll);
}

extern "C"
void ocn_tendency_update_device() {
    yakl_update_device(diag_solve::normalVelocity, GET_DPTR(normalVelocity));
    yakl_update_device(diag_solve::kineticEnergyCell, GET_DPTR(kineticEnergyCell));
    yakl_update_device(diag_solve::layerThickEdgeFlux, GET_DPTR(layerThickEdgeFlux));
    yakl_update_device(diag_solve::vertAleTransportTop, GET_DPTR(vertAleTransportTop));
    yakl_update_device(diag_solve::pressure, GET_DPTR(pressure));
    yakl_update_device(diag_solve::density, GET_DPTR(diag_density));
    yakl_update_device(diag_solve::zMid, GET_DPTR(zMid));
    yakl_update_device(diag_solve::activeTracers, GET_DPTR(activeTracers));

    yakl_update_device(diag_solve::thermExpCoeff, GET_DPTR(thermExpCoeff));
    yakl_update_device(diag_solve::salineContractCoeff, GET_DPTR(salineContractCoeff));

    YAKL_SCOPE(nEdgesAll, mesh::nEdgesAll);
    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_LOCAL_NS(tendency, tendVel);
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesAll},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        tendVel(k,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig());
}

extern "C"
void ocn_tendency_init_fluxes(d_double_2d_t * d_tend_layerThickness,
                              d_double_1d_t * d_surfaceThicknessFlux,
                              d_double_1d_t * d_surfaceThicknessFluxRunoff)
{
    int nCells = d_surfaceThicknessFlux->extent(0);
    int nVertLevels = d_tend_layerThickness->extent(0);
    //std::cerr << " nCells, nVertLevels = " << nCells << " " << nVertLevels << std::endl;
    d_double_2d_t & tend_layerThickness = *d_tend_layerThickness;
    d_double_1d_t & surfaceThicknessFlux = *d_surfaceThicknessFlux;
    d_double_1d_t & surfaceThicknessFluxRunoff = *d_surfaceThicknessFluxRunoff;
    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
        YAKL_LAMBDA(int iCell,int k)
    {
        tend_layerThickness(k,iCell) = 0.0;
    });

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nCells) ,
        YAKL_LAMBDA(int iCell)
    {
        surfaceThicknessFlux(iCell) = 0.0;
        surfaceThicknessFluxRunoff(iCell) = 0.0;
    });

}

extern "C"
void ocn_tendency_vel_hadv_coriolis(bool usePlanetVor) {

    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesAll, mesh::nEdgesAll);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(tendency, tendVel);
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, layerThicknessEdge);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeFlux);
    YAKL_LOCAL_NS(diag_solve, kineticEnergyCell);
    YAKL_LOCAL_NS(diag_solve, bTemp);
    YAKL_LOCAL_NS(diag_solve, cTemp);
    YAKL_LOCAL_NS(diag_solve, normalizedPlanetaryVorticityEdge);
    YAKL_LOCAL_NS(diag_solve, normalizedRelativeVorticityEdge);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, nEdgesOnEdge);
    YAKL_LOCAL_NS(mesh, edgesOnEdge);
    YAKL_LOCAL_NS(mesh, weightsOnEdge);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, edgeMask);
    YAKL_LOCAL_NS(mesh, dcEdge);

    auto & tmpVorticity = bTemp;
    auto & qArr = cTemp;
    auto & normRelVortEdge = normalizedRelativeVorticityEdge;
    auto & normPlanetVortEdge = normalizedPlanetaryVorticityEdge;
    
    if ( usePlanetVor ) {
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesAll},{1,nVertLevels}) ,
        YAKL_LAMBDA(int iEdge,int k) {
            tmpVorticity(k,iEdge) = normRelVortEdge(k,iEdge) +
                                 normPlanetVortEdge(k,iEdge);
        }, yakl::DefaultLaunchConfig());
    }
    else {
        yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesAll},{1,nVertLevels}) ,
        YAKL_LAMBDA(int iEdge,int k) {
            tmpVorticity(k,iEdge) = normRelVortEdge(k,iEdge);
        }, yakl::DefaultLaunchConfig());
    }

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        qArr(k,iEdge) = 0.0;
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            for ( int j = 1; j <= nEdgesOnEdge(iEdge); ++j ) {
                double eoe = edgesOnEdge(j, iEdge);
                double edgeWeight = weightsOnEdge(j, iEdge);
                double avgVorticity = 0.5 *
                              (tmpVorticity(k,iEdge) +
                               tmpVorticity(k,eoe));
                qArr(k,iEdge) = qArr(k,iEdge) +
                               edgeWeight*normalVelocity(k,eoe)*
                               avgVorticity*layerThickEdgeFlux(k,eoe);
            }
        }
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            double invLength = 1.0 / dcEdge(iEdge);
            tendVel(k,iEdge) = tendVel(k,iEdge) +
                      edgeMask(k,iEdge)* (qArr(k,iEdge) -
                       (kineticEnergyCell(k,cell2)
                      - kineticEnergyCell(k,cell1))*invLength);
        }
    }, yakl::DefaultLaunchConfig());

    //yakl_update_host(tendency::tendVel, GET_TPTR(tendVel));

}

extern "C"
void ocn_tendency_vel_vadv() {

    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(tendency, tendVel);
    YAKL_LOCAL_NS(tendency, w_dudzTopEdge);
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, vertAleTransportTop);
    YAKL_LOCAL_NS(diag_solve, layerThickEdgeFlux);
    YAKL_LOCAL_NS(diag_solve, kineticEnergyCell);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, edgeMask);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned) ,
    YAKL_LAMBDA(int iEdge) {
        w_dudzTopEdge(minLevelEdgeBot(iEdge),iEdge) = 0.0;
        w_dudzTopEdge(maxLevelEdgeTop(iEdge)+1,iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        if ( (k > minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            
            // Average w from cell center to edge
            double wAvg = 0.5 * (vertAleTransportTop(k,cell1) +
                              vertAleTransportTop(k,cell2));

            // compute dudz at vertical interface with first order derivative.
            w_dudzTopEdge(k,iEdge) = wAvg *
                                     (normalVelocity(k-1,iEdge) -
                                      normalVelocity(k,  iEdge))/
                         (0.5 * (layerThickEdgeFlux(k-1,iEdge) +
                                     layerThickEdgeFlux(k  ,iEdge)));
        }
    }, yakl::DefaultLaunchConfig());

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            tendVel(k,iEdge) = tendVel(k,iEdge) - edgeMask(k,iEdge)*
                         0.5 * (w_dudzTopEdge(k  ,iEdge) +
                                    w_dudzTopEdge(k+1,iEdge));
        }
    }, yakl::DefaultLaunchConfig());

    //yakl_update_host(tendency::tendVel, GET_TPTR(tendVel));
    //yakl::fence();
}


extern "C"
void ocn_tendency_vel_pgrad_jacobts(int indxT, int indxS, double pGradLvlWgt,
                                    double density0Inv, double gdensity0Inv)  {

    YAKL_SCOPE(nVertLevels, mesh::nVertLevels);
    YAKL_SCOPE(nEdgesOwned, mesh::nEdgesOwned);

    YAKL_LOCAL_NS(tendency, tendVel);
    YAKL_LOCAL_NS(tendency, JacobianTz);
    YAKL_LOCAL_NS(tendency, JacobianSz);

    YAKL_LOCAL_NS(diag_solve, bTemp);
    YAKL_LOCAL_NS(diag_solve, pressure);
    YAKL_LOCAL_NS(diag_solve, zMid);
    YAKL_LOCAL_NS(diag_solve, density);
    YAKL_LOCAL_NS(diag_solve, activeTracers);
    YAKL_LOCAL_NS(diag_solve, thermExpCoeff);
    YAKL_LOCAL_NS(diag_solve, salineContractCoeff);

    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, edgeMask);
    YAKL_LOCAL_NS(mesh, dcEdge);

    auto & tracers = activeTracers;
    
    auto & pGrad = bTemp;

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned) ,
    YAKL_LAMBDA(int iEdge) {
        int k = minLevelEdgeBot(iEdge);
        JacobianTz(k,iEdge) = 0.0;
        JacobianSz(k,iEdge) = 0.0;

        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        double invdcEdge = 1.0 / dcEdge(iEdge);

        pGrad(k, iEdge) = edgeMask(k,iEdge)*invdcEdge*(
                    - density0Inv*(pressure(k,cell2) -
                                   pressure(k,cell1))
                    - gdensity0Inv*0.5 *
                      (density(k,cell1)+density(k,cell2))*
                      (zMid(k,cell2)- zMid(k,cell1) ) );

        tendVel(k,iEdge) = tendVel(k,iEdge) + pGrad(k, iEdge);
//        ++k;
//        int kmax = maxLevelEdgeTop(iEdge);
//        for ( ; k <= kmax; ++k ) {
//            pGrad(k, iEdge) += gdensity0Inv*JacobianDxDs(k,iEdge)*invdcEdge
//        }
    }, yakl::DefaultLaunchConfig());


    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        if ( (k > minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            int cell1 = cellsOnEdge(1,iEdge);
            int cell2 = cellsOnEdge(2,iEdge);
            
           // eqn 2.7 in Shchepetkin and McWilliams (2003)
           // Note delta x was removed.  It must be an error in the
           // paper, ! as it makes the units incorrect.
           double Area = 0.5 * (zMid(k-1,cell1) - zMid(k,cell1) +
                             zMid(k-1,cell2) - zMid(k,cell2));

           // eqn 2.8
           double zStar = ( zMid(k-1,cell2)*zMid(k-1,cell1) -
                     zMid(k  ,cell2)*zMid(k  ,cell1) )/
                   ( zMid(k-1,cell2)-zMid(k  ,cell2) +
                     zMid(k-1,cell1)-zMid(k,cell1));

           // eqn 3.2
           double zC = 0.25 * (zMid(k,cell1) + zMid(k-1,cell1) +
                            zMid(k,cell2) + zMid(k-1,cell2));

           // eqn 4.1
           double zGamma = (1.0 - pGradLvlWgt)*zStar + pGradLvlWgt*zC;


           double TL = (tracers(indxT,k  ,cell1)*(zMid(k-1,cell1)-zGamma) +
                 tracers(indxT,k-1,cell1)*(zGamma-zMid(k  ,cell1)))/
                (zMid(k-1,cell1) - zMid(k,cell1));
           double TR = (tracers(indxT,k  ,cell2)*(zMid(k-1,cell2)-zGamma) +
                 tracers(indxT,k-1,cell2)*(zGamma-zMid(k  ,cell2)))/
                (zMid(k-1,cell2) - zMid(k,cell2));

           double SL = (tracers(indxS,k  ,cell1)*(zMid(k-1,cell1)-zGamma) +
                 tracers(indxS,k-1,cell1)*(zGamma-zMid(k  ,cell1)))/
                (zMid(k-1,cell1) - zMid(k,cell1));
           double SR = (tracers(indxS,k  ,cell2)*(zMid(k-1,cell2)-zGamma) +
                 tracers(indxS,k-1,cell2)*(zGamma-zMid(k  ,cell2)))/
                (zMid(k-1,cell2) - zMid(k,cell2));


           // eqn 2.6 in Shchepetkin and McWilliams (2003)
           JacobianTz(k,iEdge) = Area*(TL - TR);
           JacobianSz(k,iEdge) = Area*(SL - SR);
        }
    }, yakl::DefaultLaunchConfig());


    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdgesOwned) ,
    YAKL_LAMBDA(int iEdge) {
        int kmax = maxLevelEdgeTop(iEdge);
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        double invdcEdge = 1.0 / dcEdge(iEdge);
        for ( int k = minLevelEdgeBot(iEdge) + 1; k <= kmax; ++k ) {
            double alpha = 0.25 * (
                       density(k  ,cell1)*thermExpCoeff (k  ,cell1)
                     + density(k-1,cell1)*thermExpCoeff (k-1,cell1)
                     + density(k  ,cell2)*thermExpCoeff (k  ,cell2)
                     + density(k-1,cell2)*thermExpCoeff (k-1,cell2) );
            double beta  = 0.25 * (
                       density(k  ,cell1)*salineContractCoeff(k  ,cell1)
                     + density(k-1,cell1)*salineContractCoeff(k-1,cell1)
                     + density(k  ,cell2)*salineContractCoeff(k  ,cell2)
                     + density(k-1,cell2)*salineContractCoeff(k-1,cell2) );
            double jdd = -alpha*JacobianTz(k,iEdge) +
                                         beta*JacobianSz(k,iEdge);
            pGrad(k, iEdge) = pGrad(k-1, iEdge) + gdensity0Inv*jdd*invdcEdge;
        }
    }, yakl::DefaultLaunchConfig());


    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdgesOwned},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k) {
        if ( (k > minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) ) {
            tendVel(k,iEdge) += pGrad(k,iEdge);
        }
    }, yakl::DefaultLaunchConfig());

    yakl_update_host(tendency::tendVel, GET_TPTR(tendVel));
    yakl::fence();
}

