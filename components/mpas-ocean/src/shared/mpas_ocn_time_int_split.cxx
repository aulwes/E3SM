#include "mpas_ocn_diagnostics_yakl_c.hxx"
#include "mpas_ocn_yakl_c.hxx"
#include "mpas_ocn_mesh_c.hxx"
#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_diagnostics_yakl_c.hxx"

#define GET_DPTR(v) (static_cast<double *>((timeint_split::c_##v).ptr))

namespace timeint_split
{
yakl::Stream     stream1;
yakl::Stream     stream2;

d_double_1d_t   * barotropicForcing
                ;

d_double_2d_t   * normalVelocityTend,
                * normalVelocityNew,
                * normalBaroclinicVelocity,
                * normalBaroclinicVelocityCur,
                * normalBaroclinicVelocityNew
                ;

};

extern "C" ocn_yakl_type c_barotropicForcing;
extern "C" ocn_yakl_type c_normalBaroclinicVelocity;
extern "C" ocn_yakl_type c_normalVelocityTend;

extern "C"
void ocn_time_int_init()
{    
    using namespace timeint_split;
    

    normalVelocityTend = yakl_create_real("normalVelocityTend", 
                    c_normalVelocityTend.shape[0], c_normalVelocityTend.shape[1]);
    normalVelocityNew = yakl_create_real("normalVelocityNew", 
                    c_normalVelocityTend.shape[0], c_normalVelocityTend.shape[1]);
    barotropicForcing = yakl_create_real("barotropicForcing", 
                    c_barotropicForcing.shape[0]);
    normalBaroclinicVelocity = yakl_create_real("normalBaroclinicVelocity", 
                    c_normalBaroclinicVelocity.shape[0], c_normalBaroclinicVelocity.shape[1]);
    normalBaroclinicVelocityNew = yakl_create_real("normalBaroclinicVelocityNew", 
                    c_normalBaroclinicVelocity.shape[0], c_normalBaroclinicVelocity.shape[1]);
    normalBaroclinicVelocityCur = yakl_create_real("normalBaroclinicVelocityCur", 
                    c_normalBaroclinicVelocity.shape[0], c_normalBaroclinicVelocity.shape[1]);
}

extern "C"
void ocn_timeint_pre_baroclinic(double * h_ssh,
                           double * h_layerThicknessEdge)
{
    using diag_solve::stream1;

    yakl_update_device(diag_solve::ssh, h_ssh, stream1);
    yakl_update_device(diag_solve::layerThicknessEdge, h_layerThicknessEdge, stream1);
}

extern "C"
void ocn_timeint_updvel(double * h_normalVelocityTend)
{
    using diag_solve::stream1;

    yakl_update_device(timeint_split::normalVelocityTend, h_normalVelocityTend, stream1);
}

extern "C"
void ocn_timeint_postfuperp()
{
    using diag_solve::stream1;

    double * hssh = new double[diag_solve::ssh->extent(0)];
    yakl_update_host(diag_solve::ssh, hssh, stream1);
    yakl::fence();
    
    delete[] hssh;
}

extern "C"
void ocn_diag_solve_fuperp(double splitFact, double gravity, double dt, double * h_barotropicForcing,
                           double * h_normalBaroclinicVelocityCur,
                           double * h_normalBaroclinicVelocityNew, double * h_normalVelocity)
{
    using diag_solve::stream1;
    using diag_solve::stream2;

    yakl_update_device(timeint_split::normalBaroclinicVelocityCur, h_normalBaroclinicVelocityCur, stream1);
    yakl_update_device(timeint_split::normalBaroclinicVelocityNew, h_normalBaroclinicVelocityNew, stream1);

    YAKL_LOCAL_NS(diag_solve, ssh);
    YAKL_LOCAL_NS(diag_solve, normalVelocity);
    YAKL_LOCAL_NS(diag_solve, layerThicknessEdge);
    YAKL_LOCAL_NS(timeint_split, normalVelocityTend);
    YAKL_LOCAL_NS(timeint_split, barotropicForcing);
    YAKL_LOCAL_NS(timeint_split, normalBaroclinicVelocityCur);  // this is really normalBaroclinicVelocityCur
    YAKL_LOCAL_NS(timeint_split, normalBaroclinicVelocityNew);
    YAKL_LOCAL_NS(mesh, minLevelEdgeBot);
    YAKL_LOCAL_NS(mesh, maxLevelEdgeTop);
    YAKL_LOCAL_NS(mesh, nEdgesOnEdge);
    YAKL_LOCAL_NS(mesh, cellsOnEdge);
    YAKL_LOCAL_NS(mesh, dcEdge);
    YAKL_LOCAL_NS(mesh, edgesOnEdge);
    YAKL_LOCAL_NS(mesh, weightsOnEdge);
    YAKL_LOCAL_NS(mesh, fEdge);
    YAKL_SCOPE(nEdges, diag_solve::nEdges);
    YAKL_SCOPE(nVertLevels, diag_solve::nVertLevels);


    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        // initialize layerThicknessEdge to avoid divide by zero and NaN problems.
        normalVelocity(k, iEdge) = 0.0;
    }, yakl::DefaultLaunchConfig().set_stream(stream1));

    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) )
        {
          int cell1 = cellsOnEdge(1,iEdge);
          int cell2 = cellsOnEdge(2,iEdge);
          double nvel = 0.0;
          for (int j = 1; j <= nEdgesOnEdge(iEdge); ++j )
          {
               int eoe = edgesOnEdge(j,iEdge);
               nvel = nvel + weightsOnEdge(j,iEdge) * normalBaroclinicVelocityNew(k,eoe)
                                        * fEdge(eoe);
          }
          normalVelocity(k,iEdge) = nvel;
        }
    }, yakl::DefaultLaunchConfig().set_stream(stream1));
    yakl_update_host(diag_solve::normalVelocity, h_normalVelocity, stream1);



    double splgr = splitFact*gravity;
    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>(nEdges) ,
    YAKL_LAMBDA(int iEdge)
    {
        int k = minLevelEdgeBot(iEdge);
        int cell1 = cellsOnEdge(1,iEdge);
        int cell2 = cellsOnEdge(2,iEdge);
        double sshdiff = ssh(cell2) - ssh(cell1);
        double invdc = 1.0 / dcEdge(iEdge);
        double fac = splgr * sshdiff * invdc;
        
        double normalThicknessFluxSum = 0.0;
        double thicknessSum = 0.0;
        int minl = minLevelEdgeBot(iEdge);
        int maxl = maxLevelEdgeTop(iEdge);
        for ( k = minl; k <= maxl; ++k )
        {
              double uTemp = normalBaroclinicVelocityCur(k,iEdge)
                           + dt * (normalVelocityTend(k,iEdge)
                           + normalVelocity(k,iEdge)
                           + fac );
              normalThicknessFluxSum = normalThicknessFluxSum +
                             layerThicknessEdge(k,iEdge)*uTemp;
              thicknessSum = thicknessSum + layerThicknessEdge(k,iEdge);
        }
        barotropicForcing(iEdge) = 0.0;
        if ( maxl >= minl )
        {
            barotropicForcing(iEdge) = splitFact*normalThicknessFluxSum/thicknessSum/dt;
        }

    
    }, yakl::DefaultLaunchConfig().set_stream(stream1));
    yakl_update_host(timeint_split::barotropicForcing, h_barotropicForcing, stream1);

    
    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nEdges},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iEdge,int k)
    {
        if ( (k >= minLevelEdgeBot(iEdge)) && (k <= maxLevelEdgeTop(iEdge)) )
        {
          int cell1 = cellsOnEdge(1,iEdge);
          int cell2 = cellsOnEdge(2,iEdge);
          
          double uTemp = normalBaroclinicVelocityCur(k,iEdge)
                           + dt * (normalVelocityTend(k,iEdge)
                           + normalVelocity(k,iEdge)
                           + splitFact*gravity *
                             (ssh(cell2) - ssh(cell1))
                             /dcEdge(iEdge) );
          // These two steps are together here:
          // {\bf u}'_{k,n+1} =
          //    {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
          // {\bf u}'_{k,n+1/2} = \frac{1}{2}\left(
          //        {\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right)
          // so that normalBaroclinicVelocityNew is at time n+1/2
          normalBaroclinicVelocityNew(k,iEdge) = 0.5 * (
          normalBaroclinicVelocityCur(k,iEdge) + uTemp -
                 dt * barotropicForcing(iEdge));

        }
    }, yakl::DefaultLaunchConfig().set_stream(stream1));

    yakl_update_host(timeint_split::normalBaroclinicVelocityNew, h_normalBaroclinicVelocityNew, stream1);
    
    //yakl::fence();
}

