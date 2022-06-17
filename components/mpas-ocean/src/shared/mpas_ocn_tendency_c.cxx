#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"

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

    yakl::fortran::parallel_for( yakl::fortran::Bounds<1>({1,nCells}) ,
        YAKL_LAMBDA(int iCell)
    {
        surfaceThicknessFlux(iCell) = 0.0;
        surfaceThicknessFluxRunoff(iCell) = 0.0;
    });

}

