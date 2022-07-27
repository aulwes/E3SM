#include "mpas_ocn_yakl_types.hxx"
#include "mpas_ocn_yakl_c.hxx"

namespace eos_jm
{
d_double_1d_t   * p,
                * p2
                ;
                
d_double_2d_t   * tracerSalt,
                * tracerTemp,
                * density,
                * thermalExpansionCoeff,
                * salineContractionCoeff
                ;
}

extern "C" real
    rhosfc_c,
    unt0_c,
    unt1_c,
    unt2_c,
    unt3_c,
    unt4_c,
    unt5_c,
    unsqt0_c,
    unsqt1_c,
    unsqt2_c,
    uns1t0_c,
    uns1t1_c,
    uns1t2_c,
    uns1t3_c,
    uns1t4_c,
    uns2t0_c,
    bup0s1t0_c,
    bup0s1t1_c,
    bup0s1t2_c,
    bup0s1t3_c,
    bup2s1t0_c,
    bup2s1t1_c,
    bup2s1t2_c,
    bup0sqt0_c,
    bup0sqt1_c,
    bup0sqt2_c,
    bup1sqt0_c,
    bup1s0t0_c,
    bup1s0t1_c,
    bup1s0t2_c,
    bup1s0t3_c,

    bup0s0t3_c,
    bup0s0t4_c,

    bup1s1t0_c,
    bup1s1t1_c,
    bup1s1t2_c,

    bup0s0t0_c,
    bup0s0t1_c,
    bup0s0t2_c,

    bup2s0t0_c,
    bup2s0t1_c,
    bup2s0t2_c;

extern "C"
void ocn_eos_density_only(int nCells, int nVertLevels,
                          double * h_p, double * h_p2,
                          double * h_tracerSalt, double * h_tracerTemp, double * h_density)
{
    yakl::timer_start("ocn_eos_density_only");

    YAKL_LOCAL_NS(eos_jm,p);
    YAKL_LOCAL_NS(eos_jm,p2);
    YAKL_LOCAL_NS(eos_jm,tracerSalt);
    YAKL_LOCAL_NS(eos_jm,tracerTemp);
    YAKL_LOCAL_NS(eos_jm,density);

    yakl_update_device(eos_jm::p, h_p);
    yakl_update_device(eos_jm::p2, h_p2);
    yakl_update_device(eos_jm::tracerSalt, h_tracerSalt);
    yakl_update_device(eos_jm::tracerTemp, h_tracerTemp);

    YAKL_SCOPE(uns2t0, uns2t0_c);

    YAKL_SCOPE(uns1t0, uns1t0_c);
    YAKL_SCOPE(uns1t1, uns1t1_c);
    YAKL_SCOPE(uns1t2, uns1t2_c);
    YAKL_SCOPE(uns1t3, uns1t3_c);
    YAKL_SCOPE(uns1t4, uns1t4_c);

    YAKL_SCOPE(unt0, unt0_c);
    YAKL_SCOPE(unt1, unt1_c);
    YAKL_SCOPE(unt2, unt2_c);
    YAKL_SCOPE(unt3, unt3_c);
    YAKL_SCOPE(unt4, unt4_c);
    YAKL_SCOPE(unt5, unt5_c);

    YAKL_SCOPE(bup0s1t0, bup0s1t0_c);
    YAKL_SCOPE(bup0s1t1, bup0s1t1_c);
    YAKL_SCOPE(bup0s1t2, bup0s1t2_c);
    YAKL_SCOPE(bup0s1t3, bup0s1t3_c);

    YAKL_SCOPE(bup0s0t0, bup0s0t0_c);
    YAKL_SCOPE(bup0s0t1, bup0s0t1_c);
    YAKL_SCOPE(bup0s0t2, bup0s0t2_c);

    YAKL_SCOPE(bup1s1t0, bup1s1t0_c);
    YAKL_SCOPE(bup1s1t1, bup1s1t1_c);
    YAKL_SCOPE(bup1s1t2, bup1s1t2_c);

    YAKL_SCOPE(bup1s0t0, bup1s0t0_c);
    YAKL_SCOPE(bup1s0t2, bup1s0t2_c);

    YAKL_SCOPE(bup2s0t0, bup2s0t0_c);
    YAKL_SCOPE(bup2s0t1, bup2s0t1_c);
    YAKL_SCOPE(bup2s0t2, bup2s0t2_c);

    YAKL_SCOPE(bup1sqt0, bup1sqt0_c);

    YAKL_SCOPE(unsqt0, unsqt0_c);
    YAKL_SCOPE(unsqt1, unsqt1_c);
    YAKL_SCOPE(unsqt2, unsqt2_c);

    YAKL_SCOPE(bup2s1t0, bup2s1t0_c);
    YAKL_SCOPE(bup2s1t1, bup2s1t1_c);
    YAKL_SCOPE(bup2s1t2, bup2s1t2_c);

    YAKL_SCOPE(bup0sqt0, bup0sqt0_c);
    YAKL_SCOPE(bup0sqt1, bup0sqt1_c);
    YAKL_SCOPE(bup0sqt2, bup0sqt2_c);

    YAKL_SCOPE(bup0s0t3, bup0s0t3_c);
    YAKL_SCOPE(bup0s0t4, bup0s0t4_c);

    YAKL_SCOPE(bup1s0t1, bup1s0t1_c);
    YAKL_SCOPE(bup1s0t3, bup1s0t3_c);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
         auto sq  = tracerSalt(k,iCell);
         auto tq  = tracerTemp(k,iCell);

         auto sqr = sqrt(sq);
         auto t2  = tq*tq;

         ///
         /// first calculate surface (p=0) values from UNESCO eqns.
         ///

         auto work1 =      uns1t0 + uns1t1*tq +
                     (uns1t2 + uns1t3*tq + uns1t4*t2)*t2;
         auto work2 = sqr*(unsqt0 + unsqt1*tq + unsqt2*t2);

         auto rhosfc = unt1*tq + (unt2 + unt3*tq +
                                    (unt4 + unt5*tq)*t2)*t2
                          + (uns2t0*sq + work1 + work2)*sq;

         ///
         /// now calculate bulk modulus at pressure p from
         /// Jackett and McDougall formula
         ///

          auto work3 = bup0s1t0 + bup0s1t1*tq +
                (bup0s1t2 + bup0s1t3*tq)*t2 +
              p(k) *(bup1s1t0 + bup1s1t1*tq + bup1s1t2*t2) +
              p2(k)*(bup2s1t0 + bup2s1t1*tq + bup2s1t2*t2);

         auto work4 = sqr*(bup0sqt0 + bup0sqt1*tq + bup0sqt2*t2 +
                      bup1sqt0*p(k));

         auto bulkMod = bup0s0t0 + bup0s0t1*tq +
                  (bup0s0t2 + bup0s0t3*tq + bup0s0t4*t2)*t2 +
            p(k) *(bup1s0t0 + bup1s0t1*tq +
                  (bup1s0t2 + bup1s0t3*tq)*t2) +
            p2(k)*(bup2s0t0 + bup2s0t1*tq + bup2s0t2*t2) +
                   sq*(work3 + work4);


         density(k,iCell) = (unt0 + rhosfc)*bulkMod/
                            (bulkMod - p(k));
    });

    yakl_update_host(eos_jm::density, h_density);

    yakl::timer_stop("ocn_eos_density_only");
}

extern "C"
void ocn_eos_density_exp(int nCells, int nVertLevels,
                        double * h_p, double * h_p2,
                        double * h_tracerSalt, double * h_tracerTemp, double * h_density,
                        double * h_thermalExpansionCoeff, double * h_salineContractionCoeff)
{
    YAKL_LOCAL_NS(eos_jm,p);
    YAKL_LOCAL_NS(eos_jm,p2);
    YAKL_LOCAL_NS(eos_jm,tracerSalt);
    YAKL_LOCAL_NS(eos_jm,tracerTemp);
    YAKL_LOCAL_NS(eos_jm,density);
    YAKL_LOCAL_NS(eos_jm,thermalExpansionCoeff);
    YAKL_LOCAL_NS(eos_jm,salineContractionCoeff);

    yakl_update_device(eos_jm::p, h_p);
    yakl_update_device(eos_jm::p2, h_p2);
    yakl_update_device(eos_jm::tracerSalt, h_tracerSalt);
    yakl_update_device(eos_jm::tracerTemp, h_tracerTemp);

    YAKL_SCOPE(uns2t0, uns2t0_c);

    YAKL_SCOPE(uns1t0, uns1t0_c);
    YAKL_SCOPE(uns1t1, uns1t1_c);
    YAKL_SCOPE(uns1t2, uns1t2_c);
    YAKL_SCOPE(uns1t3, uns1t3_c);
    YAKL_SCOPE(uns1t4, uns1t4_c);

    YAKL_SCOPE(unt0, unt0_c);
    YAKL_SCOPE(unt1, unt1_c);
    YAKL_SCOPE(unt2, unt2_c);
    YAKL_SCOPE(unt3, unt3_c);
    YAKL_SCOPE(unt4, unt4_c);
    YAKL_SCOPE(unt5, unt5_c);

    YAKL_SCOPE(bup0s1t0, bup0s1t0_c);
    YAKL_SCOPE(bup0s1t1, bup0s1t1_c);
    YAKL_SCOPE(bup0s1t2, bup0s1t2_c);
    YAKL_SCOPE(bup0s1t3, bup0s1t3_c);

    YAKL_SCOPE(bup0s0t0, bup0s0t0_c);
    YAKL_SCOPE(bup0s0t1, bup0s0t1_c);
    YAKL_SCOPE(bup0s0t2, bup0s0t2_c);

    YAKL_SCOPE(bup1s1t0, bup1s1t0_c);
    YAKL_SCOPE(bup1s1t1, bup1s1t1_c);
    YAKL_SCOPE(bup1s1t2, bup1s1t2_c);

    YAKL_SCOPE(bup1s0t0, bup1s0t0_c);
    YAKL_SCOPE(bup1s0t2, bup1s0t2_c);

    YAKL_SCOPE(bup2s0t0, bup2s0t0_c);
    YAKL_SCOPE(bup2s0t1, bup2s0t1_c);
    YAKL_SCOPE(bup2s0t2, bup2s0t2_c);

    YAKL_SCOPE(bup1sqt0, bup1sqt0_c);

    YAKL_SCOPE(unsqt0, unsqt0_c);
    YAKL_SCOPE(unsqt1, unsqt1_c);
    YAKL_SCOPE(unsqt2, unsqt2_c);

    YAKL_SCOPE(bup2s1t0, bup2s1t0_c);
    YAKL_SCOPE(bup2s1t1, bup2s1t1_c);
    YAKL_SCOPE(bup2s1t2, bup2s1t2_c);

    YAKL_SCOPE(bup0sqt0, bup0sqt0_c);
    YAKL_SCOPE(bup0sqt1, bup0sqt1_c);
    YAKL_SCOPE(bup0sqt2, bup0sqt2_c);

    YAKL_SCOPE(bup0s0t3, bup0s0t3_c);
    YAKL_SCOPE(bup0s0t4, bup0s0t4_c);

    YAKL_SCOPE(bup1s0t1, bup1s0t1_c);
    YAKL_SCOPE(bup1s0t3, bup1s0t3_c);

    yakl::fortran::parallel_for( yakl::fortran::Bounds<2>({1,nCells},{1,nVertLevels}) ,
    YAKL_LAMBDA(int iCell,int k)
    {
        real sq  = tracerSalt(k,iCell);
        real tq  = tracerTemp(k,iCell);

        real sqr = std::sqrt(sq);
        real t2  = tq*tq;

        /***
         *** first calculate surface (p=0) values from UNESCO eqns.
         ***/

        real work1 = uns1t0 + uns1t1*tq +
                    (uns1t2 + uns1t3*tq + uns1t4*t2)*t2;
        real work2 = sqr*(unsqt0 + unsqt1*tq + unsqt2*t2);

        real rhosfc = unt1*tq + (unt2 + unt3*tq + (unt4 + unt5*tq)*t2)*t2
                        + (uns2t0*sq + work1 + work2)*sq;

        /***
        *** now calculate bulk modulus at pressure p from
        *** Jackett and McDougall formula
        ***/

        real work3 = bup0s1t0 + bup0s1t1*tq +
               (bup0s1t2 + bup0s1t3*tq)*t2 +
         p(k) *(bup1s1t0 + bup1s1t1*tq + bup1s1t2*t2) +
         p2(k)*(bup2s1t0 + bup2s1t1*tq + bup2s1t2*t2);

        real work4 = sqr*(bup0sqt0 + bup0sqt1*tq + bup0sqt2*t2 +
                          bup1sqt0*p(k));

        real bulkMod  = bup0s0t0 + bup0s0t1*tq +
                  (bup0s0t2 + bup0s0t3*tq + bup0s0t4*t2)*t2 +
            p(k) *(bup1s0t0 + bup1s0t1*tq +
                  (bup1s0t2 + bup1s0t3*tq)*t2) +
            p2(k)*(bup2s0t0 + bup2s0t1*tq + bup2s0t2*t2) +
                    sq*(work3 + work4);

        /***
        *** compute density
        ***/

        real denomk = 1.0/(bulkMod - p(k));

        density(k, iCell) = (unt0 + rhosfc)*bulkMod*denomk;

        /***
        *** compute temperature expansion coeff
        ***  by differentiating above formulae
        ***/

        real drdt0 = unt1 + 2.0*unt2*tq +
                 (3.0*unt3 + 4.0*unt4*tq +
                                   5.0*unt5*t2)*t2 +
                         (uns1t1 + 2.0*uns1t2*tq +
               (3.0*uns1t3 + 4.0*uns1t4*tq)*t2 +
                          (unsqt1 + 2.0*unsqt2*tq)*sqr )*sq;

        real dkdt  = bup0s0t1 + 2.0*bup0s0t2*tq +
                (3.0*bup0s0t3 + 4.0*bup0s0t4*tq)*t2 +
                    p(k) *(bup1s0t1 + 2.0*bup1s0t2*tq +
                                      3.0*bup1s0t3*t2) +
                    p2(k)*(bup2s0t1 + 2.0*bup2s0t2*tq) +
                       sq*(bup0s1t1 + 2.0*bup0s1t2*tq +
                                      3.0*bup0s1t3*t2 +
                   p(k)  *(bup1s1t1 + 2.0*bup1s1t2*tq) +
                   p2(k) *(bup2s1t1 + 2.0*bup2s1t2*tq) +
                           sqr*(bup0sqt1 + 2.0*bup0sqt2*tq));

        real drhodt = (denomk*(drdt0*bulkMod -
                               p(k)*(unt0+rhosfc)*dkdt*denomk));

        real invdens = 1.0 / density(k,iCell);
        thermalExpansionCoeff(k,iCell) = -drhodt * invdens;

        /***
        *** compute salinity contraction coeff
        ***  by differentiating above formulae
        ***/

        real drds0  = 2.0*uns2t0*sq + work1 + 1.5*work2;
        real dkds   = work3 + 1.5*work4;

        real drhods = denomk*(drds0*bulkMod -
                              p(k)*(unt0+rhosfc)*dkds*denomk);

        salineContractionCoeff(k,iCell) = drhods*invdens;
    });

    yakl_update_host(eos_jm::density, h_density);
    yakl_update_host(eos_jm::salineContractionCoeff, h_salineContractionCoeff);
    yakl_update_host(eos_jm::thermalExpansionCoeff, h_thermalExpansionCoeff);
}

extern "C"
void ocn_yakl_eos_jm_init(int nCells, int nVertLevels)
{
    using namespace eos_jm;
    
    p = yakl_create_real("p", nVertLevels);
    p2 = yakl_create_real("p2", nVertLevels);
    //std::cerr << "ocn_yakl_eos_jm_init: nVertLevels, nCells = " << nVertLevels << " " << nCells << std::endl;
    tracerTemp = yakl_create_real("tracerTemp", nVertLevels, nCells);
    tracerSalt = yakl_create_real("tracerSalt", nVertLevels, nCells);
    density = yakl_create_real("density", nVertLevels, nCells);
    thermalExpansionCoeff = yakl_create_real("thermalExpansionCoeff", nVertLevels, nCells);
    salineContractionCoeff = yakl_create_real("salineContractionCoeff", nVertLevels, nCells);
}