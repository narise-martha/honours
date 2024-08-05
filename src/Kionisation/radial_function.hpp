#include "Kion_functions.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/DiagramRPA0_jL.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Methods.hpp"
#include "qip/omp.hpp"
#include "Angular/Wigner369j.hpp"
#include <iostream>
#include <memory>

// This function integrates
// int func1*jL*func2
double radial_int(const std::vector<double> &func1, const std::vector<double> &func2, const std::vector<double> &jL, Grid grid){
    return NumCalc::integrate(grid.du(),0,0,grid.drdu(),func1,func2,jL);
}

// Angular integral
// C1 -> l_a = Tilde{l}_a
// C2 -> l_b = Tilde{l}_b
// m_a = m_b = 1/2
double angular(int twoj_a, int l_a, int twoj_b, int l_b, int L){

    // First part
    double first_comp = pow(-1,l_b)
    *Angular::threej_2(int(2*l_b),2*L,int(2*l_a),0,0,0)
    *Angular::cg_2(2*l_b,1,twoj_b,0,1,1)
    *Angular::cg_2(2*l_a,1,twoj_a,0,1,1);

    // Second part
    double second_comp = pow(-1,l_b-1)
    *Angular::threej_2(int(2*l_b),2*L,int(2*l_a),-2,0,2)
    *Angular::cg_2(2*l_b,1,twoj_b,2,-1,1)
    *Angular::cg_2(2*l_a,1,twoj_a,2,-1,1);

    return pow(-1,l_b)*sqrt((2*l_b + 1)*(2*l_a + 1))
    *Angular::threej_2(twoj_b,2*L,twoj_a,0,0,0)*(first_comp - second_comp);
}

// Calculates the form factor K for a single final state. 
double K_single(const DiracSpinor &Fa, const DiracSpinor &Fb, double q, int L, const auto jL){

    // Defining the quantum numbers to feed into the angular integrals
    auto l_a = Fa.l();
    auto twoj_a =  Fa.twoj();
    auto l_til_a = Angular::l_tilde_k(Fa.kappa());

    auto l_b = Fb.l();
    auto twoj_b = Fb.twoj();
    auto l_til_b = Angular::l_tilde_k(Fb.kappa());

    double C1 = angular(twoj_a,l_til_a,twoj_b,l_b,L);
    double C2 = angular(twoj_a,l_a,twoj_a,l_til_b,L);

    const auto &grid = Fa.grid();
    using namespace qip::overloads;

    // Sphereical bessel function. Need to take this out of the function
    // To reduce computation time. Generate all the Bessel functions outside of the function,
    // and then feed them into the function as an imput.
    // const auto jL = SphericalBessel::fillBesselVec(L, grid.r()*q);

    return C1*radial_int(Fb.f(),Fa.g(),jL,grid) + C2*radial_int(Fb.g(),Fa.f(),jL,grid);

}

// This function calculates the form factor K, summing over all the final states Fb
// for which the final energy is positive (results in an ionised state)
double K_total(const HF::HartreeFock *vHF, const DiracSpinor &Fa, double mass_z, double q, int L, const auto jL){
    ContinuumOrbitals cntm(vHF);
    
    // Calculating final energy state. 
    // E_f = mc^2 - |E_i| = mc^2 + E_i
    // since E_i < 0
    const auto ec = mass_z + Fa.en(); // mass_z is technically mc^2
    
    // If final energy is negative, ignore this state
    // Only want ionised final states. 
    if(ec < 0.0){
        return 0.0;
    }

    // Solving for continuum states
    cntm.solveContinuumHF(ec, 0, 4, &Fa);

    // Calculating the form factor of a single state, and then 
    // summing over the final states states
    double K = 0.0;
    for (const auto & Fb : cntm.orbitals ){
       K += K_single(Fa, Fb, q, L, jL);
    }
    return K;
}