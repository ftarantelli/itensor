#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using std::vector;

int main()
{
int N = 50; //number of sites
Real tstep = 0.02; //time step (smaller is generally more accurate)
Real ttotal = 1.0; //total time to evolve
Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form

double g = 0.5;
double h = 0.;

//Define a site set object "sites" which lets us
//easily obtain Site indices defining our Hilbert space
//and S=1/2 single-site operators
auto sites = SpinHalf(N, {"ConserveQNs=", false});

//Make initial MPS psi to be in the Neel state
auto state = InitState(sites);
for(auto j : range1(N))
    {
    state.set(j,j%2==1?"Up":"Dn");
    }
auto psi = MPS(state);

//Create a std::vector (dynamically sizeable array)
//to hold the Trotter gates
auto gates = vector<BondGate>();

//Create the gates exp(-i*tstep/2*hterm)
//and add them to gates
for(int b = 1; b <= N-1; ++b)
    {
    auto hterm = - op(sites,"Sx",b)*op(sites,"Sx",b+1);
    hterm -= g*op(sites,"Sz",b);
    hterm -= h*op(sites,"Sx",b);

    auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
    gates.push_back(g);
    }
    
//Create the gates exp(-i*tstep/2*hterm) in reverse
//order (to get a second order Trotter breakup which
//does a time step of "tstep") and add them to gates
for(int b = N-1; b >= 1; --b)
    {
    auto hterm = - op(sites,"Sx",b)*op(sites,"Sx",b+1);
    hterm -= g*op(sites,"Sz",b);
    hterm -= h*op(sites,"Sx",b);
    
    auto g = BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm);
    gates.push_back(g);
    }

//Save initial state;
auto psi0 = psi;

//Time evolve, overwriting psi when done
gateTEvol(gates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});

printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));

//Print overlap of final state with initial state
//(Will be complex so using innerC which can return complex);
Print(innerC(psi,psi0));

/*
// MEASURE
auto j = int(N/2);

psi.position(j);
//   for( auto j : range1(N) ) 
//        {


auto ket = psi(j);
auto bra = dag(prime(ket,"Site"));

auto Szjop = op(sites,"Sz",j);
//take an inner product 
auto szj = elt(bra*Szjop*ket);
printfln("Measure: %d %.12f",j,szj);

//Print();
*/

return 0;
}
