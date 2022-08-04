#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include <armadillo>

using namespace arma;
using namespace itensor;

/*
int N = 10;
auto sites = SpinHalf(N,{"ConserveQNs=",false});

auto A = randomMPS(sites); //create random product MPS

// Shift MPS gauge such that site 1 is
// the orthogonality center
A.position(1);
//Shift orthogonality center to site j
auto j = 5;
A.position(j);

// Read-only access of tensor at site j
auto Aj = A(j);

// Replace tensor at site j with
// a modified tensor.
A.set(j,2*Aj);

// Directly modify tensor at site j; "ref"
// signified that a reference to A_j tensor is returned
A.ref(j) *= -1;

// Initialize an MPS to a specific product state
auto state = InitState(sites);
for(int i = 1; i <= N; ++i)
    {
    if(i%2 == 0) state.set(i,"Up");
    else         state.set(i,"Dn");
    }
auto B = MPS(state);
*/

auto ret_kron_mps(auto mps1, auto mps2, auto sites, auto nsites, auto N){

	auto res = MPS(nsites); //MPS(2*N);
	for(int i=1; i<=N; ++i){
		mps1.position(i);
		mps2.position(i);

		res.set(2*i-1, mps1(i));
		res.set(2*i, mps2(i));

		res.ref(2*i-1) *= -1;
		res.ref(2*i) *= -1;
	}

	return(res);
}



auto ret_kron_mpo(auto mpo1, auto mpo2, auto sites, auto nsites, auto N){

	auto res = MPO(nsites); //MPO(2*N);
	for(int i=1; i<=N; ++i){
		mpo1.position(i);
		mpo2.position(i);

		res.set(2*i-1, mpo1(i));
		res.set(2*i, mpo2(i));

		res.ref(2*i-1) *= -1;
		res.ref(2*i) *= -1;
	}

	return(res);
}


auto measure(auto j, auto psi, auto sites) {

        psi.position(j);

        auto ket = psi(j);
        auto bra = dag(prime(ket,"Site"));
        
        //auto tampo = AutoMPO(sites);
        //tampo += "Sz", j;
        //auto Szjop = toMPO(tampo);
        auto Szjop = op(sites,"Sz",j);

        //take an inner product 
        auto szj = real(2.*eltC(bra*Szjop*ket));
        //std::cout << elt(bra*ket) << "   Measss\n";
        //printfln("MEASURE Sz:			%d %.12f",j,szj);
        return(szj);
}



auto tot_meas(auto psi, auto sites, int N) {

	double out(0.);
	for(int j=1; j<=N; ++j)
	   {
        psi.position(j);

        auto ket = psi(j);
        auto bra = dag(prime(ket,"Site"));
        
        //auto tampo = AutoMPO(sites);
        //tampo += "Sz", j;
        //auto Szjop = toMPO(tampo);
        auto Szjop = op(sites,"Sz",j);

        //take an inner product 
        auto szj = 2.*eltC(bra*Szjop*ket);
        out += real(szj);
        } 
        //std::cout << elt(bra*ket) << "   Measss\n";
        //printfln("MEASURE Sz:			%d %.12f",j,szj);
        return(out/float(N));
}



auto correlator(auto i, auto j, auto psi, auto sites){
//Make the operators you want to measure
	auto op_i = op(sites,"Sz",i);
	auto op_j = op(sites,"Sz",j);
//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
	psi.position(i); 
//Create the bra/dual version of the MPS psi
	auto psidag = dag(psi);
//Prime the link indices to make them distinct from
//the original ket links
	psidag.prime("Link");
//index linking i-1 to i:
	auto li_1 = leftLinkIndex(psi,i);

	auto C = prime(psi(i),li_1)*op_i;
	C *= prime(psidag(i),"Site");
	for(int k = i+1; k < j; ++k)
 	  {
 	   C *= psi(k);
 	   C *= psidag(k);
  	  }
//index linking j to j+1:
	auto lj = rightLinkIndex(psi,j);

	C *= prime(psi(j),lj)*op_j;
	C *= prime(psidag(j),"Site");

	auto result = real(4.*eltC(C)); //or eltC(C) if expecting complex

	return(result);
}

///		COMMENT

/*
auto measure(auto j, auto psi1, auto psi2, auto sites) {

        psi1.position(j);
	psi2.position(j);

        auto ket1 = psi1(j);
	auto ket2 = psi2(j);
        auto bra2 = dag(prime(ket2,"Site"));
        
        //auto tampo = AutoMPO(sites);
        //tampo += "Sz", j;
        //auto Szjop = toMPO(tampo);
        auto Szjop = op(sites,"Sz",j);

        //take an inner product 
        auto szj = real(2.*eltC(bra2*Szjop*ket1));
        //std::cout << elt(bra*ket) << "   Measss\n";
        //printfln("MEASURE Sz:			%d %.12f",j,szj);

        return(szj);
}
*/
