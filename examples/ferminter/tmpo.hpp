#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"

using namespace itensor;

auto measure(auto j, auto psi, auto sites) {

        psi.position(j);

        auto ket = psi(j);
        auto bra = dag(prime(ket,"Site"));
        
        //auto tampo = AutoMPO(sites);
        //tampo += "Sz", j;
        //auto Szjop = toMPO(tampo);
        auto NN = op(sites,"N",j);

        //take an inner product 
        auto nval = real(eltC(bra*NN*ket));
        //std::cout << elt(bra*ket) << "   Measss\n";
        //printfln("MEASURE Sz:			%d %.12f",j,szj);
        return(nval);
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
        auto Njop = op(sites,"N",j);

        //take an inner product 
        auto nj = eltC(bra*Njop*ket);
        //std::cout << elt(bra*ket) << "   Measss\n";
        //printfln("MEASURE Sz:			%d %.12f",j,szj);
        out += real(nj);
        } 
        //std::cout << elt(bra*ket) << "   Measss\n";
        //printfln("MEASURE Sz:			%d %.12f",j,szj);
        return(out);
}



auto correlator(auto i, auto j, auto psi, auto sites, auto N ){
// ###IMPORTANT###		ONLY FOR I \neq J

auto Adag_i = op(sites,"Cdag",i);
auto A_j = op(sites,"C",j);

//'gauge' the MPS to site i
//any 'position' between i and j, inclusive, would work here
psi.position(i);

auto psidag = dag(psi);
psidag.prime();

//index linking i to i-1:
auto li_1 = leftLinkIndex(psi,i);
auto Cij = prime(psi(i),li_1)*Adag_i*psidag(i);
for(int k = i+1; k < j; ++k)
    {
    Cij *= psi(k);
    Cij *= op(sites,"F",k); //Jordan-Wigner string
    Cij *= psidag(k);
    }
//index linking j to j+1:
auto lj = rightLinkIndex(psi,j);
Cij *= prime(psi(j),lj);
Cij *= A_j;
Cij *= psidag(j);

auto result = real(eltC(Cij)); //or eltC(C) if expecting complex
	/*
		auto op_i = op(sites,"Cdag",i);
	auto op_j = op(sites,"C",j);
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

	result += 4.*elt(C); //or eltC(C) if expecting complex
}
*/

	return(result);
}





///		COMMENT

    // Set the initial state to be Neel state
    //auto state = InitState(sites);
    
    /*
    state.set(1,"Up");
    for(int i = 2; i < N; i=i+2)
        {
        if((i/2)%2==1)
            {
            state.set(i,"Dn");
            state.set(i+1,"Dn");
            }
        else
            {
            state.set(i,"Up");
            state.set(i+1,"Up");
            }
        }
    state.set(N,"Up");
    */
