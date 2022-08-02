#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"

using namespace itensor;

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
	auto op_i = op(sites,"Sx",i);
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

	auto result = 4.*elt(C); //or eltC(C) if expecting complex

	return(result);
}

/*
Measuring Properties of MPS


    expect(MPS psi, SiteSet sites, string A) -> vector<Real>
    expectC(MPS psi, SiteSet sites, string A) -> vector<Cplx>

Compute the expected value ⟨ψ|A^|ψ⟩ for an operator A on every site, returning a vector of the results. The operator can be any operator name recognized by the SiteSet, such as "Sz" for the SpinHalf SiteSet. The function expectC returns a vector of complex values and should be used when complex results are expected.


    expect(MPS psi, SiteSet sites, vector<string> ops) -> vector<vector<Real>>
    expectC(MPS psi, SiteSet sites, vector<string> ops) -> vector<vector<Cplx>>

Compute the expected value of a set of operators ops for each site of the MPS psi. Similar to calling the single-operator version of expect above but is more efficient. The returned value is a vector of vectors, one for each operator, containing the expected values on each site.


    correlationMatrix(MPS psi, SiteSet sites, string A, string B) -> vector<vector<Real>> 	     
    correlationMatrixC(MPS psi, SiteSet sites, string A, string B) -> vector<vector<Cplx>>

Given an MPS psi, a SiteSet, and two strings denoting operators computes the two-point correlation function matrix Mij=⟨ψ|A^iB^j|ψ⟩ using efficient MPS techniques. Returns the matrix M.


***Multiplying MPOs***


    nmultMPO(MPO A, MPO B, Args args = Args::global()) -> MPO

Multiply MPOs A and B, returning the results MPO. MPO tensors are multiplied one at a time from left to right and the resulting tensors are compressed using the truncation parameters (such as "Cutoff" and "MaxDim") provided through the named arguments args.

For each j, MPO tensors A(j) and B(j) must share a single site index. MPO C will contain the site indices not shared by MPOs A and B. In addition, the link indices of MPO C will have the same tags as the link indices of the MPO A.

*/

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
