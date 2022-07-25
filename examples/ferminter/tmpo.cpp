#include "tmpo.hpp"
#include <fstream>




int main(int argc, char *argv[])
    {

	double pbc(-1.), t0(0.01), tend(1.), t(0.05);
	int N = 20;
	auto sites = Spinless(N, {"ConserveQNs=", false} );
	
	const double yg(2.), ymu(1.);
	double kappa(1.), jack(1.);
	
while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'N':
				N = atoi( &argv[1][1] );
			break;
		case 't':
			if(argv[1][1] == 'o')		t0 = atof( &argv[1][2] );
			else if(argv[1][1] == 'd')	tend = atof( &argv[1][2] );
			else						t = atof( &argv[1][1] );			
			break;
		case 'j':
				jack = atof( &argv[1][1] );
			break;
		case 'k':
				kappa = atof( &argv[1][1] );
			break;
		case 'p':
				pbc = atof( &argv[1][1] );
			break;
		case '-':
			//	std::sprintf(folder, "%s", &argv[1][1]);
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			exit (8);
	}
	++argv;
	--argc;
}

	auto mu = -2. + kappa * std::pow(N, -ymu) ;
	auto delta = 1.;
	auto g = jack * std::pow(N, -yg);
	
	auto ampo = AutoMPO(sites);
	for( auto n : range1(N-1) )
   	 {
   	 ampo += -delta,"Cdag",n,"Cdag",n+1;
   	 ampo += -delta,"C",n+1,"C",n;
   	 
   	 ampo += -1.,"Cdag",n,"C",n+1;
   	 ampo += -1.,"Cdag",n+1,"C",n;
   	 
   	 ampo += -mu,"Cdag",n,"C",n;
   	 ampo += -g,"Cdag",n,"C",n, "Cdag",n+1,"C",n+1;
   	 
   	 //ampo += - Cplx_i / 2., "Cdag", n, "C", n;
   	 }
   	 ampo += -mu,"Cdag",N,"C",N;

   	 ampo += -pbc*delta,"Cdag",N,"Cdag",1;
   	 ampo += -pbc*delta,"C",1,"C",N;
   	 
   	 ampo += -pbc,"Cdag",N,"C",1;
   	 ampo += -pbc,"Cdag",1,"C",N;

   	 //ampo += - Cplx_i / 2., "Cdag", N, "C", N;

	auto H = toMPO(ampo);
	//auto state = InitState(sites,"Emp");
	//for( auto n : range1(N) ) if( n%2==0 ) state.set(n,"Occ");
	auto sweeps0 = Sweeps(5); //number of sweeps is 5
	sweeps0.maxdim() = 10,20,100,100,200;
	sweeps0.cutoff() = 1E-10;
	auto [energy,psi] = dmrg(H,randomMPS(sites),sweeps0,{"Quiet=",true});
	
	std::cout << energy << "		" << tot_meas(psi, sites, N) << "	Energy\n";

    char output[80];
    std::sprintf(output, "couplingN%dk%.0fj%.0f.dat", N, kappa, jack);

    std::ofstream out_file(output, std::ios::out | std::ios::trunc);
    out_file.precision(16);
    
    auto psi1 = MPS(psi);
    auto psi2 = psi1;
    

   // auto expH1 = toExpH(ampo,-(1-1_i)/2*t0*Cplx_i);
    //auto expH2 = toExpH(ampo,-(1+1_i)/2*t0*Cplx_i);
    //printfln("Maximum bond dimension of expH1 is %d",maxLinkDim(expH1));
    auto args = Args("Method=","DensityMatrix","Cutoff=",1E-12,"MaxDim=",2000);
    for(int n = 1; n <= int(t/t0); ++n)
        {
        psi2 = applyExp(H,psi2,-t0*Cplx_i);
       // psi2.noPrime();
        //psi2 = applyMPO(expH2,psi2,args);
        psi2.noPrime().normalize();
        if(n%int(std::real(t/t0)) == 0)
            {
            //printfln("\nMaximum bond dimension at time %.1f is %d ", n*t0, maxLinkDim(psi2));
            //printfln("Energy using overlap at time %.1f is %.10f", n*t0, real(innerC(psi2,H,psi2)) );
            
            
            out_file << n*t0 << "	" << tot_meas(psi, sites, N) << "\n";
            }


        }
    return(0);
}

