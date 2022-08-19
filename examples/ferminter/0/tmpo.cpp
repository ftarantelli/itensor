#include "tmpo.hpp"
#include <fstream>




int main(int argc, char *argv[])
    {

	double pbc(-1.), t0, tend, t;
	int N = 20;

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

	auto sites = Fermion(N, {"ConserveQNs=", false} );
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
	auto state = InitState(sites,"Emp");
	for( auto n : range1(N) ) if( n%2==0 ) state.set(n,"Occ");
	auto sweeps = Sweeps(5); //number of sweeps is 5
	sweeps.maxdim() = 10,20,100,100,200;
	sweeps.cutoff() = 1E-10;
	auto [energy,psi] = dmrg(H,randomMPS(state),sweeps,"Silent");
	
	std::cout << energy << "		" << tot_meas(psi, sites, N) << "	Energy\n";

    char output[80];
    std::sprintf(output, "couplingN%dk%.0fj%.0f.dat", N, kappa, jack);

    std::ofstream out_file(output, std::ios::out | std::ios::trunc);
    out_file.precision(16);
    
    
/*
    auto psi1 = MPS(psi);
    auto psi2 = psi1;
    auto energy = real(innerC(psi1,H,psi1));
    
	double dir(1.), sum(hi);//, t00(hi*tau), tmax(4*n_cycle*hi*tau + t00);
	int events = 2*n_cycle*int((tau*hf - tau*hi)*2./J/t0) + 1;
	int frac = 10, prnt = int(2.*abs(sigma)*20);

	if( prnt > events ) prnt = events;

	hi = std::abs(hi);
    //std::cout << sum << " " << sum*N/DeltaL*2./J << " " << hi << " " << DeltaL << " " << en1 << " " << energy00 << " " << energy0 << " " << M00 << "\n";
    //exit(8);
*/   
    /*
    // QUENCH
    h = h - hf;
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += h, "Sx",i;
        }
        ampo += h, "Sx",N;
    */
/*
    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-12;
*/
    /*
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = std::pow(1, -expcutoff);//1E-12;
    sweeps.niter() = 10;
    */
/*
//std::cout << "hghgh11g\n";
    for(int n = 1; n <= nsw; ++n)
        {
        if(n < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi1,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",true,
                                      "Quiet",true});
            }
        
        // TDVP sweep
        energy = tdvp(psi1,H,-t,sweeps,{"Truncate",true,
                                        "DoNormalize",true,
                                        "Quiet",true,
                                        "NumCenter",1});
        }

    //printfln("\nEnergy after imaginary time evolution = %.10f",energy);
    //printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );

    //println("-------------------------------------MPO W^I 2nd order---------------------------------------");
//std::cout << "hghghg\n";
    auto expH1 = toExpH(ampo,-(1-1_i)/2*t0*Cplx_i);
    auto expH2 = toExpH(ampo,-(1+1_i)/2*t0*Cplx_i);
    auto expHh = toExpH(ampoh,-t0/tau*J/2.*(1+1_i)/2*t0*Cplx_i);
    //printfln("Maximum bond dimension of expH1 is %d",maxLinkDim(expH1));
    auto args = Args("Method=","DensityMatrix","Cutoff=",1E-12,"MaxDim=",2000);
    for(int n = 1; n <= events; ++n)
        {
        sum += dir * t0 / tau * J/2.;
        
        for(int i = 1; i <= N-1; ++ i)
        		ampo += -dir * (t0 / tau)*J/2., "Sx",i;
    	   ampo += -dir * (t0 / tau)*J/2., "Sx",N;
    	   
	   //h = C.time / tau;
	   //if ( int( (C.time-t00)/tau/2./hi ) % 2 == 0 )
	   if( (n+1) % (events/2/n_cycle) == 0 ) dir = -dir;

	   //ampo += ampoh * dir * (t0 / tau);
	   expH1 = toExpH(ampo,-(1-1_i)/2*t0*Cplx_i);
	   expH2 = toExpH(ampo,-(1+1_i)/2*t0*Cplx_i);
   
        psi2 = applyMPO(expH1,psi2,args);
        psi2.noPrime();
        psi2 = applyMPO(expH2,psi2,args);
        psi2.noPrime().normalize();
        //if(n%int(std::real(prnt/t0)) == 0) {
            
            //printfln("\nMaximum bond dimension at time %.1f is %d ", n*t0, maxLinkDim(psi2));
            //printfln("Energy using overlap at time %.1f is %.10f", n*t0, real(innerC(psi2,H,psi2)) );
            
            std::cout << n*t0 << "\n";
            out_file << sum*N/DeltaL*2./J << "	" << std::real(tot_meas(psi2, sites, N)/M00) << "		" << n*t0 << "\n" << std::flush;
           // }


        }
*/
    return(0);
}

