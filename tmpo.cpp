#include "tmpo.hpp"
#include <fstream>




int main(int argc, char *argv[])
    {
    
    char output[80];
    std::sprintf(output, "data.dat");//"%s/fotqu%.0fs%.0fl%d.dat", folder, upsilon*10., abs(sigma), Lsize );
    std::ofstream out_file(output, std::ios::out | std::ios::trunc);
    out_file.precision(16);
    auto expcutoff(12.);
    
    int n_cycle(1.);
    // The system will be a 2x20 ladder
    int Nx = 8;
    int Ny = 1;
    
    auto t = 0.05;
    auto tend = 1.;
    
    auto t0 = 0.01;
    
    auto J = 4.;
    auto g = 1.*J/2.;
    auto h = 0.5*J/2., hf = 0.*J/2., hi = 0.*J/2.;
    
    auto sigma(-5.), upsilon(0.1), tau(0.);
    
while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'N':
				Nx = atoi( &argv[1][1] );
			break;
		case 'g':
				g = atof( &argv[1][1] )*J/2.;
			break;
		case 'h': 
			if(argv[1][1] == 'i')		hi = atof( &argv[1][2] )*J/2.;
			else if(argv[1][1] == 'f')	hf = atof( &argv[1][2] )*J/2.;
			else						h = atof( &argv[1][1] )*J/2.;
			break;
		case 't':
			if(argv[1][1] == 'o')		t0 = atof( &argv[1][2] );
			else if(argv[1][1] == 'd')	tend = atof( &argv[1][2] );
			else						t = atof( &argv[1][1] );			
			break;
		case 'c':
				expcutoff = -atof( &argv[1][1] );
			break;
		case 'u':
				upsilon = atof( &argv[1][1] );
			break;
		case 's':
				sigma = atof( &argv[1][1] );
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

    int nsw = tend/t;
    int N = Nx*Ny;
    
        // Make N spin 1/2's
    auto sites = SpinHalf(N, {"ConserveQNs=", false});

    // Make the Hamiltonian for rung-decoupled Heisenberg ladder
    auto ampo = AutoMPO(sites);    
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += -J,"Sx",i,"Sx",i+1;
        ampo += -g, "Sz",i;
        }
        ampo += -g, "Sz",N;
        
    auto H = toMPO(ampo);
    //printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    
// #################################################################################

    auto sweeps0 = Sweeps(5); //number of sweeps is 5
    sweeps0.maxdim() = 20,50,200,200,400;
    sweeps0.cutoff() = 1E-12;

    auto psi0 = randomMPS(sites);

    auto [energy00,psi00] = dmrg(H,psi0,sweeps0,{"Quiet=",true});
    //std::cout << energy0 << "\n";
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi00;

    auto [en1,psi11] = dmrg(H,wfs,randomMPS(sites),sweeps0,{"Quiet=",true,"Weight=",20.0});
// #################################################################################

    auto DeltaL = en1 - energy00;
    
    	hi = sigma * DeltaL / N *J/2.;
	hf = - sigma * DeltaL / N *J/2.;
	tau = upsilon / std::pow(DeltaL, 2.) * N;
	h = hi;
	
	auto ampoh = AutoMPO(sites);
	for(int i = 1; i <= N-1; ++ i)
	   {
        ampo += -h, "Sx", i;
        ampoh += -1., "Sx", i;
        }
     ampoh += -1., "Sx",N;
     ampo += -h, "Sx",N;
    
   //ampo += ampoh*h;
    H = toMPO(ampo);
    
    psi0 = randomMPS(sites);

    auto [energy0,psi] = dmrg(H,psi0,sweeps0,{"Quiet=",true});
    //std::cout << energy0 << "\n";
    auto M00 = tot_meas(psi, sites, N);

    //auto obsite = 3;
    
    auto psi1 = MPS(psi);
    auto psi2 = psi1;
    auto energy = real(innerC(psi1,H,psi1));
    
	double dir(1.), sum(hi);//, t00(hi*tau), tmax(4*n_cycle*hi*tau + t00);
	int events = 2*n_cycle*int((tau*hf - tau*hi)/t0) + 1;
	int frac = 10, prnt = int(2.*abs(sigma)*20);

	if( prnt > events ) prnt = events;

	hi = std::abs(hi);
    //std::cout << sum << " " << sum*N/DeltaL*2./J << " " << hi << " " << DeltaL << " " << en1 << " " << energy00 << " " << energy0 << " " << M00 << "\n";
    //exit(8);
    
    /*
    // QUENCH
    h = h - hf;
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += h, "Sx",i;
        }
        ampo += h, "Sx",N;
    */

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-12;
    /*
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = std::pow(1, -expcutoff);//1E-12;
    sweeps.niter() = 10;
    */

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
        if(n%int(std::real(prnt/t0)) == 0)
            {
            //printfln("\nMaximum bond dimension at time %.1f is %d ", n*t0, maxLinkDim(psi2));
            //printfln("Energy using overlap at time %.1f is %.10f", n*t0, real(innerC(psi2,H,psi2)) );
            
            std::cout << n*t0 << "\n";
            out_file << sum*N/DeltaL*2./J << "	" << std::real(tot_meas(psi2, sites, N)/M00) << "		" << n*t0 << "\n" << std::flush;
            }


        }

    return 0;
    }
