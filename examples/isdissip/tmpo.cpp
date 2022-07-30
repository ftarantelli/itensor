#include "tmpo.hpp"
#include <fstream>
#include <armadillo>

using namespace arma;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> dist(0., 1.);
//double mrn; //= dist(mt);	//rand()/RAND_MAX;
//mrn = dist(mt);

int main(int argc, char *argv[])
    {
    
    // The system will be a 2x20 ladder
    int Nx = 6;
    int Ny = 1, ensem(100);
    
    auto t = 0.05;
    auto tend = 1.;
    
    auto t0 = 0.01;
    
    auto J = 4.;
    auto g = 1.*J/2., w = 1.;
    auto h = 0.5*J/2., hf = 0.*J/2., hi = 0.5*J/2.;
    
while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'N':
				Nx = atoi( &argv[1][1] );
			break;
		case 'e':
				ensem = atoi( &argv[1][1] );
			break;
		case 'g':
				g = atof( &argv[1][1] )*J/2.;
			break;
		case 'w':
				w = atof( &argv[1][1] );//*J/2.;
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

    int nsw = tend/t0;
    int N = Nx*Ny;
    
    char output[80];
    std::sprintf(output, "isdisN%d.dat", N);//"%s/fotqu%.0fs%.0fl%d.dat", folder, upsilon*10., abs(sigma), Lsize );
    std::ofstream out_file(output, std::ios::out | std::ios::trunc);
    out_file.precision(16);
    
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
        
    //auto H = toMPO(ampo);
	
    for(int i = 1; i <= N-1; ++ i)
	   {
        ampo += -hi, "Sx", i;
        }
    ampo += -hi, "Sx",N;
    
   //ampo += ampoh*h;
    auto H = toMPO(ampo);
    
    auto psi0 = randomMPS(sites);
    auto sweeps0 = Sweeps(5); //number of sweeps is 5
    sweeps0.maxdim() = 10,20,100,100,200;
    sweeps0.cutoff() = 1E-10;

    auto [energy0,psi] = dmrg(H,psi0,sweeps0,"Silent");
    //auto obsite = 3;
    
    auto psi1 = MPS(psi);
    //std::cout << measure(3, psi1, sites) << "\n";
    //exit(8);
    //auto psi2 = psi1;
    auto energy = real(innerC(psi1,H,psi1));
    

    // QUENCH
    h = hi - hf;
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += h, "Sx", i;
        }
        ampo += h, "Sx", N;
        ampo += -Cplx_i*w/2., "Id", 1; 

    H = toMPO(ampo);
	

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-12;
    /*
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = std::pow(1, -expcutoff);//1E-12;
    sweeps.niter() = 10;
    */

/*
    std::sprintf(output, "tempN%d.dat", N);//"%s/fotqu%.0fs%.0fl%d.dat", folder, upsilon*10., abs(sigma), Lsize );
    std::ofstream host(output, std::ios::out | std::ios::trunc);
    host.precision(16);

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
        energy = tdvp(psi1,H,-t0,sweeps,{"Truncate",true,
                                        "DoNormalize",true,
                                        "Quiet",true,
                                        "NumCenter",1});
	host << n*t0 << "		" << measure(3, psi1, sites) << "\n";
        }
*/
    

    //printfln("\nEnergy after imaginary time evolution = %.10f",energy);
    //printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );

    //println("-------------------------------------MPO W^I 2nd order---------------------------------------");

    auto expH1 = toExpH(ampo,-(1-1_i)/2*t0*Cplx_i);
    auto expH2 = toExpH(ampo,-(1+1_i)/2*t0*Cplx_i);
    
    //printfln("Maximum bond dimension of expH1 is %d",maxLinkDim(expH1));
    auto args = Args("Method=","Fit","Cutoff=",1E-12,"MaxDim=",2000,"IsHermitian=",false);
        //DensityMatrix		&		Fit
   auto Cop = op(sites, "Sz", 1);
   double dpm(0.), obs(0.);
   field<itensor::MPS> psi2(ensem);
   
   for(int tt=0; tt < ensem; ++tt)
   	psi2(tt) = psi1;

 for(int sweep=0; sweep<int(tend/t0); ++sweep ) {
   obs = 0.;
   if ( sweep%(int(tend/t0)/10) == 0) {
	std::cout << int(sweep/tend*t0*100.) << "%\n";
	out_file << std::flush;
   }

   for(int tt=0; tt < ensem; ++tt){
   	//if(sweep==0) psi2(tt) = psi1;
    /*
    	psi2.position(N);
    	auto psi2dag = dag(prime(psi2(N),"Site"));
    	dpm = t0*real(eltC(psi2dag*Cop*psi2(N)));
    */
	dpm = w*t0;
	//std::cout << dpm << "\n";

    	if(dist(mt) <= dpm){
    		psi2(tt).position(1);
     	auto newA = Cop*psi2(tt)(1);
		newA.noPrime();
		psi2(tt).set(1,newA);
		psi2(tt).noPrime().normalize();
		//std::cout << dpm << "kkkkkkk\n";
    	}
    	else {
    	   psi2(tt) = applyMPO(expH1,psi2(tt),args);
        psi2(tt).noPrime();
        psi2(tt) = applyMPO(expH2,psi2(tt),args);
        psi2(tt).noPrime().normalize();
    	}
    	//out_file << sweep*t0 << "		" << tot_meas(psi2, sites, N) << "\n";
	obs += measure(3, psi2(tt), sites);
    }
  out_file << sweep*t0 << "		" << obs/ensem << "\n";
  }

    return 0;
    }

