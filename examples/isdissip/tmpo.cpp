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
    
    auto t = 0.;
    auto tend = 1.;
    
    auto t0 = 0.01, ttry = 0.02, dp = 0.001, dpp = 0.002 ; // tdid <= ttry; dp << 1; dpp >= dp;

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
			else				h = atof( &argv[1][1] )*J/2.;
			break;
		case 't':
			if(argv[1][1] == 'o')		t0 = atof( &argv[1][2] );
			else if(argv[1][1] == 'd')	tend = atof( &argv[1][2] );
			break;
		case 'p':
			if(argv[1][1] == 'P')		dpp = atof( &argv[1][2] );
			else				dp = atof( &argv[1][2] );
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


    
    //printfln("Maximum bond dimension of expH1 is %d",maxLinkDim(expH1));
    auto args = Args("Method=","Fit","Cutoff=",1E-12,"MaxDim=",2000,"IsHermitian=",false);
        //DensityMatrix		&		Fit
   auto Cop = op(sites, "Sz", 1);
   double dpm(0.);
   field<itensor::MPS> psi2(ensem);
   field<double> obstime(20*int(tend/dp)+1), obs(20*int(tend/dp)+1);

   const double TRY = ttry;
   int sweep(0);
   for(int tt=0; tt < ensem; ++tt)
   	psi2(tt) = psi1;

   double delta, deltat;

 //for(int sweep=0; sweep<int(tend/t0); ++sweep ) {

 for(int tt=0; tt < ensem; ++tt){
	ttry = TRY;
	double tdid, tnext(ttry);
	sweep = 0;

	t = 0.;

	if ( tt%(ensem/10) == 0) {
		std::cout << int( float(tt)/float(ensem)*100. ) << "%\n";
		out_file << std::flush;
	}

   while(t < tend){
	ttry = tnext;
	tdid = std::min(ttry, t0*(sweep + 1) - t);

	auto expH1 = toExpH(ampo,-(1-1_i)/2*tdid*Cplx_i);
	auto expH2 = toExpH(ampo,-(1+1_i)/2*tdid*Cplx_i);

	psi1 = psi2(tt);
        psi2(tt) = applyMPO(expH1,psi2(tt),args);
        psi2(tt).noPrime();
        psi2(tt) = applyMPO(expH2,psi2(tt),args);
        psi2(tt).noPrime().normalize();

	if (1.*w*tdid < dpp){
		t += tdid;
		tnext = std::min(tnext, dp/1./w);
		dpm = 1./tdid*dist(mt);

    		if(dpm < 1.*w){
    			psi2(tt).position(1);
     			auto newA = Cop*psi2(tt)(1);
			newA.noPrime();
			psi2(tt).set(1,newA);
			psi2(tt).noPrime().normalize();
    		}
		//obs(sweep) += measure(3, psi2(tt), sites);
		//obstime(sweep) += t;
		double tempor = measure(3, psi2(tt), sites);
		delta = tempor - obs(sweep);
		obs(sweep) = obs(sweep) + delta / float(tt + 1.);

		deltat = t - obstime(sweep);
		obstime(sweep) = obstime(sweep) + deltat / float(tt + 1.);
//std::cout << dpm << "	" << t << "	" << tempor << " " << obstime(sweep) << "	" << obs(sweep) << "\n";
		++sweep;
	}
	else {
		tnext = dp / 1./w;
		psi2(tt) = psi1;
	}
    }
 }

 for(int j=0; j < sweep; ++j)
  	out_file << obstime(j) << "		" << obs(j) << "\n";

  //###ADAPTIVE###


    return 0;
}




/*
  double tempd(abs(obs(0)-mean_obs));
  int inddx(0);
  for(int tt=0; tt < ensem; ++tt){
	if(tempd > abs(obs(tt)-mean_obs) ) {
		tempd = abs(obs(tt)-mean_obs);
		inddx = tt;
	}
  }
  for(int tt=0; tt < ensem; ++tt)
	psi2(tt) = psi2(inddx);
*/
