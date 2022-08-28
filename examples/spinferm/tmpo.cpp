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

	double pbc(0.), t0(0.01), tend(1.), t(0.05);
	int N = 8;
	
	const double yg(2.), ymu(1.);
	double kappa(1.), jack(1.);
	
    int ensem(100);
	
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
		case 'e':
				ensem = atoi( &argv[1][1] );
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
	char output[80];
     std::sprintf(output, "spinfermdissN%d.dat", N);
	std::ofstream out_file(output, std::ios::out | std::ios::trunc);
	out_file.precision(16);

	int nsw = int(t/t0);
	auto sites = Electron(N, {"ConserveQNs=", false} );

	auto mu = -2. + kappa * std::pow(N, -ymu) ;
	auto delta = 1.;
	auto g = jack * std::pow(N, -yg);
	
	auto ampo = AutoMPO(sites);
	for( auto n : range1(N-1) )
   	 {
   	 //ampo += -delta,"Cdag",n,"Cdag",n+1;
   	 //ampo += -delta,"C",n+1,"C",n;
   	 
   	 ampo += -1.,"Cdagup",n,"Cup",n+1;
   	 ampo += -1.,"Cdagup",n+1,"Cup",n;

   	 ampo += -1.,"Cdagdn",n,"Cdn",n+1;
   	 ampo += -1.,"Cdagdn",n+1,"Cdn",n;
   	 
   	 //ampo += -mu,"Cdag",n,"C",n;
   	 //ampo += -g,"Cdag",n,"C",n, "Cdag",n+1,"C",n+1;
   	 
   	// ampo += - Cplx_i / 2., "Cdag", n, "C", n;
   	 }
	/*
   	 ampo += -mu,"Cdag",N,"C",N;

   	 ampo += -pbc*delta,"Cdag",N,"Cdag",1;
   	 ampo += -pbc*delta,"C",1,"C",N;
	*/

   	 ampo += -pbc,"Cdagup",N,"Cup",1;
   	 ampo += -pbc,"Cdagup",1,"Cup",N;

   	 ampo += -pbc,"Cdagdn",N,"Cdn",1;
   	 ampo += -pbc,"Cdagdn",1,"Cdn",N;

   	 //ampo += - Cplx_i / 2., "Cdag", N, "C", N;

	auto H = toMPO(ampo);
	//auto state = InitState(sites,"Emp");
	//for( auto n : range1(N) ) if( n%2==0 ) state.set(n,"Occ");
	auto sweeps0 = Sweeps(5); //number of sweeps is 5
	sweeps0.maxdim() = 10,20,100,100,200;
	sweeps0.cutoff() = 1E-10;
	auto [energy0,psi] = dmrg(H,randomMPS(sites),sweeps0,"Silent");
	
	//std::cout << energy0 << "		" << tot_meas(psi, sites, N) << "	Energy\n";

    auto psi1 = MPS(psi);

    auto energy = real(innerC(psi1,H,psi1));

    int zq = int( (N+1)/2 );

Real cutoff = 1E-8;
auto gates = vector<BondGate>();

//Create the gates exp(-i*tstep/2*hterm)
//and add them to gates
for(int b = 1; b <= N-1; ++b)
    {
    //auto hterm =  - mu*op(sites,"N",b)*op(sites, "Id", b+1 );
    
    //hterm -= delta*op(sites, "Cdag", b) * op(sites,"Cdag", b+1 );
    //hterm -= delta*op(sites, "C", b+1) * op(sites,"C", b );
    
    auto hterm = op(sites, "Cdagup", b) * op(sites,"Cup", b+1 );
    hterm -= op(sites, "Cdagup", b+1) * op(sites,"Cup", b );

    hterm -= op(sites, "Cdagdn", b) * op(sites,"Cdn", b+1 );
    hterm -= op(sites, "Cdagdn", b+1) * op(sites,"Cdn", b );
    
    //hterm -= g*op(sites, "N", b) * op(sites,"N", b+1 );
    
    //hterm -= mu*op(sites,"N",b)*op(sites, "Id", b+1 );
    //std::cout << b << "  dms;vm\n";
    //hterm -= Cplx_i /2. * op(sites,"N",b)*op(sites, "Id", b+1 );

    if ( b==zq ) hterm -= Cplx_i /2. * op(sites,"Nupdn",b)*op(sites, "Id", b+1 );

    auto g = BondGate(sites,b,b+1,BondGate::tReal,t0/2.,hterm);
    //std::cout << b << "  dms;vm\n";
    gates.push_back(g);
    //std::cout << b << "  dms;vm\n";
    }

//Create the gates exp(-i*tstep/2*hterm) in reverse
//order (to get a second order Trotter breakup which
//does a time step of "tstep") and add them to gates
for(int b = N-1; b >= 1; --b)
    {
    //auto hterm = - mu*op(sites,"N",b+1)*op(sites, "Id", b );
    
    //hterm -= delta*op(sites, "Cdag", b) * op(sites,"Cdag", b+1 );
    //hterm -= delta*op(sites, "C", b+1) * op(sites,"C", b );
    
    auto hterm = op(sites, "Cdagup", b) * op(sites,"Cup", b+1 );
    hterm -= op(sites, "Cdagup", b+1) * op(sites,"Cup", b );

    hterm -= op(sites, "Cdagdn", b) * op(sites,"Cdn", b+1 );
    hterm -= op(sites, "Cdagdn", b+1) * op(sites,"Cdn", b );    

    //hterm -= g*op(sites, "N", b) * op(sites,"N", b+1 );
    
    if ( b==zq ) hterm -= Cplx_i /2. * op(sites,"Nupdn",b)*op(sites, "Id", b+1 );
    
    auto g = BondGate(sites,b,b+1,BondGate::tReal,t0/2.,hterm);
    gates.push_back(g);
    }

//Save initial state;	psi2

//Time evolve, overwriting psi when done
//gateTEvol(gates,tend,t0,psi2,{"Cutoff=",cutoff,"Verbose=",true});


 
    auto Cop1 = op(sites, "Cup", zq);
    auto Cop2 = op(sites, "Cdn", zq);

    double dpm(0.), obs(0.);
   field<itensor::MPS> psi2(ensem);
    
 for(int sweep=0; sweep<nsw; ++sweep ) {
   obs = 0.;

   if ( sweep % (nsw/10) == 0) {
		std::cout << int( sweep/float(nsw)*100. ) << "%\n";
		out_file << std::flush;
   }

   for(int tt=0; tt < ensem; ++tt){
   	if(sweep==0) psi2(tt) = psi1;
    /*
    	psi2.position(zq);
    	auto psi2dag = dag(prime(psi2(zq),"Site"));
    	dpm = t0*real(eltC(psi2dag*Cop*psi2(zq)));
	*/
	dpm = t0*measure(zq, psi2(tt), sites);
	//std::cout << dpm/t0 << "\n";
 
    	if(dist(mt) <= dpm){
    		psi2(tt).position(zq);
     		auto newA1 = Cop1*psi2(tt)(zq);
		newA1.noPrime();
		psi2(tt).set(zq,newA1);
     		auto newA2 = Cop2*psi2(tt)(zq);
		newA2.noPrime();
		psi2(tt).set(zq,newA2);
		psi2(tt).noPrime().normalize();
		//std::cout << dpm << "kkkkkkk\n";
    	}
    	else {
    		gateTEvol(gates,t0,t0,psi2(tt),{"Cutoff=",cutoff,"Silent=",true});
    	}
    	//out_file << sweep*t0 << "		" << tot_meas(psi2, sites, N) << "\n";
	obs += tot_meas(psi2(tt), sites, zq);
    }
  out_file << sweep*t0 << "		" << obs/ensem << "\n";

//std::cout << sweep*t0 << "		" << obs/ensem << "\n";
  }

    return(0);
}

