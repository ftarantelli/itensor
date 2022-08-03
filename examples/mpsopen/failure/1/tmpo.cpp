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
    
    auto t = 0.1;
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
			else				h = atof( &argv[1][1] )*J/2.;
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
    auto lampo = AutoMPO(sites);
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += -J,"Sx",i,"Sx",i+1;
        ampo += -g, "Sz",i;
        lampo += -J,"Sx",i,"Sx",i+1;
        lampo += -g, "Sz",i;
        }
        ampo += -g, "Sz",N;
        lampo += -g, "Sz",N;
        
    //auto H = toMPO(ampo);
	
    for(int i = 1; i <= N-1; ++ i)
	   {
        ampo += -hi, "Sx", i;
        lampo += -hi, "Sx", i;
        }
    ampo += -hi, "Sx",N;
    lampo += -hi, "Sx",N;
    
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
    auto psi2 = psi1;
    auto energy = real(innerC(psi1,H,psi1));
    

    // QUENCH
    h = hi - hf;
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += h, "Sx", i;
        lampo += h, "Sx", i;
        }
        ampo += h, "Sx", N;
        lampo += h, "Sx", N;

    //H = toMPO(ampo);

    //auto lampo1 = AutoMPO(sites);
    //auto lampo2 = AutoMPO(sites);
    //lampo1 += ampo;
    //lampo2 += ampo;	
    for(int i = 1; i <= N; ++ i)
        {
        ampo += Cplx_i, "Id", i; //? - or not
	lampo += -Cplx_i, "Id", i;
        }
        ampo += 2.*w, "Sz", 1;
	lampo += 2.*w, "Sz", 1;

	ampo += -w/2., "Id", 1;
	lampo += -w/2., "Id", 1;

	ampo += -w/2., "Id", 1;
	lampo += -w/2., "Id", 1;

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-12;


    auto expH1 = toExpH(ampo,-(1-1_i)/2*t0*Cplx_i);
    auto expH2 = toExpH(ampo,-(1+1_i)/2*t0*Cplx_i);

    auto exp2H1 = toExpH(lampo,-(1-1_i)/2*t0*Cplx_i);
    auto exp2H2 = toExpH(lampo,-(1+1_i)/2*t0*Cplx_i);
    
    auto args = Args("Method=","Fit","Cutoff=",1E-12,"MaxDim=",2000,"IsHermitian=",false);
        //DensityMatrix		&		Fit
   
   double obs(0.);

   auto mampo = AutoMPO(sites);
   mampo += 2., "Sz", 3;
   auto Szz = toMPO(mampo);

 for(int sweep=0; sweep<int(tend/t0); ++sweep ) {
   
   if ( sweep%(int(tend/t0)/10) == 0) {
	std::cout << int(sweep/tend*t0*100.) << "%\n";
	out_file << std::flush;
   }
           psi1 = applyMPO(expH1,psi1,args);
           psi1.noPrime();
           psi1 = applyMPO(expH2,psi1,args);
           psi1.noPrime().normalize();

           psi2 = applyMPO(exp2H1,psi2,args);
           psi2.noPrime();
           psi2 = applyMPO(exp2H2,psi2,args);
           psi2.noPrime().normalize();
    	
	//obs = measure(3, psi1, psi2, sites);
	obs = real(innerC(psi1, Szz, psi2));	//real(innerC(psi1, psi2));	//innerC(MPS y, MPO A, MPS x)
	
        out_file << sweep*t0 << "		" << obs << "\n";

 }


return 0;
}

