#include "qtraj.h"

int Lsize, m, m1, n_sweeps, n_qtraj;
double h1;
double dt, tmax;
double dp, dpp;

int main( int argc, char *argv[]) {

Lsize = 6; m = 5; m1 = m; n_sweeps = 1; dt = 0.01; tmax = 10.;
g = 1.; h = 0.5; h1 = 0.;
dp = 0.01; dpp = 0.02;
n_qtraj = 50;
w = 1.;

while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'L':
			Lsize = atoi( &argv[1][1] );
			break;
		case 'g':
			g = atof( &argv[1][1] );	//std::stod( &char &chr, int size_chr )
			break;
		case 'h':
			if(argv[1][1] == 'i')		h = atof( &argv[1][2] );
			else if(argv[1][1] == 'f')	h1 = atof( &argv[1][2] );
			else {	h = atof( &argv[1][1] );	h1 = h;  }
			break;
		case 'm':
			if(argv[1][1] == 'f')		m1 = atoi( &argv[1][2] );
			else				m = atoi( &argv[1][1] );
			break;
		case 'n':
			//if(argv[1][1] == 'i')		n_sweeps = atoi( &argv[1][2] );
			//else if(argv[1][1] == 't')	n_time = atoi( &argv[1][2] );
			//else				n_time = atoi( &argv[1][1] );
			n_sweeps = atoi( &argv[1][1] );
			break;
		case 't':
			if(argv[1][1] == 'M')		tmax = atof( &argv[1][2] );
			else				dt = atof( &argv[1][1] );
			break;
		case 'q':
			n_qtraj = atoi( &argv[1][1] );
			break;
		case 'w':
			w = atof( &argv[1][1] );
			break;
		case 'p':
			if(argv[1][1] == 'p')	dpp = atof( &argv[1][2] );
			else			dp  = atof( &argv[1][1] );
			break;
		default:
			std::cerr << "Unlucky\n";
			exit (8);
	}
	++argv;
	--argc;
}

class chain C(Lsize);

char output[30];
std::sprintf(output, "data/dissL%dg%.0lfh%.0f_%.0lf.dat", Lsize, g, h, 10.*h );
std::ofstream out_file(output, std::ios::out | std::ios::trunc);
out_file << "#	time	Mz	\n";

C.FiniteDMRG(n_sweeps, m, out_file);

// TIME DMRG
h = h1;
// m = m1;

/*
//	QUANTUM TRAJECTORIES ADP	---	modify BondExp
C.initSToper(dt, w);
C.InitialEvolEff();

int dimO = int(tmax/dt);
VectorXcd obs = VectorXcd::Zero(dimO), vaux(dimO);
chain Caux0(C);


for(int tr=0; tr < n_qtraj; ++tr ) {
	if ( tr%(int(n_qtraj/10) + 1) == 0) std::cout << 100*tr/n_qtraj << "%\n";
	
	C = Caux0;
	C.QTraj(m, dt,  tmax, vaux);

	for(int iv=0; iv < dimO; ++iv ) {
		obs(iv) += vaux[iv]; 
	}
}

obs = obs / complex(float(n_qtraj), 0.);
for(int j=0; j < dimO; ++j) {

	out_file << float(j+1)*dt << "		" << obs(j).real() << '\n';

}
*/


//	RUNGE-KUTTA	---	modify BondExp
C.InitialEvol();

ArrayX< MatrixXcd > Uleft(Lsize-1), Uright(Lsize-1);
MatrixXcd U1=C.BondExp(dt, 'D'),U2=C.BondExp(dt, 'O'),II=MatrixXcd::Identity(4,4);

for(int i = 1; i < Lsize-2; i++) {
	if(i%2 == 0) { // Lleft		k = Lleft - 1
		Uleft[i] = II;
		Uright[i] = U2;
	} else {
		Uleft[i] = U2;
		Uright[i] = II;
	}
}
Uleft[0] = U2; Uright[Lsize-2] = U1;
Uleft[Lsize-2] = II; Uright[0] = II;

while(C.time < tmax) {
	double frac=10;
	if( frac*C.time/tmax - int(frac*C.time/tmax) < frac*dt/tmax ) 
		std::cout << int(100.*(C.time/tmax)) << "%\n";

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft[i];

	C.SweepRK(LEFT2RIGHT, m, dt, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uright[i];

	C.SweepRK(RIGHT2LEFT, m, dt, true, true, out_file);
	
	//for(int i = 0; i < Lsize-1; i++) C.texp[i] = II;
	
	//C.SweepRK(RIGHT2LEFT, m, dt, false, true, out_file);

	C.time += dt;
}


/*
//	QUANTUM TRAJECTORIES	---	modify BondExp
C.initSToper(dt, w);
C.InitialEvolEff();

int dimO = int(tmax/dt);
VectorXcd obs = VectorXcd::Zero(dimO), vaux(dimO);
VectorXd vauxt(dimO), datat = VectorXd::Zero(dimO);
chain Caux0(C);


for(int tr=0; tr < n_qtraj; ++tr ) {
	if ( tr%(int(n_qtraj/10) + 1) == 0) std::cout << 100*tr/n_qtraj << "%\n";
	
	C = Caux0;
	C.QTraj(m, dt, tmax, vaux);
	//C.QTrajAdp(m, dt, dpp, dp, tmax, vaux, vauxt);	//dpp > dp;

	for(int iv=0; iv < dimO; ++iv ) {
		obs(iv) += vaux(iv);
		//datat(iv) += vauxt(iv);
	}
}

obs = obs / complex(float(n_qtraj), 0.);
//datat = datat / float(n_qtraj);
for(int j=0; j < dimO; ++j) {

	out_file << float(j+1)*dt << "		" << obs(j).real() << '\n';
	//out_file << datat(j) << "	" << obs(j).real() << '\n';

}
*/

/*


// SUZUKI-TROTTER 4-ORDER	---	modify BondExp

C.InitialEvol();

double th = 1./( 2.-std::cbrt(2.) );
ArrayX< MatrixXcd > Uleft1(Lsize-1), Uright1(Lsize-1);
ArrayX< MatrixXcd > Uleft2(Lsize-1), Uright2(Lsize-1);

MatrixXcd U1l=C.BondExp(dt*th/2., 'O'),		U2l=C.BondExp(dt*(1.-th)/2., 'O'),
	U1r=C.BondExp(dt*th, 'O'),		U1rd=C.BondExp(dt*th, 'D'),
	U2r=C.BondExp(dt*(1.-2.*th), 'O'),	U2rd=C.BondExp(dt*(1.-2.*th), 'D'),
	II=MatrixXcd::Identity(4,4);

for(int i = 1; i < Lsize-2; i++) {
	if(i%2 == 1) { // Lleft		k = Lleft - 1
		Uleft1[i] = II;
		Uright1[i] = U1r;
		Uleft2[i] = II;
		Uright2[i] = U2r;
	} else {
		Uleft1[i] = U1l;
		Uright1[i] = II;
		Uleft2[i] = U2l;
		Uright2[i] = II;
	}
}
Uleft1[0] = U1l; Uleft1[Lsize-2] = U1l; 
Uright1[0] = II; Uright1[Lsize-2] = U1rd;
Uleft2[0] = U2l; Uleft2[Lsize-2] = U2l; 
Uright2[0] = II; Uright2[Lsize-2] = U2rd;

while(C.time < tmax) {
	double frac=10;
	if( frac*C.time/tmax - int(frac*C.time/tmax) < frac*dt/tmax ) 
		std::cout << int(100.*(C.time/tmax)) << "%\n";

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft1[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uright1[i];

	C.Sweep(RIGHT2LEFT, m, false, true, out_file);
	
	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft2[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uright2[i];

	C.Sweep(RIGHT2LEFT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft2[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uright1[i];

	C.Sweep(RIGHT2LEFT, m, false, true, out_file);
	
	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft1[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = II;

	C.Sweep(RIGHT2LEFT, m, true, true, out_file);


	C.time += dt;
}

*/
/*	SUZUKI-TROTTER 2-ORDER	---	modify BondExp
C.InitialEvol();
ArrayX< MatrixXcd > Uleft(Lsize-1), Uright(Lsize-1);
MatrixXcd U1l=C.BondExp(dt/2., 'O'),U1r=C.BondExp(dt, 'O'),U2r=C.BondExp(dt, 'D'),II=MatrixXcd::Identity(4,4);

for(int i = 1; i < Lsize-2; i++) {
	if(i%2 == 1) { // Lleft		k = Lleft - 1
		Uleft[i] = II;
		Uright[i] = U1r;
	} else {
		Uleft[i] = U1l;
		Uright[i] = II;
	}
}
Uleft[0] = U1l; Uright[Lsize-2] = U2r;
Uleft[Lsize-2] = U1l; Uright[0] = II;

while(C.time < tmax) {
	double frac=10;
	if( frac*C.time/tmax - int(frac*C.time/tmax) < frac*dt/tmax ) 
		std::cout << int(100.*(C.time/tmax)) << "%\n";

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uright[i];

	C.Sweep(RIGHT2LEFT, m, false, true, out_file);
	
	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = II;

	C.Sweep(RIGHT2LEFT, m, true, true, out_file);


	C.time += dt;
}
*/

/*	SUZUKI-TROTTER 1-ORDER	---	modify BondExp
C.InitialEvol();
ArrayX< MatrixXcd > Uleft(Lsize-1), Uright(Lsize-1);
MatrixXcd U1=C.BondExp(dt, 'D'),U2=C.BondExp(dt, 'O'),II=MatrixXcd::Identity(4,4);

for(int i = 1; i < Lsize-2; i++) {
	if(i%2 == 0) { // Lleft		k = Lleft - 1
		Uleft[i] = II;
		Uright[i] = U2;
	} else {
		Uleft[i] = U2;
		Uright[i] = II;
	}
}
Uleft[0] = U2; Uright[Lsize-2] = U1;
Uleft[Lsize-2] = II; Uright[0] = II;

while(C.time < tmax) {
	double frac=10;
	if( frac*C.time/tmax - int(frac*C.time/tmax) < frac*dt/tmax ) 
		std::cout << int(100.*(C.time/tmax)) << "%\n";

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uleft[i];

	C.Sweep(LEFT2RIGHT, m, false, true, out_file);

	for(int i = 0; i < Lsize-1; i++) C.texp[i] = Uright[i];

	C.Sweep(RIGHT2LEFT, m, true, true, out_file);

	C.time += dt;
}

*/
//delete C; C = NULL;
return(0);
}
