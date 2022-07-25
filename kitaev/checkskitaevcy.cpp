#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <math.h>

#include <armadillo>

using namespace arma;
typedef std::complex<double> complex;

const double YMU = 1., ZETA = 1., MUC = -2.;
int L(50);
double delta(1.), mu(MUC), tau;


class wire {
	public:
		//int L;
		dvec km;
		//cx_dmat gs;
		field<cx_dvec> gs;
		//field<sp_cx_dmat> hamk;
		field<cx_dmat> hamk, hamk0;
		field<cx_dmat> vpk;
		complex energy;
		//double Cval, Aval;
		//double time;
		
		wire(void) {
			//L = Lsize;
			km.set_size(int(L/2));
			hamk.set_size(int(L/2));
			hamk0.set_size(int(L/2));
			vpk.set_size(int(L/2));
			gs.set_size(int(L/2));
			for(int n=0; n<int(L/2); ++n) {
				km(n) = datum::pi/float(L)*(2.*float(n) + 1.);
				hamk(n).zeros(4,4);
				hamk0(n).zeros(4,4);
				vpk(n).zeros(4,4);
				gs(n).set_size(4);
			}
			energy = 0.;
			//time = (mu - MUC)*tau;
		}
		
		void hamkBuild(void);
		void GroundState(void);
		complex Correlations(int point1, int point2, char chr);
		inline void derivate(double t, const cx_dvec& _gsk, cx_dvec& kvec,int k);
		void RungeKuttaPsi(double _dt);
};


void wire::hamkBuild(void){

	for(int k=0; k<int(L/2); ++k) {
		vpk(k)(1,1) = - 1.;
		vpk(k)(2,2) = - 1.;
		vpk(k)(3,3) = - 2.;
		/*
		hamk(k)(0,3) = 2*delta*sin(km(k));
		hamk(k)(3,0) = 2*delta*sin(km(k));
		hamk(k)(1,1) = - (2.*cos(km(k)) + mu);
		hamk(k)(2,2) = - (2.*cos(km(k)) + mu);
		hamk(k)(3,3) = - 2.*(2.*cos(km(k)) + mu);
		*/
		hamk0(k)(0,3) = 2*delta*sin(km(k));
		hamk0(k)(3,0) = 2*delta*sin(km(k));
		hamk0(k)(1,1) = - (2.*cos(km(k)) );
		hamk0(k)(2,2) = - (2.*cos(km(k)) );
		hamk0(k)(3,3) = - 2.*(2.*cos(km(k)));
	}
//std::cout << mu << "a\n";
}

void wire::GroundState(void){

	dvec eigval;
	//cx_dvec eigval;
	cx_dmat eigvec;

	energy = 0.;
	for(int k=0; k<int(L/2); ++k) {
		/*
		//eigs_gen(eigval, eigvec, hamk(k), 1, "sr");
		eig_sym(eigval, eigvec, hamk(k));
		//std::cout << hamk(k) << '\n';
		for(int i=0; i<4; ++i)
			gs(k)(i) = (eigvec.col(0))(i);
		gs(k) = gs(k) / norm(gs(k));
		energy += eigval(0);
		//std::cout << km(k) << " " << mu << delta << " " << eigval(0) << "\n";
		//std::cout << gs(k) << "\n";	
		*/
		
		gs(k)(0) = -((-2.*1./tan(km(k)) - mu *  1./sin(km(k)) + 
  std::pow(2. + 2.*std::pow(delta,2.) + std::pow(mu,2.) + 4.* mu*cos(km(k)) + 2.*cos(2.*km(k)) - 2.*std::pow(delta,2.)*cos(2.*km(k)), 0.5) * 1./sin(km(k)) )/(2.*delta)) ;
  		gs(k)(1) = 0.; gs(k)(2) = 0.; gs(k)(3) = 1.;
  		gs(k) = gs(k) / norm(gs(k));
  		energy += -mu-2.*cos(km(k))-std::pow(2. +2.* std::pow(delta,2.)+std::pow(mu,2.)+4.* mu* cos(km(k))+2.* cos(2*km(k))-2.*std::pow(delta,2.)*cos(2.*km(k)), 0.5 );
  		
		//exit (8);
		
	}

}

complex wire::Correlations(int point1, int point2 = 0, char chr = 'C') {

	complex out(0.);
	if(point1 == 0) {
		// Particle number
		for(int k=0; k<int(L/2); ++k)
			out += ( conj(gs(k)(1))*gs(k)(1) + conj(gs(k)(2))*gs(k)(2) + 2.*conj(gs(k)(3))*gs(k)(3) );
//std::cout << "a0\n";
	}
	else{
		// c_x c_y	 --->	x-y
		if(chr == 'P'){
			//std::cout << "a1\n";
			for(int k=0; k<int(L/2); ++k)
				out += -2./float(L) * sin( km(k)*(point1-point2) ) * ( conj(gs(k)(3))*gs(k)(0) + conj(gs(k)(0))*gs(k)(3) );}
		else{	
			//std::cout << "a2\n";
			for(int k=0; k<int(L/2); ++k)
				out += 2./float(L) * cos( km(k)*(point1-point2) ) * ( conj(gs(k)(1))*gs(k)(1) + conj(gs(k)(2))*gs(k)(2) + 2.*conj(gs(k)(3))*gs(k)(3) );
		}
	
	}
	return(out);

}



inline void wire::derivate( double t, const cx_dvec& _gsk, cx_dvec& kvec, int k) {
	
		kvec = complex(0., -1.) * ( hamk(k) * _gsk );
		
}

void wire::RungeKuttaPsi(double _dt) {
	
	cx_dvec kvec1(4), kvec2(4), kvec3(4), kvec4(4), temp(4);
	double time(0.);
	for(int k=0; k<int(L/2); ++k) {
//std::cout << "aa\n";
		derivate( time         , gs(k)  , kvec1, k);
//	std::cout << "aa1\n";
		temp = gs(k) + kvec1*_dt/2.;
		derivate( time + _dt/2., temp, kvec2, k);
		temp = gs(k) + kvec2*_dt/2.;
		derivate( time + _dt/2., temp, kvec3, k);
		temp = gs(k) + kvec3*_dt;
		derivate( time + _dt   , temp, kvec4, k);

		gs(k) = gs(k) + _dt/6.*( kvec1 + 2.*kvec2 + 2.*kvec3 + kvec4 );
	
		gs(k) = gs(k) / norm(gs(k));
	}
}


int main(int argc, char *argv[]) {

	int events, k_events, point(3), n_cycle(1);
	double dt(0.02), mui(-1.), muf, dk(0.01);	//, tmax(100.), t00(0.)
	double upsilon(1.), kappa(-1), kappaf(-2.);

while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'L':
				L = atoi( &argv[1][1] );
			break;
		case 'w':
				mui = atof( &argv[1][1] );
			break;
		case 't':
				//if(argv[1][1] == 'M')	tmax = atof( &argv[1][2] );
				/*else*/					dt = atof( &argv[1][1] );
			break;
		case 'c':
				n_cycle = atoi( &argv[1][1] );
			break;
		case 'e':
				events = atoi( &argv[1][1] );
			break;
		case 'p':
				point = atoi( &argv[1][1] );
			break;
		case 'd':
				if(argv[1][1] == 'K')		dk = atof( &argv[1][2] );
				else						delta = atof( &argv[1][1] );
			break;
		case 'u':
				upsilon = atof( &argv[1][1] );
			break;
		case 'k':
				if(argv[1][1] == 'f')		kappaf = atof( &argv[1][2] );
				else						kappa = atof( &argv[1][1] );
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			std::cerr << "L - size \n m - mu\n t - delta time\n tM - time max\n c - num cycles\n e - num events\n p - point\n d - deltai\n u - upsilon\n s - sigma\n";
			exit (8);
	}
	++argv;
	--argc;
}

char output[80];
std::sprintf(output, "data/Swireu%.0fk%.0ff%.0fl%d.dat", 100*upsilon, abs(kappa), abs(kappaf), L);
//std::ofstream out_file(output, std::ios::out | std::ios::trunc);

//std::sprintf(output, "data/test.dat");
//std::ofstream out_file1(output, std::ios::out | std::ios::trunc);

//out_file.precision(16);
//out_file1.precision(16);
std::cout.precision(16);

std::sprintf(output, "# Sigma			C_{L/%d}			P_{L/%d}			A			rho-rho_c		C^o_{L/%d}		E - E_o", point, point, point) ;
//out_file << output << "\n";

tau = upsilon * std::pow(L, YMU + ZETA);

class wire K;
class wire E;
complex Cval, Aval;

mu = MUC;
E.hamkBuild();
E.GroundState();
complex rhoc = E.Correlations(0, 0, 'C');

K.hamkBuild();
E = K;

k_events = int( abs(kappaf - kappa)/dk );

for(int kitem=0; kitem <= k_events; ++kitem) {
	if ( kitem%(k_events/100) == 0) {
		//std::cout << int(inx/float(events)*100.) << "%\n";
		//out_file << std::flush;
	}

	mui = MUC + kappa * std::pow(tau, -1.*YMU / (YMU + ZETA));
	muf = - mui + 2.*MUC;

	mu = mui;
	events = 2*n_cycle*int( (tau*muf - tau*mui) / dt ) + 1;

	for(int j=0; j<int(L/2); ++j)
			K.hamk(j) = K.hamk0(j) + mu * K.vpk(j);	//***--->>>
	K.GroundState();

	double sign(1.);
	cx_dvec een(1);

	for(int sweep = 0; sweep < events; ++sweep) {
		//if ( sweep%(events/10) == 0) {
			//std::cout << int(sweep/float(events)*100.) << "%\n";
			//out_file << std::flush;
		//}
	
		//K.time += dt;	//***--->>>
		//K.time += tau*dt;
		for(int j=0; j<int(L/2); ++j)
			K.hamk(j) = K.hamk(j) + sign * dt / tau * K.vpk(j);	//***--->>>
		//K.hamk(j) = K.hamk(j) + sign * dt * K.vpk(j);
		mu += sign * dt / tau;		//***--->>>
		//E.hamk = K.hamk;
		//E.time = K.time;

		//if( sweep == events/2 - 1 ) sign = -sign;
		if( (sweep+1) % (events/2/n_cycle) == 0 ) sign = -sign;
	
		K.RungeKuttaPsi(dt);
		
		//if( mu + dt/tau > muf ) {
		if((mu+ dt/tau -MUC) * std::pow(tau, YMU/(YMU + ZETA)) >= abs(kappa)) {
			Aval = 1.;
			een(0) = 0.;
//std::cout << (mu-MUC) * std::pow(tau, YMU/(YMU + ZETA)) << "\n";		
			E.GroundState();
			Cval = K.Correlations(L/point,0,'C');
			
			for(int j=0; j<int(L/2); ++j) {
				Aval *= norm( E.gs(j).t()*K.gs(j) );
				een += K.gs(j).t() * (K.hamk(j) * K.gs(j)) ;
			}
			Aval = -std::log(Aval);
			
			std::cout << (mu-MUC) * std::pow(tau, YMU/(YMU + ZETA)) << "	" << real(Cval) << "	" << real(K.Correlations(L/point, 0, 'P')) << "	" << real(Aval) << "	" << real(K.Correlations(0, 0, 'C' ) - rhoc)/float(L) << "	" << real(E.Correlations(L/point, 0, 'C')) << "	" << real(een(0) - E.energy) << '\n';
			}
	
	}
Aval = 1.;
een(0) = 0.;
E.GroundState();
Cval = K.Correlations(L/point,0,'C');
			
for(int j=0; j<int(L/2); ++j) {
	Aval *= norm( E.gs(j).t()*K.gs(j) );
	een += K.gs(j).t() * (K.hamk(j) * K.gs(j)) ;
}
Aval = -std::log(Aval);
			
std::cout << (mu-MUC) * std::pow(tau, YMU/(YMU + ZETA)) << "	" << real(Cval) << "	" << real(K.Correlations(L/point, 0, 'P')) << "	" << real(Aval) << "	" << real(K.Correlations(0, 0, 'C' ) - rhoc)/float(L) << "	" << real(E.Correlations(L/point, 0, 'C')) << "	" << real(een(0) - E.energy) << '\n';
	
kappa -= dk;
exit(8);
}

return(0);
}
