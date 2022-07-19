#include "tdmrg.h"
#include <time.h>
#include <random>

double w;
double tmt = time(NULL);

bool do_jump(false);

// DISSIPATION WITH L_1 = \sigma^z at site = 1;

/*
class traj {
	public :
		int size;
		chain& CC;
		VectorXcd Obs;
		int tstep, tmaxm;
		
		traj(int _L) {
			size = _L;
			Obs.resize(int(tmaxm/tstep));
			chain temp(_L);
			CC = temp;
			
		}
};
*/




void chain::initSToper(double _dt, double _w) {
	int _Lsize = Lsize;
	//double th = 1./( 2.-std::cbrt(2.) );
	
	w = _w;
	
	UST.Uleft1.resize(_Lsize-1);
//	UST.Uleft2.resize(_Lsize-1);
	UST.Uright1.resize(_Lsize-1);
//	UST.Uright2.resize(_Lsize-1);

/*	
	UST.dt_ldiss[0] = _dt*th/2.; UST.dt_ldiss[1] = _dt*(1.-th)/2.;

MatrixXcd U1l=BondExp(_dt*th/2., 'O'),		U2l=BondExp(_dt*(1.-th)/2., 'O'),
	U1r=BondExp(_dt*th, 'O'),		U1rd=BondExp(_dt*th, 'D'),
	U2r=BondExp(_dt*(1.-2.*th), 'O'),	U2rd=BondExp(_dt*(1.-2.*th), 'D'),
	II=MatrixXcd::Identity(4,4);

	for(int i = 1; i < _Lsize-2; i++) {
		if(i%2 == 1) { // Lleft		k = Lleft - 1
			UST.Uleft1[i] = II;
			UST.Uright1[i] = U1r;
			UST.Uleft2[i] = II;
			UST.Uright2[i] = U2r;
		} else {
			UST.Uleft1[i] = U1l;
			UST.Uright1[i] = II;
			UST.Uleft2[i] = U2l;
			UST.Uright2[i] = II;
		}
	}
	UST.Uleft1[0] = U1l-complex(0., 1.)*_w/2.*II; UST.Uleft1[_Lsize-2] = U1l; 
	UST.Uright1[0] = II; UST.Uright1[_Lsize-2] = U1rd;
	UST.Uleft2[0] = U2l; UST.Uleft2[_Lsize-2] = U2l; 
	UST.Uright2[0] = II; UST.Uright2[_Lsize-2] = U2rd;
*/
	UST.dt_ldiss[0] = _dt; UST.dt_ldiss[1] = 0.;

	MatrixXcd U1=BondExp(_dt, 'D'), U2=BondExp(_dt, 'O'), II=MatrixXcd::Identity(4,4);

	for(int i = 1; i < _Lsize-2; i++) {
		if(i%2 == 0) { // Lleft		k = Lleft - 1
			UST.Uleft1[i] = II;
			UST.Uright1[i] = U2;
		} else {
			UST.Uleft1[i] = U2;
			UST.Uright1[i] = II;
		}
	}
	UST.Uleft1[0] = U2- complex(0., 1.)*w/2.*II; UST.Uright1[_Lsize-2] = U1;
	UST.Uleft1[_Lsize-2] = II; UST.Uright1[0] = II;
}

void chain::InitialEvolEff(void) {

	blockL[0].H = -g*Sz0 - h*Sx0-complex(0., 1.)*w/2.*MatrixXcd::Identity(2,2);
	blockR[0].H = -g*Sz0 - h*Sx0;
	
	srand(tmt);
}

void chain::BuildBlockEff(int direction, int iter) {
	MatrixXcd I2 = MatrixXcd::Identity(2,2);

	if(direction == LEFT) {
		Lleft = iter;
		blockL[iter].resize(iter+1);

		MatrixXcd HL = blockL[iter-1].H, SxL = blockL[iter-1].Sx[iter-1];
		//int diml = HL.rows();
		//MatrixXcd Ileft = MatrixXcd::Identity(diml, diml);
		MatrixXcd Ileft = blockL[iter-1].Id;

		blockL[iter].H = kron_prod(HL,I2) - g* kron_prod(Ileft,Sz0) - h* kron_prod(Ileft,Sx0) - kron_prod(SxL,Sx0); //- complex(0., 1.)*w/2.*kron_prod(Ileft, I2);
		blockL[iter].Id = kron_prod(Ileft, I2);
		for(int i=0; i<iter; ++i) {
			blockL[iter].Sx[i] = kron_prod(blockL[iter-1].Sx[i], I2);
			blockL[iter].Sz[i] = kron_prod(blockL[iter-1].Sz[i], I2);
		}
		blockL[iter].Sz[iter] = kron_prod(Ileft,Sz0);
		blockL[iter].Sx[iter] = kron_prod(Ileft,Sx0);
	} else {
	
		Rright = iter;
		blockR[iter].resize(iter+1);
		MatrixXcd HR = blockR[iter-1].H, SxR = blockR[iter-1].Sx[iter-1];
		//int dimr = HR.rows();
		//MatrixXcd Iright = MatrixXcd::Identity(dimr, dimr);
		MatrixXcd Iright = blockR[iter-1].Id;

		blockR[iter].H = kron_prod(I2, HR) - g* kron_prod(Sz0, Iright) - h* kron_prod(Sx0, Iright) - kron_prod(Sx0, SxR);// - complex(0., 1.)*w/2.*kron_prod(I2, Iright);
		blockR[iter].Id = kron_prod(I2, Iright);
		for(int i=0; i<iter; ++i) {
			blockR[iter].Sx[i] = kron_prod(I2, blockR[iter-1].Sx[i]);
			blockR[iter].Sz[i] = kron_prod(I2, blockR[iter-1].Sz[i]);
		}
		blockR[iter].Sz[iter] = kron_prod(Sz0, Iright);
		blockR[iter].Sx[iter] = kron_prod(Sx0, Iright);
	
	}
}




void chain::SweepEff(int verse, int _m, complex& _Mz, bool measure, bool evol){
	//complex Mz=complex(0.,0.);
	if(verse == LEFT2RIGHT){
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlockEff(LEFT, iter);	
			BuildBlockEff(RIGHT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				SuzukiTrotter();
			} else {
				GroundState();
				}
//			if(measure && Lleft == 1) {
//				Mz=obser();
//				outf << time << "	" << Mz.real() << '\n';
//			}

			BuildU(LEFT, _m);
		}
	} else if(verse == RIGHT2LEFT) {
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlockEff(RIGHT, iter);	
			BuildBlockEff(LEFT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				SuzukiTrotter();
			} else {
				GroundState();
				}
			if(measure && Lleft == 1) {
				_Mz = obser();
				//std::cout << _Mz << "\n";
				
			}

			BuildU(RIGHT, _m);
		}
	}
}

void chain::Jump(void) {
	int dim_r = blockR[Rright-1].Id.rows();
	assert(Lleft == 1);
	assert(gs.size() == dim_r*8);
	//MatrixXcd I_right = MatrixXcd::Identity(dim_r, dim_r);
	MatrixXcd I_right=kron_prod(MatrixXcd::Identity(2,2),blockR[Rright-1].Id);
	MatrixXcd Uev(dim_r*8,dim_r*8), SzI2(4,4);
	VectorXcd aux(gs);
	SzI2 = kron_prod(Sz0, MatrixXcd::Identity(2,2));
	Uev = kron_prod(std::sqrt(w)*SzI2, I_right);
	gs = Uev * aux;
	gs.normalize();
//std::cout << mrn  << " a0	QQ\n";
	return;
}



void chain::QTraj(int _m, double _dt, double tmaxm, VectorXcd& Obs) {

	//InitialEvolEff();
	assert(Obs.size() == int(tmaxm/_dt));
	int count_t(0);
	complex temp, temp1 = complex(1., 0.);
	MatrixXcd II = MatrixXcd::Identity(4,4);
	//double _dt(0);
	// at each cycle, we have to return to the initial conditions;
	
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(0., 1.);

	//while(time < tmaxm) {
	while( count_t < int(tmaxm/_dt) ){
		double frac = 3;
		if( frac*time/tmaxm - int(frac*time/tmaxm) < frac*_dt/tmaxm ) 
		//std::cout << "	" << int(100.*(time/tmaxm)) << "% time\n";
		std::cout << " -\n";
		
		double mrn; //= dist(mt);	//rand()/RAND_MAX;
		mrn = dist(mt); // rand() / float(RAND_MAX);
		//if( mrn < _dt && verse == LEFT2RIGHT ) {
//std::cout << mrn << "\n";			
		if( mrn < w*_dt) {
			Jump();
		} else {
		
			for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft1[i];

			SweepEff(LEFT2RIGHT, _m, temp, false, true);
		
			for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uright1[i];
		
			SweepEff(RIGHT2LEFT, _m, temp1, true, true);
		}
		time += _dt;
		//std::cout << temp1 << "\n";
		//Mz += temp1;
		Obs[count_t] = temp1;
		count_t += 1;
//std::cout << count_t-1 << " " << int(time/_dt) << " " << int(tmaxm/_dt) <<" \n";

	}
	//Mz = Mz / complex(float(_n_qtraj), 0.);
	std::cout << '\n' ;
	//outf << time << "	" << Mz.real() << '\n';
}

//	QTrajectoriesAdaptive

void chain::QTrajAdp(int _m, double _dt, double dpp, double dp, double tmaxm, VectorXcd& Obs, VectorXd& dataT){

	//InitialEvolEff();
	assert(dp < dpp);
	assert(dataT.size() == int(tmaxm/_dt));
	assert(Obs.size() == int(tmaxm/_dt));
	int count_t(0);
	complex temp, temp1 = complex(1., 0.);
	MatrixXcd II = MatrixXcd::Identity(4,4);
	// double _dt(0);
	// at each cycle, we have to return to the initial conditions;
	double dt_try(_dt);
	std::random_device rd;
	std::mt19937 mt(rd());
	// while(time < tmaxm) {
	while( count_t < int(tmaxm/_dt) ){
		double frac = 3;
		if( frac*time/tmaxm - int(frac*time/tmaxm) < frac*_dt/tmaxm ) 
		//std::cout << "	" << int(100.*(time/tmaxm)) << "% time\n";
		std::cout << " -\n";
		//do_jump = false;

		double dt_did = std::min(_dt, dt_try);
		
		if(dpp <= w*dt_did) {
		
			dt_try = dp / w ; // dpp > dp
		
		} else {
		
			initSToper(dt_did, w);
		
			for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft1[i];

			SweepEff(LEFT2RIGHT, _m, temp, false, true);
		
			for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uright1[i];
		
			SweepEff(RIGHT2LEFT, _m, temp1, false, true);

			dt_try = std::min(dp / w, dt_try);
			std::uniform_real_distribution<double> dist(0., 1./dt_did);

			double mrn; //= dist(mt);	//rand()/RAND_MAX;
			mrn = dist(mt); // rand() / float(RAND_MAX);

			if(mrn < w)
				Jump();
			time += dt_did;

			assert(Lleft == 1);
			temp1 = obser();
			//std::cout << temp1 << "\n";
			//Mz += temp1;
			Obs[count_t] = temp1;
			dataT[count_t] = time;
			count_t += 1;
		}
//std::cout << count_t-1 << " " << int(time/_dt) << " " << int(tmaxm/_dt) <<" \n";

	}
	//Mz = Mz / complex(float(_n_qtraj), 0.);
	std::cout << '\n' ;
	//outf << time << "	" << Mz.real() << '\n';
	
}

//	RUNGE-KUTTA

void chain::derivate( double t, MatrixXcd& rho, MatrixXcd& k) {
	int dim = rho.cols();
	k = MatrixXcd::Zero(dim, dim);
		if(Lleft == 1) {
			int dim_r = blockR[Rright].Id.rows();
			//MatrixXcd I_right = MatrixXcd::Identity(dim_r, dim_r);
			MatrixXcd I_right = blockR[Rright].Id;
			MatrixXcd Ham(dim_r*4,dim_r*4);
			Ham = kron_prod(texp[0], I_right);
			assert(dim == Ham.rows());
			k += complex(0., 1.) * ( rho*Ham - Ham*rho );
			MatrixXcd L1 = kron_prod(Sz0, MatrixXcd::Identity(2,2));
			L1 = kron_prod(L1, I_right);
			assert(L1.cols() == dim);
			k += w*(L1*(rho*L1.adjoint()) - rho);//0.5*(L1.adjoint()*L1)*rho - 0.5 * rho * (L1.adjoint()*L1) );
		}
	//} else {
		if(Rright == 1) {
			int dim_l = blockL[Lleft].Id.rows();
			//MatrixXcd I_left = MatrixXcd::Identity(dim_l, dim_l);
			MatrixXcd I_left = blockL[Lleft].Id;
			MatrixXcd Ham(dim_l*4,dim_l*4);
			Ham = kron_prod(I_left, texp[Lsize-2]);
			assert(dim == Ham.rows());
			k += complex(0., 1.) * ( rho*Ham - Ham*rho );
		}
	
	//int diml=blockL[Lleft].Id.rows()/2, dimr=blockR[Rright].Id.rows()/2;

	//MatrixXcd I_left = MatrixXcd::Identity(diml, diml);
	//MatrixXcd I_right = MatrixXcd::Identity(dimr, dimr);
	MatrixXcd I_left = blockL[Lleft-1].Id;
	MatrixXcd I_right = blockR[Rright-1].Id;
	MatrixXcd Ham;

	Ham = kron_prod(texp[Lleft], I_right);
	Ham = kron_prod(I_left, Ham);
	
	assert( rho.cols() == Ham.rows() );

	k += complex(0., 1.) * ( rho*Ham - Ham*rho );

}

void chain::RungeKutta(double _dt) {

	//int dim1=blockL[Lleft].H.rows(), dim2=blockR[Rright].H.rows();
	int dim;

	MatrixXcd rho;
	rho = gs * gs.adjoint();
	dim = rho.cols();
	MatrixXcd k1(dim,dim), k2(dim,dim), k3(dim,dim), k4(dim,dim), temp;

	derivate( time         , rho            , k1);
	temp = rho + k1*_dt/2.;
	derivate( time + _dt/2., temp, k2);
	temp = rho + k2*_dt/2.;
	derivate( time + _dt/2., temp, k3);
	temp = rho + k3*_dt;
	derivate( time + _dt   , temp, k4);

	rho = rho + _dt/6.*( k1 + 2.*k2 + 2.*k3 + k4 );
	//time = time + _dt;
	for(int i=0; i<dim; ++i) {
		gs(i) = rho(i, 0);
	}
	gs.normalize();
}


void chain::Blockderivate(double t, const MatrixXcd& rho, MatrixXcd& k, const MatrixXcd& Ham, const MatrixXcd& L1) {

	int dim = rho.cols();
	k = MatrixXcd::Zero(dim, dim);

	k += complex(0., 1.) * ( rho*Ham - Ham*rho );
	k += w*( L1*(rho*L1.adjoint()) - 0.5*(L1.adjoint()*L1)*rho - 0.5 * rho * (L1.adjoint()*L1) );;

}

void chain::BlockRungeKutta(double _dt) {

	double tstep = _dt / 2. / float(Lsize -3 );
	int dim = gs.size();
	
	MatrixXcd HL = blockL[Lleft].H, HR = blockR[Rright].H, Ham, rho;
	MatrixXcd SxL = blockL[Lleft].Sx[Lleft], SxR = blockR[Rright].Sz[Rright];
	MatrixXcd IL = blockL[Lleft].Id, IR = blockR[Rright].Id, L1;

	Ham = kron_prod(HL, IR) + kron_prod(IL, HR) - kron_prod(SxL, SxR);
	L1 = kron_prod(blockL[Lleft].Sz[0], IR);
	rho = gs * gs.adjoint();
	assert(dim == Ham.cols());
	assert(L1.cols() == dim);
	
	MatrixXcd k1(dim,dim), k2(dim,dim), k3(dim,dim), k4(dim,dim), temp;

	Blockderivate( time           , rho , k1, Ham, L1);
	temp = rho + k1*tstep/2.;
	Blockderivate( time + tstep/2., temp, k2, Ham, L1);
	temp = rho + k2*tstep/2.;
	Blockderivate( time + tstep/2., temp, k3, Ham, L1);
	temp = rho + k3*tstep;
	Blockderivate( time + tstep   , temp, k4, Ham, L1);

	rho = rho + tstep/6.*( k1 + 2.*k2 + 2.*k3 + k4 );
	//time = time + tstep;
	for(int i=0; i<dim; ++i) {
		gs(i) = rho(i, 0);
	}
	gs.normalize();

}


void chain::SweepRK(int verse, int _m, double _dt, bool measure, bool evol, std::ofstream& outf){
	complex Mz;
	if(verse == LEFT2RIGHT){
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlock(LEFT, iter);	
			BuildBlock(RIGHT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				RungeKutta(_dt);
				//BlockRungeKutta(_dt);
			} else {
				GroundState();
				}
//			if(measure && Lleft == 1) {
//				Mz=obser();
//				outf << time << "	" << Mz.real() << '\n';
//			}

			BuildU(LEFT, _m);
		}
	} else if(verse == RIGHT2LEFT) {
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlock(RIGHT, iter);	
			BuildBlock(LEFT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				RungeKutta(_dt);
				//BlockRungeKutta(_dt);
			} else {
				GroundState();
				}
			if(measure && Lleft == 1) {
				Mz=obser();
				if ((Mz.real() > 1.) || (Mz.real() < 0.)){
					std::cerr << "Overflow Error\n";}//exit (8);}

				outf << time + _dt << "	" << Mz.real() << '\n';
				//std::cout << Mz << '\n';
			}

			BuildU(RIGHT, _m);
		}
	}
}

//	QTrajectoriesEff
/*
//void chain::QTraj(int _m, int _n_qtraj, int _dt,  std::ofstream& outf) {
void chain::QTraj(int _m, double _dt, double tmaxm, VectorXcd& Obs) {

	//InitialEvolEff();
	assert(Obs.size() == int(tmaxm/_dt));
	int count_t(0);
	complex temp, temp1 = complex(1., 0.);
	MatrixXcd II = MatrixXcd::Identity(4,4);
	//double _dt(0);
	// at each cycle, we have to return to the initial conditions;

	//while(time < tmaxm) {
	while( count_t < int(tmaxm/_dt) ){
		double frac = 5;
		if( frac*time/tmaxm - int(frac*time/tmaxm) < frac*_dt/tmaxm ) 
		//std::cout << "	" << int(100.*(time/tmaxm)) << "% time\n";
		std::cout << " -\n";
//std::cout << "b\n"; 
		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft1[i];
//std::cout << UST.Uleft1[Lsize-2] << " " <<"a00\n"; exit (8);
		SweepEff(LEFT2RIGHT, _m, UST.dt_ldiss[0], temp, false, true);
//std::cout << "a01\n";
		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uright1[i];
//std::cout << "a02\n";
		SweepEff(RIGHT2LEFT, _m, _dt, temp, false, true);
//std::cout << "a03\n";
		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft2[i];

		SweepEff(LEFT2RIGHT, _m, UST.dt_ldiss[1], temp, false, true);

		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uright2[i];

		SweepEff(RIGHT2LEFT, _m, _dt, temp, false, true);

		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft2[i];

		SweepEff(LEFT2RIGHT, _m, UST.dt_ldiss[1], temp, false, true);

		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uright1[i];

		SweepEff(RIGHT2LEFT, _m, _dt, temp, false, true);
	
		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft1[i];

		SweepEff(LEFT2RIGHT, _m, UST.dt_ldiss[0], temp, false, true);

		for(int i = 0; i < Lsize-1; i++) texp[i] = II;

		SweepEff(RIGHT2LEFT, _m, _dt, temp1, true, true);
//std::cout << "a3\n";
		time += _dt;
		//std::cout << temp1 << "\n";
		//Mz += temp1;
		Obs[count_t] = temp1;
		count_t += 1;
//std::cout << count_t-1 << " " << int(time/_dt) << " " << int(tmaxm/_dt) <<" \n";

	}
	//Mz = Mz / complex(float(_n_qtraj), 0.);
	std::cout << '\n' ;
	//outf << time << "	" << Mz.real() << '\n';
}
*/


/*
void chain::SuzukiTrotterEff(int verse, double _dt) {
	
	//if(verse == LEFT2RIGHT) {
		if(Lleft == 1) {
			std::random_device rd;
			std::mt19937 mt(rd());
			std::uniform_real_distribution<double> dist(0., 1.);
			
			double mrn; //= dist(mt);	//rand()/RAND_MAX;
			mrn = dist(mt); // rand() / float(RAND_MAX);
			//if( mrn < _dt && verse == LEFT2RIGHT ) {
//std::cout << mrn << "\n";			
			if( mrn < w*_dt && verse == LEFT2RIGHT ) {
				int dim_r = blockR[Rright].Id.rows();
				assert(gs.size() == dim_r*4);
			//MatrixXcd I_right = MatrixXcd::Identity(dim_r, dim_r);
				MatrixXcd I_right = blockR[Rright].Id;
				MatrixXcd Uev(dim_r*4,dim_r*4), SzI2(4,4);
				VectorXcd aux(gs);
				SzI2 = kron_prod(Sz0, MatrixXcd::Identity(2,2));
				Uev = kron_prod(std::sqrt(w)*SzI2, I_right);
				gs = Uev * aux;
				gs.normalize();
				do_jump = true;
//std::cout << mrn  << " a0	QQ\n";
				return;
			} else {
// the Heff terms are inserted in the LEFT_OPERATOR of Suzuki-Trotter algorithm
			///////////////////////////////////////////////////////////
				int dim_r = blockR[Rright].Id.rows();
				assert(gs.size() == dim_r*4);
			//MatrixXcd I_right = MatrixXcd::Identity(dim_r, dim_r);
				MatrixXcd I_right = blockR[Rright].Id;
				MatrixXcd Uev(dim_r*4,dim_r*4);
				VectorXcd aux(gs);
				Uev = kron_prod(texp[0], I_right);
				gs = Uev * aux;
				gs.normalize();
//std::cout << mrn << "  b1\n";
			}
		}
	//} else {
		if(Rright == 1) {
			int dim_l = blockL[Lleft].Id.rows();
			assert(gs.size() == dim_l*4);
			//MatrixXcd I_left = MatrixXcd::Identity(dim_l, dim_l);
			MatrixXcd I_left = blockL[Lleft].Id;
			MatrixXcd Uev(dim_l*4,dim_l*4);
			VectorXcd aux(gs);
			Uev = kron_prod(I_left, texp[Lsize-2]);
			gs = Uev * aux;
			gs.normalize();
		}
	//}
	
	int diml=blockL[Lleft].Id.rows()/2, dimr=blockR[Rright].Id.rows()/2;
	VectorXcd vaux(gs);
	MatrixXcd temp;
	assert(gs.size()==diml*4*dimr);

	//MatrixXcd I_left = MatrixXcd::Identity(diml, diml);
	//MatrixXcd I_right = MatrixXcd::Identity(dimr, dimr);
	MatrixXcd I_left = blockL[Lleft-1].Id;
	MatrixXcd I_right = blockR[Rright-1].Id;
	assert( I_left.rows() * I_right.rows() * 4 == gs.size() );

	temp = kron_prod(texp[Lleft], I_right);
	temp = kron_prod(I_left, temp);
	
	gs = temp * vaux;
	
	gs.normalize();
	
}

void chain::SweepEff(int verse, int _m, double _dt, complex& _Mz, bool measure, bool evol){
	//complex Mz=complex(0.,0.);
	if(verse == LEFT2RIGHT){
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlockEff(LEFT, iter);	
			BuildBlockEff(RIGHT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				SuzukiTrotterEff(verse, _dt);
				if(do_jump == true) {
					for(int i = 0; i < Lsize-1; i++)
						texp[i] = MatrixXcd::Identity(4,4);
				}
			} else {
				GroundState();
				}
//			if(measure && Lleft == 1) {
//				Mz=obser();
//				outf << time << "	" << Mz.real() << '\n';
//			}

			BuildU(LEFT, _m);
		}
	} else if(verse == RIGHT2LEFT) {
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlockEff(RIGHT, iter);	
			BuildBlockEff(LEFT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				SuzukiTrotterEff(verse, _dt);
			} else {
				GroundState();
				}
			if(measure && Lleft == 1) {
				_Mz = obser();
				//std::cout << _Mz << "\n";
				
			}

			BuildU(RIGHT, _m);
		}
	}
}


void chain::QTraj(int _m, double _dt, double tmaxm, VectorXcd& Obs) {

	//InitialEvolEff();
	assert(Obs.size() == int(tmaxm/_dt));
	int count_t(0);
	complex temp, temp1 = complex(1., 0.);
	MatrixXcd II = MatrixXcd::Identity(4,4);
	//double _dt(0);
	// at each cycle, we have to return to the initial conditions;

	//while(time < tmaxm) {
	while( count_t < int(tmaxm/_dt) ){
		double frac = 3;
		if( frac*time/tmaxm - int(frac*time/tmaxm) < frac*_dt/tmaxm ) 
		//std::cout << "	" << int(100.*(time/tmaxm)) << "% time\n";
		std::cout << " -\n";
		do_jump = false;
		
		for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uleft1[i];

		SweepEff(LEFT2RIGHT, _m, _dt, temp, false, true);
		
		if(do_jump == true) {
//std::cout << "c0\n";
			assert(Rright == 1);
			for(int i = 0; i < Lsize-1; i++) texp[i] = II;
		} else {
//std::cout << "d0\n";
			for(int i = 0; i < Lsize-1; i++) texp[i] = UST.Uright1[i];
		}
		
		SweepEff(RIGHT2LEFT, _m, _dt, temp1, true, true);

		time += _dt;
		//std::cout << temp1 << "\n";
		//Mz += temp1;
		Obs[count_t] = temp1;
		count_t += 1;
//std::cout << count_t-1 << " " << int(time/_dt) << " " << int(tmaxm/_dt) <<" \n";

	}
	//Mz = Mz / complex(float(_n_qtraj), 0.);
	std::cout << '\n' ;
	//outf << time << "	" << Mz.real() << '\n';
}
*/	
