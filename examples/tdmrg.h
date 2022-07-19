#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <complex>
#include <math.h>
#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
//#include <Eigen/Geometry>		// m.reshaped
#include <Spectra/SymEigsSolver.h>


using namespace Eigen;
typedef std::complex<double> complex;

const int LEFT=-1, RIGHT=1;
const int LEFT2RIGHT=2, RIGHT2LEFT=-2;

double g, h;

//////////////////////////////////////////////////////////////////////////////////

MatrixXcd kron_prod( const MatrixXcd& bl1, const MatrixXcd& bl2){

	MatrixXcd res(bl1.rows()*bl2.rows(), bl1.cols()*bl2.cols());
	for(int i=0; i<bl1.rows(); ++i){
		for(int j=0; j<bl1.cols(); ++j){
			//.block(x_initpoint, y_initpoint, x_size, y_size)
			res.block(i*bl2.rows(), j*bl2.cols(), bl2.rows(), bl2.cols())=bl1(i, j)*bl2;
			//std::cout << res.block(i,j,i+2,j+2) << '\n';
		}
	}
	return(res);
}

//////////////////////////////////////////////////////////////////////////////////

class Bblock {
	public:
		MatrixXcd H, U, Id;
		ArrayX< MatrixXcd > Sx, Sz;

		Bblock(void) { resize(1); }
		Bblock(int _size) { resize(_size); };
		//Bblock update(int direction); // char == int direction;
		void reduction(const MatrixXcd& _U); // U==vec;
		Bblock& resize(int n) {
			Sz.resize(n);
			Sx.resize(n);
			return *this;
		}
		//MatrixXcd printA(void);
};
	
void Bblock::reduction(const MatrixXcd& _U){
	assert(H.cols() ==_U.rows());
	H = _U.adjoint()*( H * _U );
	Id = _U.adjoint()*( Id * _U );
	int Ssize = Sx.size();
	assert(Ssize == Sz.size());
	for(int i=0; i < Ssize; ++i) {
		Sx[i] = _U.adjoint()*( Sx[i] * _U );
		Sz[i] = _U.adjoint()*( Sz[i] * _U );
	}
	U = _U;
}

///////////////////////////////////////////////////////////////////////////////////
class SToper {
	public:
		double dt_ldiss[2];
		ArrayX< MatrixXcd > Uleft1, Uright1;
		//ArrayX< MatrixXcd > Uleft2, Uright2;
		
		//SToper(int _Lsize, double _dt);
};
///////////////////////////////////////////////////////////////////////////////////

class chain {
	public:
		MatrixXcd Sx0;
		MatrixXcd Sz0;
		ArrayX<Bblock > blockL;
		ArrayX<Bblock > blockR;
		VectorXcd gs;
		ArrayX< MatrixXcd > texp;
		SToper UST;	
		complex energy;
		int Lsize, Rright, Lleft;
		double time;
	
		chain(int _L) {
			Sx0.resize(2,2);
			Sz0.resize(2,2);
			Sx0 <<  0, 1,
		  		1, 0;
			Sz0 <<  1, 0,
		  	        0, -1;
		  	blockL.resize(_L);
		  	blockR.resize(_L);
			blockL[0].H = -g*Sz0 - h*Sx0;	
			blockL[0].U = MatrixXcd::Identity(2,2);
			blockR[0].H = -g*Sz0 - h*Sx0;	
			blockR[0].U = MatrixXcd::Identity(2,2);
			blockL[0].Sx[0] = Sx0;
			blockR[0].Sz[0] = Sz0;
			blockR[0].Sx[0] = Sx0;
			blockL[0].Sz[0] = Sz0;
			blockR[0].Id = MatrixXcd::Identity(2,2);
			blockL[0].Id = MatrixXcd::Identity(2,2);	
		  	texp.resize(_L-1);
		  	for(int j=0; j<_L-2; ++j) {
		  		texp[j].resize(4,4);
		  		}
		  	Lsize = _L;
		  	time = 0.;
		  }
		  void BuildBlock(int direction, int iter);
		  void GroundState(void);
		  void Truncation(int direction, const MatrixXcd& _U);
		  void BuildU(int direction, int _m);
		  void gs_ReBuild(int verse);
		  complex obser(void);
		  MatrixXcd BondExp(double _tstep, const char& chr, double _w);
		  void InitialEvol(void);
		  void SuzukiTrotter(void);
		  void Sweep(int verse, int _m, bool measure, bool evol, std::ofstream& outf);
		  void FiniteDMRG(int n_sweeps, int _m, std::ofstream& outf);
		  //////////////////////////////////////////////////
		  void BuildBlockEff(int direction, int iter);
		  void InitialEvolEff(void);
		  void initSToper(double _dt, double _w);
		  //void SuzukiTrotterEff(int verse, double _dt);
		  void Jump(void);
		  void SweepEff(int verse, int _m, complex& _Mz, bool measure, bool evol);
		  void QTraj(int _m, double _dt, double tmaxm, VectorXcd& Obs);
		  //////////////////////////////////////////////////
		  void derivate(double t, MatrixXcd& rho, MatrixXcd& k);
		  void RungeKutta(double _dt);
		  void SweepRK(int verse, int _m, double _dt, bool measure, bool evol, std::ofstream& outf);
		  //////////////////////////////////////////////////
		  void Blockderivate(double t, const MatrixXcd& rho, MatrixXcd& k, const MatrixXcd& Ham, const MatrixXcd& L1);
		  void BlockRungeKutta(double _dt);
		  //////////////////////////////////////////////////
		  void QTrajAdp(int _m, double _dt, double dpp, double dp, double tmaxm, VectorXcd& Obs, VectorXd& dataT);

		  //void QTraj(int _m, int _n_qtraj, std::ofstream& outf);
		  //std::ofstream& outf
		  //friend class SToper;
};


void chain::BuildBlock(int direction, int iter) {
	MatrixXcd I2 = MatrixXcd::Identity(2,2);

	if(direction == LEFT) {
		Lleft = iter;
		blockL[iter].resize(iter+1);

		MatrixXcd HL = blockL[iter-1].H, SxL = blockL[iter-1].Sx[iter-1];
		//int diml = HL.rows();
		//MatrixXcd Ileft = MatrixXcd::Identity(diml, diml);
		MatrixXcd Ileft = blockL[iter-1].Id;

		blockL[iter].H = kron_prod(HL,I2) - g* kron_prod(Ileft,Sz0) - h* kron_prod(Ileft,Sx0) - kron_prod(SxL,Sx0);
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

		blockR[iter].H = kron_prod(I2, HR) - g* kron_prod(Sz0, Iright) - h* kron_prod(Sx0, Iright) - kron_prod(Sx0, SxR);
		blockR[iter].Id = kron_prod(I2, Iright);
		for(int i=0; i<iter; ++i) {
			blockR[iter].Sx[i] = kron_prod(I2, blockR[iter-1].Sx[i]);
			blockR[iter].Sz[i] = kron_prod(I2, blockR[iter-1].Sz[i]);
		}
		blockR[iter].Sz[iter] = kron_prod(Sz0, Iright);
		blockR[iter].Sx[iter] = kron_prod(Sx0, Iright);
	
	}
}
		
		

void chain::GroundState(void) {
		
	int diml=blockL[Lleft].H.rows(), dimr=blockR[Rright].H.rows();		
	int dimt = diml*dimr;
	MatrixXcd Hsuper(dimt, dimt);
	//MatrixXcd IR=MatrixXcd::Identity(dimr, dimr), IL=MatrixXcd::Identity(diml, diml);
	MatrixXcd IL=blockL[Lleft].Id, IR=blockR[Rright].Id;
	MatrixXcd HL=blockL[Lleft].H, HR=blockR[Rright].H;
	MatrixXcd SxL=blockL[Lleft].Sx[Lleft], SxR=blockR[Rright].Sx[Rright];
	
	assert( dimt == SxL.rows()*SxR.rows() );
	
	Hsuper = kron_prod(HL, IR) + kron_prod(IL, HR) - kron_prod(SxL,SxR);	
		
	int lz_iter = 8;
	using namespace Spectra;
		MatrixXd Hsupe(Hsuper.rows(), Hsuper.cols());
		Hsupe = Hsuper.real();
		DenseSymMatProd< double > op(Hsupe);
		SymEigsSolver< DenseSymMatProd< double > > eigs(op, 1, lz_iter);
		// Initialize and compute
		eigs.init();
		int nconv;	nconv = eigs.compute(SortRule::SmallestAlge);
		energy = eigs.eigenvalues()(0);
		//Hsupe = (eigs.eigenvectors());
		//psi = Hsupe.block(0,0,Hsupe.rows(),1);
		gs = eigs.eigenvectors();
		gs.normalize();
}

void chain::Truncation(int direction, const MatrixXcd& _U) {

	if(direction == LEFT)
		blockL[Lleft].reduction(_U);
	else if(direction == RIGHT)
		blockR[Rright].reduction(_U);
	else {
		blockL[Lleft].reduction(_U);
		blockR[Rright].reduction(_U);
	}

}		
		
void chain::BuildU(int direction, int _m) {		
	int dim1=blockL[Lleft].H.rows(), dim2=blockR[Rright].H.rows();
	int dim;
	MatrixXcd rho1(dim1, dim2), rho;
	
	assert(gs.size() == dim1*dim2);
	
	rho1 = gs.reshaped<RowMajor>(dim1, dim2);
	if ( direction == RIGHT ) {
		rho = rho1.adjoint() * rho1;
	}
	else {
		rho = rho1 * rho1.adjoint();
	}
	MatrixXcd vec(rho.rows(), rho.cols());
	/*
	/////////////////////////////////////////////////////////////////////////
	std::ofstream outdata;
	outdata.open("temp.dat", std::ios::trunc);
	outdata.precision(9);
	for(int i=0; i<rho.rows();++i){
		for(int j=0; j<rho.cols(); ++j){
			outdata << rho(i, j).real();
			if(rho(i, j).imag()>=0.) outdata << '+';
			outdata << rho(i, j).imag() << "j    ";
		}
		outdata << '\n';
	}
	outdata.close();
	
	//std::cout << rho << '\n';
	system("python3 eigen.py");
	std::ifstream indata;
	indata.open("temp.dat");
	double a, b;
	for(int i=0; i<rho.rows();++i){
		for(int j=0; j<rho.cols(); ++j) {
			indata >> a >> b;
			vec(i,j) = complex(a, b); }
	}
	indata.close();
	assert(rho.cols()==vec.cols());		assert(vec.rows() == rho.cols());
	//std::cout << vec << '\n'; exit (8);
	/////////////////////////////////////////////////////////////////////////
	*/

	SelfAdjointEigenSolver<MatrixXcd> eigenrho(rho);

	if(eigenrho.info() != Success) abort();

	vec = eigenrho.eigenvectors();

	//std::cout << eigenrho.eigenvalues();
	if ( direction == LEFT) {
		dim = std::min(_m, dim1);
	}
	else if ( direction == RIGHT) {
		dim = std::min(_m,dim2);
	}
	else {
		dim = std::min(_m, dim1);
	}

	MatrixXcd _U;
	_U = vec.block(0, vec.cols()-dim, vec.rows(), dim);
	
	Truncation(direction, _U);
		
}


void chain::gs_ReBuild(int verse) {

	if(verse == LEFT2RIGHT) {
		if(Lleft == 1) {
			assert(Rright == Lsize-Lleft-2);
			assert(gs.size() == blockL[Lleft].H.rows()*blockR[Rright].H.rows());
			return;
		}
		else {
			MatrixXcd I2 = MatrixXcd::Identity(2,2);
			assert(Rright == Lsize-Lleft-2);
			MatrixXcd vec_mem=blockR[Rright].U;
			MatrixXcd vec_old=blockL[Lleft-1].U;
			
			vec_old = kron_prod( vec_old.adjoint(), I2 );
			vec_old = kron_prod( vec_old, vec_mem);
			assert(vec_old.cols() == gs.size());
			gs = vec_old * gs;
			assert(gs.size() ==  blockL[Lleft].H.rows()*blockR[Rright].H.rows());
		}
		
	} else {
		if(Rright == 1) {
			assert(Lleft == Lsize-Rright-2);
			assert(gs.size() == blockR[Rright].H.rows()*blockL[Lleft].H.rows());
			return;
		}
		else {
			assert(Lleft == Lsize-Rright-2);
			MatrixXcd I2 = MatrixXcd::Identity(2,2);
			MatrixXcd vec_mem=blockL[Lleft].U;
			MatrixXcd vec_old=blockR[Rright-1].U;
			
			vec_old = kron_prod( I2, vec_old.adjoint() );
			vec_old = kron_prod( vec_mem, vec_old );
			assert(vec_old.cols() == gs.size());
			gs = vec_old * gs;
			assert(gs.size() ==  blockL[Lleft].H.rows()*blockR[Rright].H.rows());
		
		}
	}
}

inline complex chain::obser(void) {
	complex output;
	int diml = blockL[Lleft-1].Id.rows();
	//MatrixXcd IL=MatrixXcd::Identity(diml, diml);
	MatrixXcd IL = kron_prod( blockL[Lleft-1].Id, MatrixXcd::Identity(2,2) );
	MatrixXcd SzR = kron_prod( Sz0, blockR[Rright-1].Id);
 
	assert( Lleft == Lsize - Rright - 2 );
	assert( diml*2*2*blockR[Rright-1].Id.rows() == gs.size() );
	// SzR = blockR[Rright].Sz[Rright]
	output = gs.adjoint() * (kron_prod(IL, SzR) * gs);
	return(output);
}

/////////////////////////////////////////////////////////////////////////////////
MatrixXcd chain::BondExp(double _tstep, const char& chr, double _w = 0.) {

	MatrixXcd Uij(4,4);
	MatrixXcd Hij(4,4), U(4,4), Ut(4,4), aux(4,4);
	VectorXcd wk(4);
	MatrixXcd I2 = MatrixXcd::Identity(2,2);
/*	// SUZUKI 2&4 ORDER
	if(chr=='D')
		Hij = - g*kron_prod(I2,Sz0) - h*kron_prod(I2,Sx0);
	else
		Hij = -g*kron_prod(Sz0,I2)-h*kron_prod(Sx0,I2)-kron_prod(Sx0,Sx0);
*/
	// SUZUKI 1 ORDER
	if(chr=='D')
		Hij = - g*kron_prod(Sz0,I2) - h*kron_prod(Sx0,I2) - g*kron_prod(I2,Sz0) - h*kron_prod(I2,Sx0) - kron_prod(Sx0,Sx0);
	/*else if (chr == 'Q')
		Hij =  g*kron_prod(Sz0,I2) - h*kron_prod(Sx0,I2) - kron_prod(Sx0,Sx0) - complex(0., 1.)*_w/2.*kron_prod(I2,I2);*/ //DON'T WORK
	else
		Hij = - g*kron_prod(Sz0,I2) - h*kron_prod(Sx0,I2) - kron_prod(Sx0,Sx0);
		
return Hij;	// -> for <<<RUNGE-KUTTA method>>>
		
	// PRIMORDIAL EVOLUTION OPERATOR
//	Uij = MatrixXcd::Identity(4,4) - complex(0., _tstep) * Hij;

/*
	SelfAdjointEigenSolver<MatrixXcd> eigenhijcd(Hij);
	wk = eigenhijcd.eigenvalues();
	U = eigenhijcd.eigenvectors();

	Ut = U.adjoint();

	Hij = MatrixXcd::Zero(4,4);

	for(int i = 0; i < 4; i++)	Hij(i,i) = std::exp(complex(0.,0.) - complex(0.,1.)*wk(i)*complex(_tstep, 0.) );

	aux = Hij * Ut;
	Uij = U * aux;
*/

//return Uij;	// -> for <<<QUANTUM TRAJECTORIES>>> algorithm
}
/////////////////////////////////////////////////////////////////////////////////

void chain::InitialEvol(void) {

	blockL[0].H = -g*Sz0 - h*Sx0;
	blockR[0].H = -g*Sz0 - h*Sx0;
}

void chain::SuzukiTrotter(void) {

	//if(verse == LEFT2RIGHT) {
		if(Lleft == 1) {
			int dim_r = blockR[Rright].Id.rows();
			assert(gs.size() == dim_r*4);
			//MatrixXcd I_right = MatrixXcd::Identity(dim_r, dim_r);
			MatrixXcd I_right = blockR[Rright].Id;
			MatrixXcd Uev(dim_r*4,dim_r*4);
			VectorXcd aux(gs);
			Uev = kron_prod(texp[0], I_right);
			gs = Uev * aux;
			gs.normalize();
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



void chain::Sweep(int verse, int _m, bool measure, bool evol, std::ofstream& outf){
	complex Mz;
	if(verse == LEFT2RIGHT){
		for(int iter=1; iter < Lsize-2; ++iter) {
			BuildBlock(LEFT, iter);	
			BuildBlock(RIGHT, Lsize - iter - 2);
			
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
			BuildBlock(RIGHT, iter);	
			BuildBlock(LEFT, Lsize - iter - 2);
			
			if(evol) {
				gs_ReBuild(verse);
				SuzukiTrotter();
			} else {
				GroundState();
				}
			if(measure && Lleft == 1) {
				Mz=obser();
				if ((Mz.real() > 1.) || (Mz.real() < 0.)){
					std::cerr << "Overflow Error\n";}//exit (8);}

				outf << time << "	" << Mz.real() << '\n';
				std::cout << Mz << '\n';
				std::cout << energy << "   eNeRgY \n";
			}

			BuildU(RIGHT, _m);
		}
	}
}






void chain::FiniteDMRG(int n_sweeps, int _m, std::ofstream& outf) {
	assert(n_sweeps > 0);
	// INFINITE DMRG
	for (int iter=1; iter <= Lsize/2; ++iter) {
		BuildBlock(LEFT, iter);
		BuildBlock(RIGHT, iter);

		GroundState();
		BuildU(0, _m);
	}
	// FINITE DMRG
	int first_iter = Lsize/2 + 1;

	for(int iter=first_iter; iter < Lsize-2; ++iter) {
		BuildBlock(LEFT, iter);	
		BuildBlock(RIGHT, Lsize - iter - 2);
		GroundState();
		BuildU(LEFT, _m);
	}
	
	bool measure;
	measure = true;
	
	Sweep(RIGHT2LEFT, _m, measure, false, outf);

	for (int sweep=1; sweep <= n_sweeps-1; ++sweep) {
		Sweep(LEFT2RIGHT, _m, false, false, outf);
		Sweep(RIGHT2LEFT, _m, true, false, outf);
	}
}
