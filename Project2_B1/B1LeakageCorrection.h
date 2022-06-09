#pragma once

#include "define.h"

class B1LeakageCorrection
{
private:
	//variable from input file;
	int nGrp;
	
	VectorXd upperEB;
	VectorXd xs_abs;
	VectorXd xs_scat;
	VectorXd nu_fis;
	VectorXd chi;
	
	VectorXd inf_spec;
	
	VectorXd xs_total;
	VectorXd xs_trans;
	VectorXd smgg_TrC;
	VectorXd xs_scat_TrC;

	MatrixXd sm_p0;
	MatrixXd sm_p1;
private:
	double absB; // |B|  은 그때 그때 필요할때 Bsq를 통해 구한다.
	double Bsq;	//B^2 을 들고다닌다.
	double source;	// psi

	VectorXd alpha; 

	VectorXd pre_phi;
	VectorXd phi;  

	MatrixXd invDMat;	//D^(-1),  inverse of D matrix
	MatrixXd DMat;	//D^(-1),  inverse of D matrix

	MatrixXd MMat;

	bool bConverge;

	VectorXd DifCoeffi;

	double pre_Bsq;
	double pre_source;

	MatrixXd AMat;
	MatrixXd FMat;



public:
	B1LeakageCorrection();
	~B1LeakageCorrection();

	int b1man(string inName);

	int readInput(string inName);
	int allocateVariable();

	int calAlpha();
	int setMatrixD(); // set invDMat and DMat
	int setMatrixM(); // set MMat

	double maxeig();

	int changeBuckling();

};

