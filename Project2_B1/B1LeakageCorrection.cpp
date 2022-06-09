#include "B1LeakageCorrection.h"
#include "define.h"


B1LeakageCorrection::B1LeakageCorrection()
{
	bConverge = true;
}


B1LeakageCorrection::~B1LeakageCorrection()
{

}

int B1LeakageCorrection::b1man(string inName)
{
	int i, j;
	ofstream fout1("output.txt");
	readInput(inName);

	for (i = 0; i < nGrp; i++) {
		MMat(i, i) = xs_total(i);
	}
	for (i = 0; i < nGrp; i++) {
		for (j = 0; j < nGrp; j++) {
			MMat(i,j) = MMat(i, j)-1.0*sm_p0(i,j);
		}
	}
	FMat = chi*nu_fis.transpose();
	AMat = MMat.inverse()*FMat;

	double keff0=maxeig();
	double kinf=keff0;
	
	Bsq = 0;
	for (i = 0; i < nGrp; i++) phi(i) = 1.;

	double Bsqd = 0;
	double keffd, keff = keff0;
	Bsq = 0.00000000001;
	double Bsqn;

	int itr = 0;
	keffd = keff0;

	do {
		itr = itr + 1;
		calAlpha();
		setMatrixD();
		MMat = -sm_p0 + Bsq*DMat;
		for (i = 0; i < nGrp; i++) {
			MMat(i, i) += xs_total(i);
		}
		AMat = MMat.inverse()*FMat;
		keff=maxeig();
		cout  << itr << '\t' << Bsq << '\t' << keff << endl;
		fout1 << itr << '\t'<<Bsq << '\t' << keff << endl;
		
		Bsqn = (1 - 1 / keff) / (1 / keff - 1 / keffd)*(Bsq - Bsqd) + Bsq;
		
		Bsqd = Bsq; Bsq = Bsqn;
		keffd = keff;
		bConverge = (abs(keff - 1.)) > 1e-5;
	} while (bConverge);
	VectorXd Current;
	Current = VectorXd::Zero(nGrp);

	absB = sqrt(abs(Bsq));
	Current = DMat*absB*phi;


	DifCoeffi = (DMat*phi);
	for (i = 0; i < nGrp; i++) {
		DifCoeffi(i) /= phi(i);
	}

	fout1 << endl;
	fout1.setf(ios::scientific);
	fout1 << setprecision(5);
	int width=12;
	fout1 << "Grp"<< setw(width) << "upperEB" << setw(width) << "inf_spec " << setw(width+10) << "normalizedFlux/delu" << setw(width+10)<<"normalizedJ/delu"<<endl;
	width = 15;
	for (i = 0; i < nGrp-1; i++) {
		fout1 << setw(2) << i << setw(width) << upperEB(i) << setw(width) << inf_spec(i)/inf_spec.sum() / log(upperEB(i) / upperEB(i + 1)) << setw(width) << phi(i) / (phi.sum()) / log(upperEB(i) / upperEB(i + 1)) <<  setw(width) << Current(i) / (Current.sum()) / log(upperEB(i) / upperEB(i + 1)) <<endl;
	}
	i = nGrp - 1;
	fout1 << setw(2) <<i<< setw(width)  << upperEB(i) << setw(width) << inf_spec(i) / inf_spec.sum() / log(10) << setw(width) << phi(i) / (phi.sum()) /log(10)<< setw(width) << Current(i) / (Current.sum()) / log(10) << endl;

	fout1 << endl;
	fout1 << "DifCoeff" << endl;
	for (i = 0; i < nGrp ; i++) {
		fout1 << setw(width) <<DifCoeffi(i)<<setw(width)<<1./xs_trans(i)/3<< endl;
	}
	return 0;
}

int B1LeakageCorrection::readInput(string inName)
{
	ifstream fin(inName);
	string str, str1;
	int i, j, k;

	do {
		fin >> str;
	} while (str != "Scattering_TrC");
	fin >> str;
	do {
		str1 = str;
		for (i = 0; i < 11; i++)	fin >> str;
	} while (str != "from");
	nGrp = atof(str1.c_str());
	fin.close();

	fin.open(inName);
	do {
		fin >> str;
	} while (str != "Scattering_TrC");

	allocateVariable();

	for (i = 0; i < nGrp; i++) {
		fin >> str;

		fin >> str;	upperEB(i) = atof((str.c_str()));
		fin >> str; xs_abs(i) = atof(str.c_str());
		fin >> str; xs_scat(i) = atof(str.c_str());
		fin >> str; nu_fis(i) = atof(str.c_str());
		fin >> str; chi(i) = atof(str.c_str());
		fin >> str; inf_spec(i) = atof(str.c_str());
		fin >> str; xs_total(i) = atof(str.c_str());
		fin >> str; xs_trans(i) = atof(str.c_str());
		fin >> str; smgg_TrC(i) = atof(str.c_str());
		fin >> str; xs_scat_TrC(i) = atof(str.c_str());
	}
	fin >> str >> str >> str >> str; fin >> str;
	while (!fin.eof()) {
		i = atoi(str.c_str());
		fin >> str;
		j = atoi(str.c_str());
		fin >> str;		sm_p0(j - 1, i - 1) = atof(str.c_str());
		fin >> str;		sm_p1(j - 1, i - 1) = atof(str.c_str());
		fin >> str;
	}
	return 0;
}

int B1LeakageCorrection::allocateVariable()
{


	//upperEB.resize(nGrp);
	upperEB = VectorXd::Zero(nGrp);
	xs_abs = VectorXd::Zero(nGrp);
	xs_scat = VectorXd::Zero(nGrp);
	nu_fis = VectorXd::Zero(nGrp);
	chi = VectorXd::Zero(nGrp);

	inf_spec = VectorXd::Zero(nGrp);

	xs_total = VectorXd::Zero(nGrp);
	xs_trans = VectorXd::Zero(nGrp);
	smgg_TrC = VectorXd::Zero(nGrp);
	xs_scat_TrC = VectorXd::Zero(nGrp);

	sm_p0 = MatrixXd::Zero(nGrp, nGrp);
	sm_p1 = MatrixXd::Zero(nGrp, nGrp);

	//
	alpha = VectorXd::Zero(nGrp);

	pre_phi = VectorXd::Zero(nGrp);
	phi = VectorXd::Zero(nGrp);

	invDMat = MatrixXd::Zero(nGrp, nGrp); 	//D^(-1),  inverse of D matrix
	DMat = MatrixXd::Zero(nGrp, nGrp); 	//D,   D matrix
	MMat = MatrixXd::Zero(nGrp, nGrp);   //M    M matrix

	DifCoeffi = VectorXd::Zero(nGrp);

	AMat = MatrixXd::Zero(nGrp, nGrp);
	FMat = MatrixXd::Zero(nGrp, nGrp);


	return 0;
}

int B1LeakageCorrection::calAlpha()
{
	int i;
	double A00;

	absB = sqrt(abs(Bsq));

	if (Bsq > 0) {
		for (i = 0; i < nGrp; i++) {
			A00 = atan(absB / xs_total(i)) / absB;
			alpha(i) = Bsq*A00 / 3 / xs_total(i) / (1 - xs_total(i) * A00);
		}
	}
	else if (Bsq < 0) {
		for (i = 0; i < nGrp; i++) {
			A00 = atan(absB / xs_total(i)) / absB;
			alpha(i) = (log((xs_total(i) + absB) / (xs_total(i) - absB))) / (2 * absB);
		}
	}
	else for (i = 0; i < nGrp; i++) alpha(i) = 0.;

	return 0;
}

int B1LeakageCorrection::setMatrixD()
{
	int i, j;

	invDMat = -sm_p1;
	for (i = 0; i < nGrp; i++) {
		invDMat(i, i) += alpha(i)*xs_total(i);
	}
	invDMat = 3 * invDMat;
	DMat = invDMat.inverse();

	return 0;
}

int B1LeakageCorrection::setMatrixM()
{
	int i, j;
	MMat = -sm_p0 + Bsq*DMat;
	for (i = 0; i < nGrp; i++) {
		MMat(i, i) += xs_total(i);
	}
	return 0;
}

double B1LeakageCorrection::maxeig()
{
	int i = 0;
	VectorXd xd;
	for (i = 0; i < nGrp; i++) phi(i) = 1.;
	double lam = 1, lamd, sigma;
	for (i = 0; i<1000; i++)
	{
		xd = phi;
		phi = AMat*xd / lam;
		//	   x=A*x/lam;
		lamd = lam;
		lam = lamd*sqrt((phi.dot(phi) / (phi.dot(xd))));
		if (abs(lam - lamd)<1e-7)
			break;
	}
	sigma = 1 - lamd / lam;
	lam = lam / (1 - sigma);
	if (i == 1000)
		std::cout << " eigenvalue not converged within 1000 iterations" << std::endl;
	return lam;
}

int B1LeakageCorrection::changeBuckling()
{
	Bsq = Bsq + (Bsq - pre_Bsq) / (source - pre_source)*(1 - source);

	return 0;
}