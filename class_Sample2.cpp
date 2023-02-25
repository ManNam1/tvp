#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

class classSample {
	std::string sNam;
	bool qPhas;
	std::vector<std::string> uNam;
	std::vector<double> uCom;
	std::vector<std::string> iNam;
	std::vector<double> iCom;
	std::string lNam;
	double uPlsCN = 0, uPlsMW = 0, uPlsSG = 0;
	double spcH2 = 0, spcHE = 0, spcN2 = 0, spcCO2 = 0, spcH2S = 0, spcNeo = 0;
	double scnC7 = 0, scnC8 = 0, scnC9 = 0, scnC10 = 0;
	double nCtot = 0;
	double iPlsMW = 0, iPlsSG = 0, iPlsAL = 0;
	int nCom = 0, nUsr = 0;
	std::vector<double> sCom;

	std::vector<double> BubT, BubP, DewT, DewP, FlsT, FlsP;
	double pCri  = 0, tCri =  0;
	double TBmin = 0, TBmax = 0, PBmax = 0;
	double TDmin = 0, TDmax = 0, PDmax = 0;
	double pBar  = 0, tHrm =  0;

public:
	classSample(std::string sNam) {
		this->sNam = sNam;
		this->qPhas = false;
	}

	// Internal Names & Moles
	void setIntName(std::string iNam) { this->iNam.push_back(iNam); }
	void setIntMoles(double iCom) { this->iCom.push_back(iCom); }

	// Sample Long Name
	void setLongName(std::string lNam) { this->lNam = lNam; }

	// User Sample Plus Fraction Carbon Number, Mole Weight and Specific Gravity   
	void setUserPlus(double uPlsCN, double uPlsMW, double uPlsSG) {
		this->uPlsCN = uPlsCN;
		this->uPlsMW = uPlsMW;
		this->uPlsSG = uPlsSG;
	}

	//-- Special Components ----------------------------------------------
	void setSpecial(double spcH2, double spcHE, double spcN2, double spcCO2, double spcH2S, double spcNeo) {
		this->spcH2 = spcH2;
		this->spcHE = spcHE;
		this->spcN2 = spcN2;
		this->spcCO2 = spcCO2;
		this->spcH2S = spcH2S;
		this->spcNeo = spcNeo;
	}

	//-- Number of Components in SCN 7 thru 10 ---------------------------
	void setSCN(double scnC7, double scnC8, double scnC9, double scnC10) {
		this->scnC7 = scnC7;
		this->scnC8 = scnC8;
		this->scnC9 = scnC9;
		this->scnC10 = scnC10;
	}

	//-- Total Number of User Specified Components -----------------------
	void setTotCom(int nCtot) {
		this->nCtot = nCtot;
	}

	//-- Set User Names & Compositions -----------------------------------
	void setUserName(std::string uNam) {
		this->uNam.push_back(uNam);
	}
	void setUserComp(double uCom) {
		this->uCom.push_back(uCom);
	}

	//-- Set Int Plus Fraction Properties --------------------------------
	void setIntPlusMW(double MW) { iPlsMW = MW; }
	void setIntPlusSG(double SG) { iPlsSG = SG; }
	void setIntPlusAL(double AL) { iPlsAL = AL; }

	// Number of Components and Composition
	void setIntComp(int nCom, int nUsr) {
		this->nCom = nCom;
		this->nUsr = nUsr;
		this->sCom.resize(nCom); // assuming sCom is a vector in C++
	}

	//-- Internal Name & Moles -------------------------------------------
	void setInternal(std::string iNam, double iMol) {
		this->iNam.push_back(iNam);
		this->iCom.push_back(iMol);
	}

	//-- Set and Get a Component's Moles ---------------------------------
	void sZI(int iC, double zI) {
		sCom[iC] = zI;
	}

	double gZI(int iC) {
		return sCom[iC];
	}

	double gTot() {
		double sTot = 0.0;
		for (int iC = 0; iC < nCom; iC++) {
			sTot += sCom[iC];
		}
		return sTot;
	}

	// Values from the Michelsen Approx Phase Envelope Calculation
	void setTempPsat(
		std::vector<double>& BubT, std::vector<double>& BubP,
		std::vector<double>& DewT, std::vector<double>& DewP,
		std::vector<double>& FlsT, std::vector<double>& FlsP,
		double pCri, double tCri)
	{
		qPhas = true;
		this->BubT = BubT;
		this->BubP = BubP;
		this->DewT = DewT;
		this->DewP = DewP;
		this->FlsT = FlsT;
		this->FlsP = FlsP;

		int nBub = BubT.size();
		int nDew = DewT.size();
		int nFls = FlsT.size();

		// Bubble Point Line
		if (nBub > 0) {
			double TBmin = *std::min_element(BubT.begin(), BubT.end());
			double TBmax = *std::max_element(BubT.begin(), BubT.end());
			double PBmax = *std::max_element(BubP.begin(), BubP.end());
			this-> TBmin = TBmin;
			this-> TBmax = TBmax;
			this-> PBmax = PBmax;
		}
		else {
			this->TBmin = std::numeric_limits<double>::quiet_NaN();
			this->TBmax = std::numeric_limits<double>::quiet_NaN();
			this->PBmax = std::numeric_limits<double>::quiet_NaN();
		}

		// Dew Point Line
		if (nDew > 0) {
			double TDmin = *std::min_element(DewT.begin(), DewT.end());
			double TDmax = *std::max_element(DewT.begin(), DewT.end());
			double PDmax = *std::max_element(DewP.begin(), DewP.end());
			this-> TDmin = TDmin;
			this-> TDmax = TDmax;
			this-> PDmax = PDmax;
		}
		else {
			this->TDmin = std::numeric_limits<double>::quiet_NaN();
			this->TDmax = std::numeric_limits<double>::quiet_NaN();
			this->PDmax = std::numeric_limits<double>::quiet_NaN();
		}

		// Critical Point
		this->pCri = pCri;
		this->tCri = tCri;

		// CricondenBar
		if (!std::isnan(PBmax) && !std::isnan(PDmax)) {
			pBar = std::max(PBmax, PDmax);
		}
		else if (!std::isnan(PBmax)) {
			pBar = PBmax;
		}
		else if (!std::isnan(PDmax)) {
			pBar = PDmax;
		}
		else {
			pBar = std::numeric_limits<double>::quiet_NaN();
		}

		// CricondenTherm
		if (!std::isnan(TDmax)) {
			tHrm = TDmax;
		}
		else if (!std::isnan(TBmax)) {
			tHrm = TBmax;
		}
		else {
			tHrm = std::numeric_limits<double>::quiet_NaN();
		}
	}


};


