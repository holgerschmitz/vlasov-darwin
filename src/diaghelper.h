#include <fstream>

#ifndef DIAGHELPER_H
#define DIAGHELPER_H

//using namespace std;

class Potential;

class ES_EFieldEnergy {
    private:
        int t;
        Potential *pField;
        double Energy;
        ofstream Output;
    public:
		void Init (Potential *pField_);
		void Execute (); 
};

#endif
