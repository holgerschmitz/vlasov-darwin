#include <matrix.H>

class Helmholtz {
    protected:
        double dx[2];

        int gama,nu1,nu2;
        double epsilon;
    public:
        Helmholtz();
        void solve(   NumMatrix<double,2> &u,
	                  NumMatrix<double,2> &f,
                      NumMatrix<double,2> &lambda);

        double Dx() { 
            return dx[0];
        }
        double Dy() { 
            return dx[1];
        }
    private:
        void gauss( NumMatrix<double,2> &u,     
	                NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda);   
                     
        void defect(NumMatrix<double,2> &u,
                    NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda,
                    NumMatrix<double,2> &fn,
                    NumMatrix<double,2> &lambdan);
                                      
        void prolongate(NumMatrix<double,2> &u,
	                    NumMatrix<double,2> &un);
        
        void mgi(   NumMatrix<double,2> &u,     
	                NumMatrix<double,2> &f,
                    NumMatrix<double,2> &lambda);    
                    
        void normalize( NumMatrix<double,2> &f);
                
        double distance(NumMatrix<double,2> &u,
	                    NumMatrix<double,2> &f);
                        
        void boundary(NumMatrix<double,2> &u);
};
