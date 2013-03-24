#include <iostream>
#include <cmath>
#include <cstdlib>

 double logHypergeometricProb(double* logFacs, int a, int b, int c, int d);
 void initLogFacs(double* logFacs, int n);

double logFisherExactTest(int a, int b, int c, int d) {
		int n = a + b + c + d;
		double* logFacs = new double[n+1];
		initLogFacs(logFacs, n);
		 double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d);
		 double pFraction = 0;
		 for(int x=0; x <= n; ++x) {
				 if ( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) {
						 double l = logHypergeometricProb(logFacs,x,a+b-x,a+c-x,d-a+x);
						 if ( l <= logpCutoff ) pFraction += exp(l - logpCutoff);
				 }
		 }
		 double logpValue = logpCutoff + log(pFraction);
		 delete [] logFacs;
		 return logpValue;
}

void initLogFacs(double* logFacs, int n) {
		logFacs[0] = 0;
		for(int i=1; i < n+1; ++i) 
				logFacs[i] = logFacs[i-1] + log((double)i);
}

double logHypergeometricProb(double* logFacs, int a, int b, int c, int d) {
		return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d]
				- logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}
