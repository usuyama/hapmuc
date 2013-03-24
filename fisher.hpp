#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm> // for min and max

using namespace boost::math;
using namespace std;

double fisher_test(unsigned a, unsigned b, unsigned c, unsigned d) {
		unsigned N = a + b + c + d;
		unsigned r = a + c;
		unsigned n = c + d;
		unsigned max_for_k = min(r, n);	
		unsigned min_for_k = (unsigned)max(0, int(r + n - N));
		hypergeometric_distribution<> hgd(r, n, N);	
		double cutoff = pdf(hgd, c);
		double tmp_p = 0.0;
		for(int k = min_for_k;k < max_for_k + 1;k++) {
				double p = pdf(hgd, k);
				if(p <= cutoff) tmp_p += p;
		}
		return tmp_p;
}
