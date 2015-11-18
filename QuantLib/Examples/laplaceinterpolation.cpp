#include <ql/math/matrix.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/experimental/math/laplaceinterpolation.hpp>

using namespace QuantLib;

int main() {

    const int N = 25;
    const double delPerc = 0.9;

    Matrix sample(N,N);

    for(Size i=0;i<N;++i) {
        for(Size j=0;j<N;++j) {
            sample[i][j] = sin((double)i/N*4.0)*cos((double)j/N*4.0);
        }
    }

    // std::cout << "Original Matrix:" << std::endl;
    // std::cout << sample << std::endl;

    MersenneTwisterUniformRng mt(42);

    Matrix sample2(sample), del(N,N);

    for(Size l=0;l<(N*N*delPerc);++l) {
        Size i = static_cast<Size>(mt.nextReal() * N);
        Size j = static_cast<Size>(mt.nextReal() * N);
        sample2[i][j] = Null<Real>();
        del[i][j] = 1;
    }

    // std::cout << "Randomly erased points:" << std::endl;
    // std::cout << sample << std::endl;

    // std::cout << "Reconstructed matrix:" << std::endl;
    laplaceInterpolation(sample2,1E-6);
    // std::cout << sample << std::endl;

    for(Size i=0;i<N;++i) {
        for(Size j=0;j<N;++j) {
            std::cout << i << " " << j << " " << sample[i][j] << " " << sample2[i][j] << " " << del[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;

}
