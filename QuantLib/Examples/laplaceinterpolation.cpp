#include <ql/math/matrix.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/experimental/math/laplaceinterpolation.hpp>

using namespace QuantLib;

int main() {

    int N = 100;
    double delPerc = 0.9;
    char *arg1 = getenv("DELPERC");
    if (arg1 != NULL)
        delPerc = atof(arg1);
    char *arg2 = getenv("N");
    if (arg2 != NULL)
        N = atoi(arg2);

    Matrix sample(N, N);

    for (Size i = 0; i < N; ++i) {
        for (Size j = 0; j < N; ++j) {
            sample(i, j) =
                sin((double)(i + 1) / N * 4.0) * cos((double)(j + 1) / N * 4.0);
        }
    }

    // std::cout << "Original Matrix:" << sample.size1() << "x" <<
    // sample.size2() << std::endl;
    // std::cout << sample << std::endl;

    MersenneTwisterUniformRng mt(42);

    Matrix sample2(sample), del(N, N);

    for (Size l = 0; l < (N * N * delPerc); ++l) {
        Size i = static_cast<Size>(mt.nextReal() * N);
        Size j = static_cast<Size>(mt.nextReal() * N);
        if (sample2[i][j] != Null<Real>()) {
            sample2[i][j] = Null<Real>();
            del[i][j] = 1;
        } else {
            --l;
        }
    }

    // std::cout << "Randomly erased points:" << std::endl;
    // std::cout << sample << std::endl;

    // std::cout << "Reconstructed matrix:" << std::endl;
    laplaceInterpolation(sample2, 1E-6);
    // std::cout << sample << std::endl;

    for (Size i = 0; i < N; ++i) {
        for (Size j = 0; j < N; ++j) {
            std::cout << i << " " << j << " " << sample[i][j] << " "
                      << sample2[i][j] << " " << del[i][j] << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
