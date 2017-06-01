#include <iostream>
#include <cstdlib>
#include <ctime>

#include <GeNuSys/element_traits.h>

#include <GeNuSys/vector.h>
#include <GeNuSys/sparse_vector.h>
#include <GeNuSys/matrix.h>
#include <GeNuSys/sparse_matrix.h>

#include <GeNuSys/p_norm.h>
#include <GeNuSys/frobenius_norm.h>
#include <GeNuSys/operator_norm.h>

#include <GeNuSys/linalg_algorithms.h>

#include <GeNuSys/polynomial.h>
#include <GeNuSys/lehmer_schur.h>

#include <GeNuSys/number_system.h>
#include <GeNuSys/radix_properties.h>
#include <GeNuSys/digit_set.h>
#include <GeNuSys/smith_hash.h>
#include <GeNuSys/numsys_traits.h>
#include <GeNuSys/simultaneous.h>

GeNuSys::LinAlg::Matrix<int> randomMatrix(unsigned int rows, unsigned int cols, unsigned int range)
{
    GeNuSys::LinAlg::Matrix<int> result(rows, cols);
    for (unsigned int i = 0; i < rows; ++i)
    {
        for (unsigned int j = 0; j < cols; ++j)
        {
            result.set(i, j, rand() % (2 * range + 1) - range);
        }
    }
    return result;
}

GeNuSys::LinAlg::Vector<int> randomVec(unsigned int length, unsigned int range)
{
    GeNuSys::LinAlg::Vector<int> result(length);
    for (unsigned int i = 0; i < length; ++i)
    {
        result.set(i, rand() % (2 * range + 1) - range);
    }
    return result;
}

template<typename ElementType2>
void exampleSmithNormalForm(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    std::cout << "===============================" << std::endl;
    std::cout << "SMITH NORMAL FORM EXAMPLE" << std::endl;

    std::cout << "M = " << std::endl << mat << std::endl;

    GeNuSys::LinAlg::SmithNormalForm<int> snf = GeNuSys::LinAlg::Algorithms::getSmithNormalForm(mat);

    std::cout << "U = " << std::endl << snf.U << std::endl;

    std::cout << "V = " << std::endl << snf.V << std::endl;

    std::cout << "S = " << std::endl << snf.S << std::endl;

    std::cout << "U * SM * V = " << std::endl << (snf.U * mat * snf.V) << std::endl;

    std::cout << "UInv = " << std::endl << GeNuSys::LinAlg::Algorithms::invert(snf.U) << std::endl;

    std::cout << "===============================" << std::endl;
}

template<typename ElementType2>
void exampleDecomposeQR(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    std::cout << "===============================" << std::endl;
    std::cout << "QR DECOMPOSE EXAMPLE" << std::endl;

    std::cout << "M = " << std::endl << mat << std::endl;

    GeNuSys::LinAlg::QR<ElementType2> qr = GeNuSys::LinAlg::Algorithms::decomposeQR(mat);

    std::cout << "Q = " << std::endl << qr.Q << std::endl;

    std::cout << "R = " << std::endl << qr.R << std::endl;

    std::cout << "QR = " << std::endl << qr.Q* qr.R << std::endl;

    std::cout << "===============================" << std::endl;
}

template<typename ElementType2>
void exampleDecomposeLU(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    std::cout << "===============================" << std::endl;
    std::cout << "LU DECOMPOSE EXAMPLE" << std::endl;

    std::cout << "M = " << std::endl << mat << std::endl;

    GeNuSys::LinAlg::LU<ElementType2> lu = GeNuSys::LinAlg::Algorithms::decomposeLU(mat);

    std::cout << "L = " << std::endl << lu.L << std::endl;

    std::cout << "U = " << std::endl << lu.U << std::endl;

    std::cout << "LU = " << std::endl << lu.L* lu.U << std::endl;

    std::cout << "===============================" << std::endl;
}

template<typename ElementType2>
void exampleHessenbergForm(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    std::cout << "===============================" << std::endl;
    std::cout << "HESSENBERG FORM EXAMPLE" << std::endl;

    std::cout << "M = " << std::endl << mat << std::endl;

    GeNuSys::LinAlg::HessenbergForm<ElementType2> hf = GeNuSys::LinAlg::Algorithms::getHessenbergForm(mat);

    std::cout << "A = " << std::endl << hf.A << std::endl;

    std::cout << "Q = " << std::endl << hf.Q << std::endl;

    std::cout << "QTQ = " << std::endl << hf.QT* hf.Q << std::endl;

    std::cout << "QTAQ = " << std::endl << hf.QT* hf.A* hf.Q << std::endl;

    std::cout << "===============================" << std::endl;
}

template<typename ElementType2>
void exampleSchurForm(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    typedef typename GeNuSys::ElementTraits<typename GeNuSys::ElementTraits<ElementType2>::ComplexType>::RealType ComplexRealType;

    std::cout << "===============================" << std::endl;
    std::cout << "SCHUR FORM EXAMPLE" << std::endl;

    std::cout << "M = " << std::endl << mat << std::endl;

    GeNuSys::LinAlg::SchurForm<ElementType2> schur = GeNuSys::LinAlg::Algorithms::getSchurForm(mat);

    std::cout << "QTMQ = " << std::endl << schur.QT* GeNuSys::LinAlg::Matrix<ComplexRealType>(mat) * schur.Q << std::endl;

    std::cout << "U = " << std::endl << schur.U << std::endl;

    std::cout << "Q = " << std::endl << schur.Q << std::endl;

    std::cout << "QTQ = " << std::endl << schur.QT* schur.Q << std::endl;

    std::cout << "QUQT = " << std::endl << schur.Q* schur.U* schur.QT << std::endl;

    std::cout << "DIFF = " << GeNuSys::LinAlg::PNorm<00>::norm(GeNuSys::LinAlg::Matrix<ComplexRealType>(mat) - schur.Q * schur.U * schur.QT) << std::endl;

    std::cout << "===============================" << std::endl;
}

template<typename ElementType2>
void exampleJordanForm(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    typedef typename GeNuSys::ElementTraits<typename GeNuSys::ElementTraits<ElementType2>::ComplexType>::RealType ComplexRealType;

    std::cout << "===============================" << std::endl;
    std::cout << "JORDAN FORM EXAMPLE" << std::endl;

    std::cout << "M = " << std::endl << mat << std::endl;

    GeNuSys::LinAlg::JordanForm<ElementType2> jordan = GeNuSys::LinAlg::Algorithms::getJordanForm(mat);

    std::cout << "PMinvP = " << std::endl << jordan.P* GeNuSys::LinAlg::Matrix<ComplexRealType>(mat) * jordan.invP << std::endl;

    std::cout << "P = " << std::endl << jordan.P << std::endl;

    std::cout << "J = " << std::endl << jordan.J << std::endl;

    std::cout << "invPJP = " << std::endl << jordan.invP* jordan.J* jordan.P << std::endl;

    std::cout << GeNuSys::LinAlg::PNorm<00>::norm(GeNuSys::LinAlg::Matrix<ComplexRealType>(mat) - jordan.invP * jordan.J * jordan.P) << std::endl;

    std::cout << "===============================" << std::endl;
}

template<typename ElementType2>
void exampleNumSys(const GeNuSys::LinAlg::Matrix<ElementType2>& mat)
{
    GeNuSys::LinAlg::Matrix<ElementType2> M = mat;

    GeNuSys::NumSys::RadixProperties<ElementType2> props(M);

    std::vector<GeNuSys::LinAlg::SparseVector<ElementType2>> jcds = GeNuSys::NumSys::DigitSet::getJCanonical(props, 0);
    std::cout << "JCanonical DigitSet:" << std::endl;
    for (unsigned int i = 0; i < jcds.size(); ++i)
    {
        std::cout << "Digit " << i << " : " << jcds[i] << std::endl;
    }

    std::vector<GeNuSys::LinAlg::SparseVector<ElementType2>> sds = GeNuSys::NumSys::DigitSet::getJSymmetric(props, 0);
    std::cout << "JSymmetric DigitSet:" << std::endl;
    for (unsigned int i = 0; i < sds.size(); ++i)
    {
        std::cout << "Digit " << i << " : " << sds[i] << std::endl;
    }

    std::vector<GeNuSys::LinAlg::Vector<ElementType2>> ads = GeNuSys::NumSys::DigitSet::getAdjoint(props);
    std::cout << "Adjoint DigitSet:" << std::endl;
    for (unsigned int i = 0; i < ads.size(); ++i)
    {
        std::cout << "Digit " << i << " : " << ads[i] << std::endl;
    }

    std::vector<GeNuSys::LinAlg::Vector<ElementType2>> dds = GeNuSys::NumSys::DigitSet::getDense(props);
    std::cout << "Dense DigitSet:" << std::endl;
    for (unsigned int i = 0; i < dds.size(); ++i)
    {
        std::cout << "Digit " << i << " : " << dds[i] << std::endl;
    }

    GeNuSys::NumSys::NumberSystem <
    ElementType2,
    GeNuSys::LinAlg::Vector,
    GeNuSys::LinAlg::SparseMatrix,
    GeNuSys::LinAlg::OperatorNorm<typename GeNuSys::ElementTraits<ElementType2>::RationalType>
    > numSys(props, dds, props.getOperatorNorm());
    std::cout << "NUMSYS CREATED!" << std::endl;
    std::cout << "Volume: " << GeNuSys::NumSys::Traits::getVolume(props.getInverse(), dds) << std::endl;
    clock_t tStart = clock();

    // GeNuSys::NumSys::BitVector bitVec(16000000000L);

    std::vector<std::vector<GeNuSys::LinAlg::Vector<ElementType2>>> cycles = numSys.getCycles();
    for (unsigned int i = 0; i < cycles.size(); ++i)
    {
        std::cout << "CYCLE " << i << " : " << std::endl;
        for (unsigned int j = 0; j < cycles[i].size(); ++j)
        {
            std::cout << "  " << cycles[i][j] << std::endl;
        }
    }

    printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

void exampleSimultaneous(const int i, const int j, bool impr)
{
    GeNuSys::LinAlg::Matrix<int> M = GeNuSys::NumSys::Simultaneous::getSimultaneous(i, j);

    std::cout << "===============================" << std::endl;

    std::cout << "SIMULTANEOUS (" << i << ", " << j << ")" << std::endl;

    std::cout << "M = " << M << std::endl;

    //exampleJordanForm(GeNuSys::LinAlg::Algorithms::invert(M));

    if (GeNuSys::LinAlg::Algorithms::det(M) == 0)
    {
        std::cout << "NOT INVERTIBLE" << std::endl;
        return;
    }

    GeNuSys::NumSys::RadixProperties<int> props(M);

    if (props.getOperatorNorm().norm(props.getInverse()) >= 1.0)
    {
        std::cout << "NOT EXPANSIVE" << std::endl;
        return;
    }


    std::vector<GeNuSys::LinAlg::Vector<int>> digits = GeNuSys::NumSys::Simultaneous::/*readDigitSet(); */getDigitSet(props,/* props.getOperatorNorm()*/GeNuSys::LinAlg::PNorm<2>(), true);

    GeNuSys::LinAlg::Matrix<int> imprM;
    std::vector<GeNuSys::LinAlg::Vector<int>> imprDigits;

    if (impr)
    {
        GeNuSys::LinAlg::Matrix<int> T = GeNuSys::NumSys::Traits::findBasisTransformation(props.getInverse(), digits, 15, 5, 2);

        imprM = T * M * GeNuSys::LinAlg::Traits::template convertUnsafe<double, int>(GeNuSys::LinAlg::Algorithms::invert(T));
        for (unsigned int i = 0; i < digits.size(); ++i)
        {
            imprDigits.push_back(T * digits[i]);
        }
    }
    else
    {
        imprM = M;
        for (unsigned int i = 0; i < digits.size(); ++i)
        {
            imprDigits.push_back(digits[i]);
        }
    }

    GeNuSys::NumSys::RadixProperties<int> imprProps(imprM);


    std::cout << "Digits = " << std::endl;
    std::cout << "[ " << digits[0][0];
    for (unsigned int i = 1; i < digits.size(); ++i)
    {
        std::cout << ", " << digits[i][0];
    }
    std::cout << " ]," << std::endl;

    std::cout << "[ " << digits[0][1];
    for (unsigned int i = 1; i < digits.size(); ++i)
    {
        std::cout << ", " << digits[i][1];
    }
    std::cout << " ]" << std::endl;

    std::vector<int> lowerBound, upperBound;
    GeNuSys::NumSys::Traits::getBounds(imprProps.getInverse(), imprDigits, lowerBound, upperBound);

    std::cout << "LOWER : [";
    for (unsigned int i = 0; i < lowerBound.size(); ++i)
    {
        std::cout << " " << lowerBound[i];
    }
    std::cout << " ]" <<  std::endl;

    std::cout << "UPPER : [";
    for (unsigned int i = 0; i < upperBound.size(); ++i)
    {
        std::cout << " " << upperBound[i];
    }
    std::cout << " ]" <<  std::endl;

    std::cout << "VOL = " <<  GeNuSys::NumSys::Traits::getVolume(imprProps.getInverse(), imprDigits) << std::endl;



    GeNuSys::NumSys::NumberSystem<int, GeNuSys::LinAlg::Vector, GeNuSys::LinAlg::SparseMatrix, GeNuSys::LinAlg::OperatorNorm<double>> numSys(imprProps, imprDigits, imprProps.getOperatorNorm());

    clock_t tStart = clock();

    std::vector<std::vector<GeNuSys::LinAlg::Vector<int>>> cycles = numSys.getCycles();
    std::cout << std::endl;
    for (unsigned int i = 0; i < cycles.size(); ++i)
    {
        std::cout << "CYCLE " << i << " : " << std::endl;
        for (unsigned int j = 0; j < cycles[i].size(); ++j)
        {
            std::cout << "  " << cycles[i][j] << std::endl;
        }
    }

    printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    std::cout << "===============================" << std::endl;
}

int main()
{
    srand(time(NULL));
    std::cout << std::fixed;
    std::cout.precision(6);


    GeNuSys::LinAlg::Matrix<int> smithMat(3, 3);
    smithMat.set(0, 0, 9);
    smithMat.set(0, 1, -36);
    smithMat.set(0, 2, 30);
    smithMat.set(1, 0, -36);
    smithMat.set(1, 1, 192);
    smithMat.set(1, 2, -180);
    smithMat.set(2, 0, 30);
    smithMat.set(2, 1, -180);
    smithMat.set(2, 2, 180);

    GeNuSys::LinAlg::Matrix<int> cnsMat(11, 11);
    for (unsigned int i = 1; i < 11; ++i)
    {
        cnsMat.set(i, i - 1, 1);
    }
    cnsMat.set(0, 10, -2);
    cnsMat.set(1, 10, 1);

    GeNuSys::LinAlg::Matrix<int> cnsMatSmall(4, 4);
    for (unsigned int i = 1; i < 4; ++i)
    {
        cnsMatSmall.set(i, i - 1, 1);
    }
    cnsMatSmall.set(0, 3, -2);
    cnsMatSmall.set(1, 3, 1);

    GeNuSys::LinAlg::Matrix<int> cnsMat2(8, 8);
    for (unsigned int i = 1; i < 8; ++i)
    {
        cnsMat2.set(i, i - 1, 1);
    }
    cnsMat2.set(0, 7, -2);
    cnsMat2.set(1, 7, 1);

    GeNuSys::LinAlg::Matrix<int> jordanMat(6, 6);
    jordanMat.set(0, 4, -1);
    jordanMat.set(0, 5, -1);

    jordanMat.set(1, 1, -8);
    jordanMat.set(1, 2, 4);
    jordanMat.set(1, 3, -3);
    jordanMat.set(1, 4, 1);
    jordanMat.set(1, 5, -3);

    jordanMat.set(2, 0, -3);
    jordanMat.set(2, 1, 13);
    jordanMat.set(2, 2, -8);
    jordanMat.set(2, 3, 6);
    jordanMat.set(2, 4, 2);
    jordanMat.set(2, 5, 9);

    jordanMat.set(3, 0, -2);
    jordanMat.set(3, 1, 14);
    jordanMat.set(3, 2, -7);
    jordanMat.set(3, 3, 4);
    jordanMat.set(3, 4, 2);
    jordanMat.set(3, 5, 10);

    jordanMat.set(4, 0, 1);
    jordanMat.set(4, 1, -18);
    jordanMat.set(4, 2, 11);
    jordanMat.set(4, 3, -11);
    jordanMat.set(4, 4, 2);
    jordanMat.set(4, 5, -6);

    jordanMat.set(5, 0, -1);
    jordanMat.set(5, 1, 19);
    jordanMat.set(5, 2, -11);
    jordanMat.set(5, 3, 10);
    jordanMat.set(5, 4, -2);
    jordanMat.set(5, 5, 7);

    exampleJordanForm(jordanMat);

    /*
    exampleSimultaneous(0, -1, false);
    exampleSimultaneous(-1, 0, false);
    exampleSimultaneous(-1, -1, false);
    */

    //int i, j;
    //std::cin >> i;
    //std::cin >> j;
    exampleNumSys(GeNuSys::NumSys::Simultaneous::getSimultaneous(2, 1));
    exampleSimultaneous(2, 1, false);
    exampleSimultaneous(2, 1, true);

    /*
    exampleSimultaneous(0, 1, false);
    exampleSimultaneous(1, 1, false);
    exampleSimultaneous(2, 0, false);
    exampleSimultaneous(2, 1, false);
    exampleSimultaneous(2, 2, false);
    exampleSimultaneous(1, 2, false);
    exampleSimultaneous(0, 2, false);

    for (int i = -3; i <= 3; ++i) {
        for (int j = -3; j <= 3; ++j) {
            exampleSimultaneous(i, j, false);
        }
    }
    */

    exampleNumSys(cnsMatSmall);

    exampleSmithNormalForm(smithMat);
    //exampleDecomposeQR<int>(randomMatrix(5, 5, 10));
    //exampleDecomposeLU<int>(randomMatrix(5, 5, 10));
    exampleHessenbergForm<int>(cnsMat);
    exampleSchurForm<int>(cnsMat);
    exampleSchurForm<double>(GeNuSys::LinAlg::Algorithms::invert(cnsMat));

    exampleSchurForm<long long>(jordanMat);

    //for (unsigned int i = 0; i < 10; ++i) {
    //    exampleSchurForm<long long>(randomMatrix(40, 40, 10));
    //}

    //GeNuSys::LinAlg::Matrix<double> invM = GeNuSys::LinAlg::Algorithms::invert(cnsMat);
    //GeNuSys::LinAlg::OperatorNorm<double> opNorm(invM);

    //std::cout << opNorm.norm(invM) << std::endl;

    /*
    GeNuSys::NumSys::SmithHash<int, GeNuSys::LinAlg::Matrix> hash(GeNuSys::LinAlg::Algorithms::getSmithNormalForm(cnsMat));

    std::cout << GeNuSys::LinAlg::OperatorNorm<double>::norm(GeNuSys::LinAlg::Algorithms::invert(cnsMat)) << std::endl;

    GeNuSys::LinAlg::Matrix<double> invM = GeNuSys::LinAlg::Algorithms::invert(cnsMat);
    */

    //std::cout << invM << std::endl;

    //std::cout << invM * GeNuSys::LinAlg::Matrix<double>(cnsMat) << std::endl;

    std::complex<double> a(2.0, 0.0);
    double dd = 2.00001;
    GeNuSys::ElementTraits<double>::asTypeUnsafe<int>(dd);

    return 0;
}
