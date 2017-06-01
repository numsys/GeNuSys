#include <iostream>
#include <vector>

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

int main()
{
    std::vector<std::pair<int, GeNuSys::LinAlg::Matrix<long long>>> mats;

    /* Generated with (in Sage):
    x = ZZ['x'].0
    for j in range(2,30):
        p=cyclotomic_polynomial(j)
        m=companion_matrix(p(x+euler_phi(j)))
    */

    //Matrices with very large values (not fitting in long long) are omitted.

    //mats.push_back(std::make_pair(2,GeNuSys::LinAlg::Matrix<long long>{1,1,std::vector<long long>{-2}}));
    mats.push_back(std::make_pair(3, GeNuSys::LinAlg::Matrix<long long> {2, 2, std::vector<long long>{0, -7, 1, -5}}));
    mats.push_back(std::make_pair(4, GeNuSys::LinAlg::Matrix<long long> {2, 2, std::vector<long long>{0, -5, 1, -4}}));
    mats.push_back(std::make_pair(5, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -341, 1, 0, 0, -313, 0, 1, 0, -109, 0, 0, 1, -17}}));
    mats.push_back(std::make_pair(6, GeNuSys::LinAlg::Matrix<long long> {2, 2, std::vector<long long>{0, -3, 1, -3}}));
    mats.push_back(std::make_pair(7, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -55987, 1, 0, 0, 0, 0, -54121, 0, 1, 0, 0, 0, -21835, 0, 0, 1, 0, 0, -4705, 0, 0, 0, 1, 0, -571, 0, 0, 0, 0, 1, -37}}));
    mats.push_back(std::make_pair(8, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -257, 1, 0, 0, -256, 0, 1, 0, -96, 0, 0, 1, -16}}));
    mats.push_back(std::make_pair(9, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -46873, 1, 0, 0, 0, 0, -46764, 0, 1, 0, 0, 0, -19458, 0, 0, 1, 0, 0, -4321, 0, 0, 0, 1, 0, -540, 0, 0, 0, 0, 1, -36}}));
    mats.push_back(std::make_pair(10, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -205, 1, 0, 0, -215, 0, 1, 0, -85, 0, 0, 1, -15}}));
    mats.push_back(std::make_pair(11, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -11111111111, 1, 0, 0, 0, 0, 0, 0, 0, 0, -10987654321, 0, 1, 0, 0, 0, 0, 0, 0, 0, -4890260631, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1289971041, 0, 0, 0, 1, 0, 0, 0, 0, 0, -223336551, 0, 0, 0, 0, 1, 0, 0, 0, 0, -26518161, 0, 0, 0, 0, 0, 1, 0, 0, 0, -2186871, 0, 0, 0, 0, 0, 0, 1, 0, 0, -123681, 0, 0, 0, 0, 0, 0, 0, 1, 0, -4591, 0, 0, 0, 0, 0, 0, 0, 0, 1, -101}}));
    mats.push_back(std::make_pair(12, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -241, 1, 0, 0, -248, 0, 1, 0, -95, 0, 0, 1, -16}}));
    mats.push_back(std::make_pair(13, GeNuSys::LinAlg::Matrix<long long> {12, 12, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -9726655034461, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -9652968253897, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4391062241797, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1210663993297, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -225325359853, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -29823734809, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -2878513429, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -204130465, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -10556029, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -388201, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -9637, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -145}}));
    mats.push_back(std::make_pair(14, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -39991, 1, 0, 0, 0, 0, -40943, 0, 1, 0, 0, 0, -17479, 0, 0, 1, 0, 0, -3983, 0, 0, 0, 1, 0, -511, 0, 0, 0, 0, 1, -35}}));
    mats.push_back(std::make_pair(15, GeNuSys::LinAlg::Matrix<long long> {8, 8, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, -14709241, 1, 0, 0, 0, 0, 0, 0, -14960831, 0, 1, 0, 0, 0, 0, 0, -6656664, 0, 0, 1, 0, 0, 0, 0, -1692257, 0, 0, 0, 1, 0, 0, 0, -268839, 0, 0, 0, 0, 1, 0, 0, -27329, 0, 0, 0, 0, 0, 1, 0, -1736, 0, 0, 0, 0, 0, 0, 1, -63}}));
    mats.push_back(std::make_pair(16, GeNuSys::LinAlg::Matrix<long long> {8, 8, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, -16777217, 1, 0, 0, 0, 0, 0, 0, -16777216, 0, 1, 0, 0, 0, 0, 0, -7340032, 0, 0, 1, 0, 0, 0, 0, -1835008, 0, 0, 0, 1, 0, 0, 0, -286720, 0, 0, 0, 0, 1, 0, 0, -28672, 0, 0, 0, 0, 0, 1, 0, -1792, 0, 0, 0, 0, 0, 0, 1, -64}}));
    mats.push_back(std::make_pair(18, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -46441, 1, 0, 0, 0, 0, -46548, 0, 1, 0, 0, 0, -19422, 0, 0, 1, 0, 0, -4319, 0, 0, 0, 1, 0, -540, 0, 0, 0, 0, 1, -36}}));
    mats.push_back(std::make_pair(20, GeNuSys::LinAlg::Matrix<long long> {8, 8, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, -16519105, 1, 0, 0, 0, 0, 0, 0, -16582640, 0, 1, 0, 0, 0, 0, 0, -7278975, 0, 0, 1, 0, 0, 0, 0, -1824800, 0, 0, 0, 1, 0, 0, 0, -285761, 0, 0, 0, 0, 1, 0, 0, -28624, 0, 0, 0, 0, 0, 1, 0, -1791, 0, 0, 0, 0, 0, 0, 1, -64}}));
    mats.push_back(std::make_pair(21, GeNuSys::LinAlg::Matrix<long long> {12, 12, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8177824843189, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8238594109103, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3803964767172, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1064441620177, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -201046346351, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -27001783368, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -2644229953, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -190237152, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -9979307, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -372241, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -9372, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -143}}));
    mats.push_back(std::make_pair(22, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -9090909091, 1, 0, 0, 0, 0, 0, 0, 0, 0, -9173553719, 0, 1, 0, 0, 0, 0, 0, 0, 0, -4166040571, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1121269039, 0, 0, 0, 1, 0, 0, 0, 0, 0, -198066451, 0, 0, 0, 0, 1, 0, 0, 0, 0, -23993959, 0, 0, 0, 0, 0, 1, 0, 0, 0, -2018731, 0, 0, 0, 0, 0, 0, 1, 0, 0, -116479, 0, 0, 0, 0, 0, 0, 0, 1, 0, -4411, 0, 0, 0, 0, 0, 0, 0, 0, 1, -99}}));
    mats.push_back(std::make_pair(24, GeNuSys::LinAlg::Matrix<long long> {8, 8, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, -16773121, 1, 0, 0, 0, 0, 0, 0, -16775168, 0, 1, 0, 0, 0, 0, 0, -7339648, 0, 0, 1, 0, 0, 0, 0, -1834976, 0, 0, 0, 1, 0, 0, 0, -286719, 0, 0, 0, 0, 1, 0, 0, -28672, 0, 0, 0, 0, 0, 1, 0, -1792, 0, 0, 0, 0, 0, 0, 1, -64}}));
    mats.push_back(std::make_pair(26, GeNuSys::LinAlg::Matrix<long long> {12, 12, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8230246567621, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8283004558439, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3820896027325, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1068266933903, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -201613539829, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -27059454071, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -2648302189, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -190434335, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -9985573, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -372359, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -9373, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -143}}));
    mats.push_back(std::make_pair(28, GeNuSys::LinAlg::Matrix<long long> {12, 12, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8854610100337, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -8864787813096, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4067280159839, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1130865760560, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -212215332241, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -28316182968, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -2754698687, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -196867680, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -10257841, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -380040, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -9503, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -144}}));

    /* Generated with (in Sage):
    x = ZZ['x'].0
    for j in range(2,13):
        p=cyclotomic_polynomial(j)
        for k in range(2,euler_phi(j)):
            m=companion_matrix(p(x+k))
    */

    mats.push_back(std::make_pair(502, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -31, 1, 0, 0, -49, 0, 1, 0, -31, 0, 0, 1, -9}}));
    mats.push_back(std::make_pair(503, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -121, 1, 0, 0, -142, 0, 1, 0, -64, 0, 0, 1, -13}}));
    mats.push_back(std::make_pair(702, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -127, 1, 0, 0, 0, 0, -321, 0, 1, 0, 0, 0, -351, 0, 0, 1, 0, 0, -209, 0, 0, 0, 1, 0, -71, 0, 0, 0, 0, 1, -13}}));
    mats.push_back(std::make_pair(703, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -1093, 1, 0, 0, 0, 0, -2005, 0, 1, 0, 0, 0, -1549, 0, 0, 1, 0, 0, -643, 0, 0, 0, 1, 0, -151, 0, 0, 0, 0, 1, -19}}));
    mats.push_back(std::make_pair(704, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -5461, 1, 0, 0, 0, 0, -7737, 0, 1, 0, 0, 0, -4589, 0, 0, 1, 0, 0, -1457, 0, 0, 0, 1, 0, -261, 0, 0, 0, 0, 1, -25}}));
    mats.push_back(std::make_pair(705, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -19531, 1, 0, 0, 0, 0, -22461, 0, 1, 0, 0, 0, -10791, 0, 0, 1, 0, 0, -2771, 0, 0, 0, 1, 0, -401, 0, 0, 0, 0, 1, -31}}));
    mats.push_back(std::make_pair(802, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -17, 1, 0, 0, -32, 0, 1, 0, -24, 0, 0, 1, -8}}));
    mats.push_back(std::make_pair(803, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -82, 1, 0, 0, -108, 0, 1, 0, -54, 0, 0, 1, -12}}));
    mats.push_back(std::make_pair(902, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -73, 1, 0, 0, 0, 0, -204, 0, 1, 0, 0, 0, -246, 0, 0, 1, 0, 0, -161, 0, 0, 0, 1, 0, -60, 0, 0, 0, 0, 1, -12}}));
    mats.push_back(std::make_pair(903, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -757, 1, 0, 0, 0, 0, -1485, 0, 1, 0, 0, 0, -1224, 0, 0, 1, 0, 0, -541, 0, 0, 0, 1, 0, -135, 0, 0, 0, 0, 1, -18}}));
    mats.push_back(std::make_pair(904, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -4161, 1, 0, 0, 0, 0, -6192, 0, 1, 0, 0, 0, -3852, 0, 0, 1, 0, 0, -1281, 0, 0, 0, 1, 0, -240, 0, 0, 0, 0, 1, -24}}));
    mats.push_back(std::make_pair(905, GeNuSys::LinAlg::Matrix<long long> {6, 6, std::vector<long long>{0, 0, 0, 0, 0, -15751, 1, 0, 0, 0, 0, -18825, 0, 1, 0, 0, 0, -9390, 0, 0, 1, 0, 0, -2501, 0, 0, 0, 1, 0, -375, 0, 0, 0, 0, 1, -30}}));
    mats.push_back(std::make_pair(1002, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -11, 1, 0, 0, -23, 0, 1, 0, -19, 0, 0, 1, -7}}));
    mats.push_back(std::make_pair(1003, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -61, 1, 0, 0, -86, 0, 1, 0, -46, 0, 0, 1, -11}}));
    mats.push_back(std::make_pair(1102, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -2047, 1, 0, 0, 0, 0, 0, 0, 0, 0, -9217, 0, 1, 0, 0, 0, 0, 0, 0, 0, -18943, 0, 0, 1, 0, 0, 0, 0, 0, 0, -23297, 0, 0, 0, 1, 0, 0, 0, 0, 0, -18943, 0, 0, 0, 0, 1, 0, 0, 0, 0, -10625, 0, 0, 0, 0, 0, 1, 0, 0, 0, -4159, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1121, 0, 0, 0, 0, 0, 0, 0, 1, 0, -199, 0, 0, 0, 0, 0, 0, 0, 0, 1, -21}}));
    mats.push_back(std::make_pair(1103, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -88573, 1, 0, 0, 0, 0, 0, 0, 0, 0, -280483, 0, 1, 0, 0, 0, 0, 0, 0, 0, -401041, 0, 0, 1, 0, 0, 0, 0, 0, 0, -340762, 0, 0, 0, 1, 0, 0, 0, 0, 0, -190474, 0, 0, 0, 0, 1, 0, 0, 0, 0, -73162, 0, 0, 0, 0, 0, 1, 0, 0, 0, -19552, 0, 0, 0, 0, 0, 0, 1, 0, 0, -3589, 0, 0, 0, 0, 0, 0, 0, 1, 0, -433, 0, 0, 0, 0, 0, 0, 0, 0, 1, -31}}));
    mats.push_back(std::make_pair(1104, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -1398101, 1, 0, 0, 0, 0, 0, 0, 0, 0, -3378745, 0, 1, 0, 0, 0, 0, 0, 0, 0, -3679725, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2377905, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1009605, 0, 0, 0, 0, 1, 0, 0, 0, 0, -294249, 0, 0, 0, 0, 0, 1, 0, 0, 0, -59613, 0, 0, 0, 0, 0, 0, 1, 0, 0, -8289, 0, 0, 0, 0, 0, 0, 0, 1, 0, -757, 0, 0, 0, 0, 0, 0, 0, 0, 1, -41}}));
    mats.push_back(std::make_pair(1105, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -12207031, 1, 0, 0, 0, 0, 0, 0, 0, 0, -23803711, 0, 1, 0, 0, 0, 0, 0, 0, 0, -20904541, 0, 0, 1, 0, 0, 0, 0, 0, 0, -10887146, 0, 0, 0, 1, 0, 0, 0, 0, 0, -3723526, 0, 0, 0, 0, 1, 0, 0, 0, 0, -873806, 0, 0, 0, 0, 0, 1, 0, 0, 0, -142486, 0, 0, 0, 0, 0, 0, 1, 0, 0, -15941, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1171, 0, 0, 0, 0, 0, 0, 0, 0, 1, -51}}));
    mats.push_back(std::make_pair(1106, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -72559411, 1, 0, 0, 0, 0, 0, 0, 0, 0, -118513705, 0, 1, 0, 0, 0, 0, 0, 0, 0, -87151915, 0, 0, 1, 0, 0, 0, 0, 0, 0, -37996945, 0, 0, 0, 1, 0, 0, 0, 0, 0, -10876387, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2135737, 0, 0, 0, 0, 0, 1, 0, 0, 0, -291355, 0, 0, 0, 0, 0, 0, 1, 0, 0, -27265, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1675, 0, 0, 0, 0, 0, 0, 0, 0, 1, -61}}));
    mats.push_back(std::make_pair(1107, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -329554457, 1, 0, 0, 0, 0, 0, 0, 0, 0, -462945547, 0, 1, 0, 0, 0, 0, 0, 0, 0, -292750473, 0, 0, 1, 0, 0, 0, 0, 0, 0, -109740282, 0, 0, 0, 1, 0, 0, 0, 0, 0, -27004818, 0, 0, 0, 0, 1, 0, 0, 0, 0, -4558170, 0, 0, 0, 0, 0, 1, 0, 0, 0, -534444, 0, 0, 0, 0, 0, 0, 1, 0, 0, -42981, 0, 0, 0, 0, 0, 0, 0, 1, 0, -2269, 0, 0, 0, 0, 0, 0, 0, 0, 1, -71}}));
    mats.push_back(std::make_pair(1108, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -1227133513, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1512003793, 0, 1, 0, 0, 0, 0, 0, 0, 0, -838567321, 0, 0, 1, 0, 0, 0, 0, 0, 0, -275667617, 0, 0, 0, 1, 0, 0, 0, 0, 0, -59484649, 0, 0, 0, 0, 1, 0, 0, 0, 0, -8803697, 0, 0, 0, 0, 0, 1, 0, 0, 0, -905017, 0, 0, 0, 0, 0, 0, 1, 0, 0, -63809, 0, 0, 0, 0, 0, 0, 0, 1, 0, -2953, 0, 0, 0, 0, 0, 0, 0, 0, 1, -81}}));
    mats.push_back(std::make_pair(1109, GeNuSys::LinAlg::Matrix<long long> {10, 10, std::vector<long long>{0, 0, 0, 0, 0, 0, 0, 0, 0, -3922632451, 1, 0, 0, 0, 0, 0, 0, 0, 0, -4303999495, 0, 1, 0, 0, 0, 0, 0, 0, 0, -2125515925, 0, 0, 1, 0, 0, 0, 0, 0, 0, -622149130, 0, 0, 0, 1, 0, 0, 0, 0, 0, -119528830, 0, 0, 0, 0, 1, 0, 0, 0, 0, -15749614, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1441378, 0, 0, 0, 0, 0, 0, 1, 0, 0, -90469, 0, 0, 0, 0, 0, 0, 0, 1, 0, -3727, 0, 0, 0, 0, 0, 0, 0, 0, 1, -91}}));
    mats.push_back(std::make_pair(1202, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -13, 1, 0, 0, -28, 0, 1, 0, -23, 0, 0, 1, -8}}));
    mats.push_back(std::make_pair(1203, GeNuSys::LinAlg::Matrix<long long> {4, 4, std::vector<long long>{0, 0, 0, -73, 1, 0, 0, -102, 0, 1, 0, -53, 0, 0, 1, -12}}));
    for (const auto& mat : mats)
    {
        std::cout << "=======================================" << std::endl;
        GeNuSys::NumSys::RadixProperties<long long> props(mat.second);
        std::cout << "n=" << mat.first << std::endl;
        std::cout << "Abs det:" << props.getAbsDet() << std::endl;
        if (props.getAbsDet() < 100000)
        {
            std::cout << "Computing digit set..." << std::endl;
            std::vector<GeNuSys::LinAlg::SparseVector<long long>> sds = GeNuSys::NumSys::DigitSet::getJSymmetric(props, 0);
            if (GeNuSys::NumSys::Traits::getVolume(props.getInverse(), sds) > 1000000000)
            {
                //Skip if volume too large so this example can run in a few minutes
                std::cout << "Large volume, skipped" << std::endl;
                continue;
            }
            std::cout << "Creating number system object..." << std::endl;
            std::cout << "Volume:" << GeNuSys::NumSys::Traits::getVolume(props.getInverse(), sds) << std::endl;
            GeNuSys::LinAlg::Matrix<long long> T = GeNuSys::NumSys::Traits::findBasisTransformation(props.getInverse(), sds, 15, 5, 2);
            auto imprM = T * mat.second * GeNuSys::LinAlg::Traits::template convertUnsafe<typename GeNuSys::ElementTraits<long long int>::RationalType, long long>(GeNuSys::LinAlg::Algorithms::invert(T));
            std::vector<GeNuSys::LinAlg::SparseVector<long long>> imprDigits;
            for (unsigned int i = 0; i < sds.size(); ++i)
            {
                imprDigits.push_back(T * sds[i]);
            }
            GeNuSys::NumSys::RadixProperties<long long> imprProps(imprM);
            std::cout << "Reduced volume:" << GeNuSys::NumSys::Traits::getVolume(imprProps.getInverse(), imprDigits) << std::endl;
            GeNuSys::NumSys::NumberSystem<long long, GeNuSys::LinAlg::SparseVector, GeNuSys::LinAlg::Matrix, GeNuSys::LinAlg::OperatorNorm<typename GeNuSys::ElementTraits<long long>::RationalType>> numSys(imprProps, imprDigits, imprProps.getOperatorNorm());
            auto cycles = numSys.getCycles();
            for (unsigned int i = 0; i < cycles.size(); ++i)
            {
                std::cout << "CYCLE " << i << " : " << std::endl;
                for (unsigned int j = 0; j < cycles[i].size(); ++j)
                {
                    std::cout << "  " << cycles[i][j] << std::endl;
                }
            }
        }
        else
        {
            //Skip if digit set too large so this example can run in a few minutes
            std::cout << "Large digit set, skipped" << std::endl;
        }
    }

    return 0;
}
