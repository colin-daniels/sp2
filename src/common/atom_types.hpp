#ifndef SP2_ATOM_TYPES_HPP
#define SP2_ATOM_TYPES_HPP

#include <cstdint>
#include <type_traits>

namespace sp2 {

namespace detail {
// "magic" ids that ensure that (id[i] ^ id[j]) is unique for all i and j
constexpr std::uint32_t atom_ids[] = {
        0,     1,     2,     3,     4,     8,    12,    16,    32,    48,    64,
       85,   106,   128,   192,   256,   278,   329,   356,   401,   512,   537,
      594,   645,   768,  1024,  1050,  1094,  1121,  1161,  1412,  1556,  1666,
     1952,  2048,  2085,  2136,  2274,  2315,  2600,  2628,  3072,  3680,  4096,
     4134,  4234,  4244,  4385,  4552,  4673,  4988,  5160,  5392,  5645,  6304,
     6402,  7258,  8192,  8233,  8338,  8519,  8710,  8909,  9004,  9296,  9482,
     9761, 10257, 10343, 10428, 10624, 11557, 11914, 11993, 12288, 12366, 13346,
    14990, 14993, 16384, 16518, 16552, 16613, 16645, 16930, 16968, 17872, 18450,
    18754, 18840, 19571, 19660, 21136, 21505, 22211, 22536, 22614, 24644, 24840,
    26143, 28022, 28161, 28721, 29496, 29700, 31216, 31405, 31490, 32768, 32922,
    32929, 33048, 33290, 33797, 34624, 34836, 35016, 35632, 36288, 36742, 37857,
    38299, 38400, 38918, 39476, 39758, 41056, 41085, 41218, 41553, 43959, 45065,
    45118, 45232, 48279, 49152, 50227, 50320, 50720, 50962, 51333, 52387, 53572,
    54372, 55773, 58376, 59424, 59866, 60259, 62221, 62528, 64799
};
} // namespace detail

using namespace detail;

enum class bond_types : std::uint32_t;
enum class atom_types : std::uint32_t
{
    NONE = 0,
    // periodic table
    H  =   atom_ids[1],    He =   atom_ids[2],    Li =   atom_ids[3],
    Be =   atom_ids[4],    B  =   atom_ids[5],    C  =   atom_ids[6],
    N  =   atom_ids[7],    O  =   atom_ids[8],    F  =   atom_ids[9],
    Ne =  atom_ids[10],    Na =  atom_ids[11],    Mg =  atom_ids[12],
    Al =  atom_ids[13],    Si =  atom_ids[14],    P  =  atom_ids[15],
    S  =  atom_ids[16],    Cl =  atom_ids[17],    Ar =  atom_ids[18],
    K  =  atom_ids[19],    Ca =  atom_ids[20],    Sc =  atom_ids[21],
    Ti =  atom_ids[22],    V  =  atom_ids[23],    Cr =  atom_ids[24],
    Mn =  atom_ids[25],    Fe =  atom_ids[26],    Co =  atom_ids[27],
    Ni =  atom_ids[28],    Cu =  atom_ids[29],    Zn =  atom_ids[30],
    Ga =  atom_ids[31],    Ge =  atom_ids[32],    As =  atom_ids[33],
    Se =  atom_ids[34],    Br =  atom_ids[35],    Kr =  atom_ids[36],
    Rb =  atom_ids[37],    Sr =  atom_ids[38],    Y  =  atom_ids[39],
    Zr =  atom_ids[40],    Nb =  atom_ids[41],    Mo =  atom_ids[42],
    Tc =  atom_ids[43],    Ru =  atom_ids[44],    Rh =  atom_ids[45],
    Pd =  atom_ids[46],    Ag =  atom_ids[47],    Cd =  atom_ids[48],
    In =  atom_ids[49],    Sn =  atom_ids[50],    Sb =  atom_ids[51],
    Te =  atom_ids[52],    I  =  atom_ids[53],    Xe =  atom_ids[54],
    Cs =  atom_ids[55],    Ba =  atom_ids[56],    La =  atom_ids[57],
    Ce =  atom_ids[58],    Pr =  atom_ids[59],    Nd =  atom_ids[60],
    Pm =  atom_ids[61],    Sm =  atom_ids[62],    Eu =  atom_ids[63],
    Gd =  atom_ids[64],    Tb =  atom_ids[65],    Dy =  atom_ids[66],
    Ho =  atom_ids[67],    Er =  atom_ids[68],    Tm =  atom_ids[69],
    Yb =  atom_ids[70],    Lu =  atom_ids[71],    Hf =  atom_ids[72],
    Ta =  atom_ids[73],    W  =  atom_ids[74],    Re =  atom_ids[75],
    Os =  atom_ids[76],    Ir =  atom_ids[77],    Pt =  atom_ids[78],
    Au =  atom_ids[79],    Hg =  atom_ids[80],    Tl =  atom_ids[81],
    Pb =  atom_ids[82],    Bi =  atom_ids[83],    Po =  atom_ids[84],
    At =  atom_ids[85],    Rn =  atom_ids[86],    Fr =  atom_ids[87],
    Ra =  atom_ids[88],    Ac =  atom_ids[89],    Th =  atom_ids[90],
    Pa =  atom_ids[91],    U  =  atom_ids[92],    Np =  atom_ids[93],
    Pu =  atom_ids[94],    Am =  atom_ids[95],    Cm =  atom_ids[96],
    Bk =  atom_ids[97],    Cf =  atom_ids[98],    Es =  atom_ids[99],
    Fm = atom_ids[100],    Md = atom_ids[101],    No = atom_ids[102],
    Lr = atom_ids[103],    Rf = atom_ids[104],    Db = atom_ids[105],
    Sg = atom_ids[106],    Bh = atom_ids[107],    Hs = atom_ids[108],
    Mt = atom_ids[109],    Ds = atom_ids[110],    Rg = atom_ids[111],
    Cn = atom_ids[112],    Nh = atom_ids[113],    Fl = atom_ids[114],
    Mc = atom_ids[115],    Lv = atom_ids[116],    Ts = atom_ids[117],
    Og = atom_ids[118]
};

constexpr bond_types btype(atom_types a, atom_types b)
{
    using ut = std::underlying_type_t<atom_types>;
    return static_cast<bond_types>(
        static_cast<ut>(a) ^ static_cast<ut>(b)
    );
}

using atype = atom_types;

} // namespace sp2

#endif // SP2_ATOM_TYPES_HPP
