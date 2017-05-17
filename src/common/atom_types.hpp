#ifndef SP2_ATOM_TYPES_HPP
#define SP2_ATOM_TYPES_HPP

#include "common/enums.hpp"

#include <cstdint>
#include <type_traits>

namespace sp2 {

/// atom type enum, uses special ids to ensure (atom_type::A ^ atom_type::B)
/// is unique, which is how bond types are calculated
enum class atom_types : std::uint32_t;

/// convenience using-declaration for referencing atom types
using atype = atom_types;

/// get atomic number from input atom type
constexpr unsigned int atomic_number(atom_types type);


/// (mostly) opaque bond type enum, underlying value calculated by
///     sp2::btype(atom_type, atom_type)
/// only commonly used bond types are explicitly listed due to 118^2 being a
/// lot of enums to have in source
enum class bond_types : std::uint32_t;

/// determine bond type from two atom types, order does not matter
constexpr bond_types btype(atom_types a, atom_types b);

////////////////////////////////////////////////////////////////////////////////
// Implementation details                                                     //
////////////////////////////////////////////////////////////////////////////////

namespace detail {
/// ids such that for all i and j, (id[i] ^ id[j]) is unique && id[i] < 2^16,
/// do not depend on these values being constant throughout versions
constexpr std::uint16_t xor_ids[] = {
        0,     1,     2,     4,     8,    15,    16,    32,    51,    64,
       85,   106,   128,   150,   171,   219,   237,   247,   256,   279,
      297,   455,   512,   537,   557,   594,   643,   803,   863,   998,
     1024,  1051,  1070,  1112,  1169,  1333,  1345,  1620,  1866,  2048,
     2076,  2085,  2185,  2372,  2456,  2618,  2800,  2873,  3127,  3284,
     3483,  3557,  3763,  4096,  4125,  4135,  4174,  4435,  4459,  4469,
     4497,  4752,  5255,  5732,  5804,  5915,  6100,  6369,  6907,  7069,
     8192,  8263,  8351,  8422,  8458,  8571,  8750,  8858,  9124,  9314,
     9500, 10026, 10455, 10556, 11778, 11885, 11984, 13548, 14007, 14514,
    14965, 15125, 15554, 16384, 16457, 16517, 16609, 16771, 16853, 17022,
    17453, 17891, 18073, 18562, 18980, 19030, 19932, 20075, 20745, 21544,
    22633, 23200, 24167, 25700, 26360, 26591, 26776, 28443, 28905, 29577,
    32705, 32768, 32844, 32949, 33042, 33251, 33354, 33582, 33930, 34082,
    34583, 34933, 35379, 36604, 36686, 37304, 37492, 39816, 39894, 39936,
    40542, 42236, 42731, 43283, 45004, 45762, 46146, 48698, 50969, 51520,
    54582, 55959, 56663, 57422, 60693, 61593, 62498
};
} // namespace detail

/// somewhat pointless to have xor_ids in detail given we immediately bring it
/// back into this namespace but it's just to highlight that xor_ids shouldn't
/// be used elsewhere
using detail::xor_ids;

enum class atom_types : std::uint32_t
{
    NONE = 0,
    // periodic table
    H  =   xor_ids[1], He =   xor_ids[2], Li =   xor_ids[3],
    Be =   xor_ids[4], B  =   xor_ids[5], C  =   xor_ids[6],
    N  =   xor_ids[7], O  =   xor_ids[8], F  =   xor_ids[9],
    Ne =  xor_ids[10], Na =  xor_ids[11], Mg =  xor_ids[12],
    Al =  xor_ids[13], Si =  xor_ids[14], P  =  xor_ids[15],
    S  =  xor_ids[16], Cl =  xor_ids[17], Ar =  xor_ids[18],
    K  =  xor_ids[19], Ca =  xor_ids[20], Sc =  xor_ids[21],
    Ti =  xor_ids[22], V  =  xor_ids[23], Cr =  xor_ids[24],
    Mn =  xor_ids[25], Fe =  xor_ids[26], Co =  xor_ids[27],
    Ni =  xor_ids[28], Cu =  xor_ids[29], Zn =  xor_ids[30],
    Ga =  xor_ids[31], Ge =  xor_ids[32], As =  xor_ids[33],
    Se =  xor_ids[34], Br =  xor_ids[35], Kr =  xor_ids[36],
    Rb =  xor_ids[37], Sr =  xor_ids[38], Y  =  xor_ids[39],
    Zr =  xor_ids[40], Nb =  xor_ids[41], Mo =  xor_ids[42],
    Tc =  xor_ids[43], Ru =  xor_ids[44], Rh =  xor_ids[45],
    Pd =  xor_ids[46], Ag =  xor_ids[47], Cd =  xor_ids[48],
    In =  xor_ids[49], Sn =  xor_ids[50], Sb =  xor_ids[51],
    Te =  xor_ids[52], I  =  xor_ids[53], Xe =  xor_ids[54],
    Cs =  xor_ids[55], Ba =  xor_ids[56], La =  xor_ids[57],
    Ce =  xor_ids[58], Pr =  xor_ids[59], Nd =  xor_ids[60],
    Pm =  xor_ids[61], Sm =  xor_ids[62], Eu =  xor_ids[63],
    Gd =  xor_ids[64], Tb =  xor_ids[65], Dy =  xor_ids[66],
    Ho =  xor_ids[67], Er =  xor_ids[68], Tm =  xor_ids[69],
    Yb =  xor_ids[70], Lu =  xor_ids[71], Hf =  xor_ids[72],
    Ta =  xor_ids[73], W  =  xor_ids[74], Re =  xor_ids[75],
    Os =  xor_ids[76], Ir =  xor_ids[77], Pt =  xor_ids[78],
    Au =  xor_ids[79], Hg =  xor_ids[80], Tl =  xor_ids[81],
    Pb =  xor_ids[82], Bi =  xor_ids[83], Po =  xor_ids[84],
    At =  xor_ids[85], Rn =  xor_ids[86], Fr =  xor_ids[87],
    Ra =  xor_ids[88], Ac =  xor_ids[89], Th =  xor_ids[90],
    Pa =  xor_ids[91], U  =  xor_ids[92], Np =  xor_ids[93],
    Pu =  xor_ids[94], Am =  xor_ids[95], Cm =  xor_ids[96],
    Bk =  xor_ids[97], Cf =  xor_ids[98], Es =  xor_ids[99],
    Fm = xor_ids[100], Md = xor_ids[101], No = xor_ids[102],
    Lr = xor_ids[103], Rf = xor_ids[104], Db = xor_ids[105],
    Sg = xor_ids[106], Bh = xor_ids[107], Hs = xor_ids[108],
    Mt = xor_ids[109], Ds = xor_ids[110], Rg = xor_ids[111],
    Cn = xor_ids[112], Nh = xor_ids[113], Fl = xor_ids[114],
    Mc = xor_ids[115], Lv = xor_ids[116], Ts = xor_ids[117],
    Og = xor_ids[118]
};


/// specialization of enum_map<> template variable to enable to/from string
/// conversions for atom_types
template<>
constexpr enum_map_t<atom_types> enum_map<atom_types> = {
    {atype::NONE, "NONE"},
    // periodic table
    {atype::H,  "H"},  {atype::He, "He"}, {atype::Li, "Li"}, {atype::Be, "Be"},
    {atype::B,  "B"},  {atype::C,  "C"},  {atype::N,  "N"},  {atype::O,  "O"},
    {atype::F,  "F"},  {atype::Ne, "Ne"}, {atype::Na, "Na"}, {atype::Mg, "Mg"},
    {atype::Al, "Al"}, {atype::Si, "Si"}, {atype::P,  "P"},  {atype::S,  "S"},
    {atype::Cl, "Cl"}, {atype::Ar, "Ar"}, {atype::K,  "K"},  {atype::Ca, "Ca"},
    {atype::Sc, "Sc"}, {atype::Ti, "Ti"}, {atype::V,  "V"},  {atype::Cr, "Cr"},
    {atype::Mn, "Mn"}, {atype::Fe, "Fe"}, {atype::Co, "Co"}, {atype::Ni, "Ni"},
    {atype::Cu, "Cu"}, {atype::Zn, "Zn"}, {atype::Ga, "Ga"}, {atype::Ge, "Ge"},
    {atype::As, "As"}, {atype::Se, "Se"}, {atype::Br, "Br"}, {atype::Kr, "Kr"},
    {atype::Rb, "Rb"}, {atype::Sr, "Sr"}, {atype::Y,  "Y"},  {atype::Zr, "Zr"},
    {atype::Nb, "Nb"}, {atype::Mo, "Mo"}, {atype::Tc, "Tc"}, {atype::Ru, "Ru"},
    {atype::Rh, "Rh"}, {atype::Pd, "Pd"}, {atype::Ag, "Ag"}, {atype::Cd, "Cd"},
    {atype::In, "In"}, {atype::Sn, "Sn"}, {atype::Sb, "Sb"}, {atype::Te, "Te"},
    {atype::I,  "I"},  {atype::Xe, "Xe"}, {atype::Cs, "Cs"}, {atype::Ba, "Ba"},
    {atype::La, "La"}, {atype::Ce, "Ce"}, {atype::Pr, "Pr"}, {atype::Nd, "Nd"},
    {atype::Pm, "Pm"}, {atype::Sm, "Sm"}, {atype::Eu, "Eu"}, {atype::Gd, "Gd"},
    {atype::Tb, "Tb"}, {atype::Dy, "Dy"}, {atype::Ho, "Ho"}, {atype::Er, "Er"},
    {atype::Tm, "Tm"}, {atype::Yb, "Yb"}, {atype::Lu, "Lu"}, {atype::Hf, "Hf"},
    {atype::Ta, "Ta"}, {atype::W,  "W"},  {atype::Re, "Re"}, {atype::Os, "Os"},
    {atype::Ir, "Ir"}, {atype::Pt, "Pt"}, {atype::Au, "Au"}, {atype::Hg, "Hg"},
    {atype::Tl, "Tl"}, {atype::Pb, "Pb"}, {atype::Bi, "Bi"}, {atype::Po, "Po"},
    {atype::At, "At"}, {atype::Rn, "Rn"}, {atype::Fr, "Fr"}, {atype::Ra, "Ra"},
    {atype::Ac, "Ac"}, {atype::Th, "Th"}, {atype::Pa, "Pa"}, {atype::U,  "U"},
    {atype::Np, "Np"}, {atype::Pu, "Pu"}, {atype::Am, "Am"}, {atype::Cm, "Cm"},
    {atype::Bk, "Bk"}, {atype::Cf, "Cf"}, {atype::Es, "Es"}, {atype::Fm, "Fm"},
    {atype::Md, "Md"}, {atype::No, "No"}, {atype::Lr, "Lr"}, {atype::Rf, "Rf"},
    {atype::Db, "Db"}, {atype::Sg, "Sg"}, {atype::Bh, "Bh"}, {atype::Hs, "Hs"},
    {atype::Mt, "Mt"}, {atype::Ds, "Ds"}, {atype::Rg, "Rg"}, {atype::Cn, "Cn"},
    {atype::Nh, "Nh"}, {atype::Fl, "Fl"}, {atype::Mc, "Mc"}, {atype::Lv, "Lv"},
    {atype::Ts, "Ts"}, {atype::Og, "Og"}
};

constexpr unsigned int atomic_number(atom_types type)
{
    using ut = std::underlying_type_t<atom_types>;
    auto xor_id = static_cast<ut>(type);

    std::size_t low = 0,
        high = std::extent<decltype(detail::xor_ids)>::value;

    while (low < high)
    {
        auto mid = (high + low) / 2;
        if (xor_id < detail::xor_ids[mid])
            high = mid;
        else if (xor_id > detail::xor_ids[mid])
            low = mid + 1;
        else
            return static_cast<unsigned int>(mid);
    }

    return 0;
}

constexpr bond_types btype(atom_types a, atom_types b)
{
    using ut = std::underlying_type_t<atom_types>;
    return static_cast<bond_types>(
        (static_cast<ut>(a) ^ static_cast<ut>(b)) ?
            (static_cast<ut>(a) ^ static_cast<ut>(b)) : static_cast<ut>(a)
    );
}

namespace detail {
constexpr auto btype_ut(atom_types a, atom_types b)
{
    using ut = std::underlying_type_t<bond_types>;
    return static_cast<ut>(btype(a, b));
}
} // namespace detail

enum class bond_types : std::uint32_t
{
    // common bond types listed for convenience
    NONE = detail::btype_ut(atype::NONE, atype::NONE),
    CC =   detail::btype_ut(atype::C, atype::C),
    CH =   detail::btype_ut(atype::C, atype::H),
    HH =   detail::btype_ut(atype::H, atype::H)
};

} // namespace sp2

SPECIALIZE_ENUM_HASH(sp2::atom_types)
SPECIALIZE_ENUM_HASH(sp2::bond_types)

#endif // SP2_ATOM_TYPES_HPP
