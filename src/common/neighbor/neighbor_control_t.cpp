#ifdef SP2_NEIGHBOR

#include "neighbor_control_t.hpp"
#include <common/math/mat3x3_t.hpp>
#include <common/util/random.hpp>

#include <gtest/gtest.h>
sp2::mat3x3_t gen_lattice(sp2::util::rng_t &rng, double r_min = 1, double r_max = 100)
{
    auto gen_vec = [&]{
        return sp2::random_vec3(rng.get_gen()) * rng.rand(r_min, r_max);
    };

    sp2::vec3_t lattice_vecs[3] = {
        gen_vec(), gen_vec(), gen_vec()
    };

    sp2::mat3x3_t lattice;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice[j][i] = lattice_vecs[i][j];

    return lattice;
}

TEST(lattice, all)
{
    auto rng = []{
        sp2::util::rng_t rng;
        rng.seed_random();
        return rng;
    }();

    constexpr int nt = 0,
        n_points = 1e2;

    std::vector<sp2::vec3_t> points(n_points);
    for (int t = 0; t < nt; ++t)
    {
        auto mat = gen_lattice(rng);
        for (auto &p : points)
            p = sp2::vec3_t(rng.rand(0.0, 1.0), rng.rand(0.0, 1.0), rng.rand(0.0, 1.0));

        for (auto &a : points)
        {
            for (auto &b : points)
            {
                auto delta = b - a;
                for (auto &coord : delta)
                    coord -= std::round(coord);

                double round_dist = (mat * delta).mag_sq();

                double brute_force_dist = std::numeric_limits<double>::max();
                for (int i = -1; i <= 1; ++i)
                    for (int j = -1; j <= 1; ++j)
                        for (int k = -1; k <= 1; ++k)
                            brute_force_dist = std::min(brute_force_dist,
                                (mat * (delta + sp2::vec3_t(i, j, k))).mag_sq());

                ASSERT_DOUBLE_EQ(round_dist, brute_force_dist);
            }
        }

        std::cout << "t: " << t << std::endl;
    }
}
#include <regex>

namespace sp2
{

/// seitz matrix type (4x4 transformation matrix)
using seitz_mat_t = scpp::array<double, 4, 4>;

enum class crystal_system : int
{
    invalid,
    triclinic,
    monoclinic,
    othorhombic,
    tetragonal,
    trigonal,
    hexagonal,
    cubic
};

enum class translational_symmetry : int
{
    invalid,
    /// P - Primitive
        primitive,
    /// A - Base centered on A faces only
        base_centered_a,
    /// B - Base centered on B faces only
        base_centered_b,
    /// C - Base centered on C faces only
        base_centered_c,
    /// I - Body centered (from the German "Innenzentriert")
        body_centered,
    /// F - Face centered (from the German "Flächenzentriert")
        face_centered,
    /// R - Rhombohedral
        rhombohedral
};

constexpr crystal_system get_crystal_sys(unsigned int space_group_number)
{
    // if (space_group_number == 0 || space_group_number > 230)
    //     return 0; // error

    // space groups are arranged such that all groups corresponding to
    // a particular crystal system are in contiguous blocks
    if (space_group_number < 3)
        return crystal_system::triclinic;

    if (space_group_number < 16)
        return crystal_system::monoclinic;

    if (space_group_number < 75)
        return crystal_system::othorhombic;

    if (space_group_number < 143)
        return crystal_system::tetragonal;

    if (space_group_number < 168)
        return crystal_system::trigonal;

    if (space_group_number < 195)
        return crystal_system::hexagonal;

    return crystal_system::cubic;
}

/// generic axis class used when generating symmetry transformations
class axis
{
private:
    int data[3];

    /// private constructor, the static member functions:
    ///     a(), b(), c(), and none()
    /// should be used to obtain transformation axes instead
    constexpr axis(int x, int y, int z) :
        data{x, y, z} {}
public:
    constexpr axis() : data{} {}

    constexpr axis(axis&&) = default;
    constexpr axis(const axis&) = default;

    constexpr axis& operator=(axis&&) = default;
    constexpr axis& operator=(const axis&) = default;

    /// read-only
    constexpr int operator[](std::size_t i) const {
        return data[i];}

    constexpr static axis a() {
        return axis(1, 0, 0);}

    constexpr static axis b() {
        return axis(0, 1, 0);}

    constexpr static axis c() {
        return axis(0, 0, 1);}

    constexpr static axis none() {
        return axis(0, 0, 0);}

    /// unary minus
    constexpr axis operator-() const {
        return axis(-data[0], -data[1], -data[2]);}

    constexpr axis& operator+=(const axis &rhs)
    {
        for (int i = 0; i < 3; ++i)
            data[i] += rhs.data[i];

        return *this;
    }

    constexpr axis& operator-=(const axis &rhs)
    {
        for (int i = 0; i < 3; ++i)
            data[i] -= rhs.data[i];

        return *this;
    }

    constexpr axis operator+(const axis &rhs) const {
        return axis(*this) += rhs;}

    constexpr axis operator-(const axis &rhs) const {
        return axis(*this) -= rhs;}

    /// unique conversion to integer so that we can use this in switch()
    /// note: an axis with any value > 1 is invalid (e.g a() + a())
    constexpr operator uint8_t() const
    {
        uint8_t val = 0;
        for (int i = 0; i < 3; ++i)
        {
            // Bitwise encoding:
            //   - The first bit represents whether the nth axis is nonzero.
            //   - The second bit represents the sign of the nth axis
            //     Note: it is always zero if the first bit is zero
            constexpr uint8_t is_positive = 1 | (1 << 1),
                is_negative = 1 | (0 << 1);

            if (data[i] < 0)
                val |= is_negative;
            else if (data[i] > 0)
                val |= is_positive;
            else
            {
                // data[i] == 0, do nothing
            }

            // shift val by two to store next axis
            if (i + 1 < 3)
                val <<= 2;
        }

        return val;
    }

    constexpr operator vec3_t() const
    {
        return {data[0], data[1], data[2]};
    }
};


/// Get default rotation axis based on what number rotation this is and what
/// the previous rotation's (if any) order was (set to zero or omit second
/// argument for n/a). returns axis::none() if no default is defined. The
/// default axis should be overwritten if it is explicitly specified later.
constexpr axis default_axis(int rotation_num,
    int rotation_order,
    int previous_order = 0,
    axis previous_axis = axis()
)
{
    // if rotation_order == 2 and ' or "
    // if prev_axis
    //  - x:
    //    - ' b - c
    //    - " b + c
    //  - y:
    //    - ' a - c
    //    - " a + c
    //  - z:
    //    - ' a - b
    //    - " a + b

    // if cubic and rotation_order == 3, a + b + c is implied

    // For most Hall symbols the rotation axes applicable to each N are implied
    // and an explicit axis symbol A is not needed. The rules for default axis
    // directions are:
    //
    switch (rotation_num)
    {
        // 1. The first rotation has an axis direction of:
        //        c
    case 1:
        return axis::c();

        // 2. The second rotation (if N is 2) has an axis direction of:
        //        a      if preceded by an N of 2 or 4
        //        a - b  if preceded by an N of 3 or 6
    case 2:
        if (previous_order == 2 || previous_order == 4)
            return axis::a();
        else if (previous_order == 3 || previous_order == 6)
            return axis::a() - axis::b();
        // no default
        break;

        // 3. The third rotation (N is always 3) has an axis direction of:
        //        a + b + c
    case 3:
        if (previous_order == 3)
            return axis::a() + axis::b() + axis::c();
        // no default
        break;

    default:
        // no default
        break;
    }

    // the axis must be specified if none of the conditions match,
    // so we return a value of axis::none() to indicate this
    return axis::none();
}

constexpr seitz_mat_t origin_shift(seitz_mat_t matrix, vec3_t translation)
{
    // The origin-shift translation vector V has the construction (va vb vc),
    // where va, vb and vc denote the shifts in 12ths parallel to the cell
    // edges a, b and c, respectively. va/12, vb/12 and vc/12 are the
    // coordinates of the unshifted origin in the shifted basis system. The
    // shifted Seitz matrices Sn' are derived from the unshifted matrices Sn
    // with the transformation
    //
    //       (1 0 0 va/12)        (1 0 0 -va/12)
    // Sn' = (0 1 0 vb/12) * Sn * (0 1 0 -vb/12)
    //       (0 0 1 vc/12)        (0 0 1 -vc/12)
    //       (0 0 0   1  )        (0 0 0    1  )

    seitz_mat_t temp = {};
    for (int i = 0; i < 4; ++i)
    {
        temp[i][i] = 1;
        if (i < 3)
            temp[3][i] = translation[i] / 12;
    }

    // first multiplication (left)
    matrix = temp * matrix;

    // reverse sign on offsets
    for (int i = 0; i < 3; ++i)
        temp[3][i] *= -1;

    // second multiplication (right)
    return matrix * temp;
}

/// returns translations
std::vector<vec3_t> get_lattice_translations(char lattice_symbol)
{
    // The lattice symbol L specifies one or more Seitz matrices which are
    // needed to generate the space-group symmetry elements. For
    // noncentrosymmetric lattices the rotation matrices are for 1.
    // For centrosymmetric lattices the lattice symbols are preceded by a minus
    // sign `-', rotations are 1 and -1, and the total number of generator
    // matrices implied by each symbol is twice the number of implied lattice
    // translations.
    switch (lattice_symbol)
    {
    case 'p':
        // P - Primitive
        return {
            vec3_t(0, 0, 0)
        };
    case 'a':
        // A - Base centered on A faces only
        return {
            vec3_t(0, 0, 0),
            vec3_t(0, 1, 1) / 2
        };
    case 'b':
        // B - Base centered on B faces only
        return {
            vec3_t(0, 0, 0),
            vec3_t(1, 0, 1) / 2
        };
    case 'c':
        // C - Base centered on C faces only
        return {
            vec3_t(0, 0, 0),
            vec3_t(1, 1, 0) / 2
        };
    case 'i':
        // I - Body centered (from the German "Innenzentriert")
        return {
            vec3_t(0, 0, 0),
            vec3_t(1, 1, 1) / 2
        };
    case 'f':
        // F - Face centered (from the German "Flächenzentriert")
        return {
            vec3_t(0, 0, 0),
            vec3_t(0, 1, 1) / 2,
            vec3_t(1, 0, 1) / 2,
            vec3_t(1, 1, 0) / 2
        };
    case 'r':
        // R - Rhombohedral
        return {
            vec3_t(0, 0, 0),
            vec3_t(2, 1, 1) / 3,
            vec3_t(1, 2, 2) / 3
        };
        // TODO S, T?
    case 's':
    case 't':
    default:
        throw std::invalid_argument("Unknown lattice symbol '" +
                                    std::string(1, lattice_symbol) + "' in Hall Notation.");
    }
}

constexpr vec3_t get_alpha_translation(char translation_symbol)
{
    // The symbol T specifies the translation elements of a Seitz matrix.
    // Alphabetical symbols specify translations along a fixed direction.
    translation_symbol = static_cast<char>(std::tolower(translation_symbol));

    switch (translation_symbol)
    {
    case 'a': return vec3_t(1, 0, 0) / 2;
    case 'b': return vec3_t(0, 1, 0) / 2;
    case 'c': return vec3_t(0, 0, 1) / 2;
    case 'n': return vec3_t(1, 1, 1) / 2;
    case 'u': return vec3_t(1, 0, 0) / 4;
    case 'v': return vec3_t(0, 1, 0) / 4;
    case 'w': return vec3_t(0, 0, 1) / 4;
    case 'd': return vec3_t(1, 1, 1) / 4;
    default:
        throw std::invalid_argument(
            "Unknown alphabetical translation symbol '"
            + std::string(1, translation_symbol)
            + "' in Hall Notation.");
    }
}

constexpr vec3_t get_subscript_translation(axis rotation_axis,
    int rotation_order, int subscript_value)
{

    return static_cast<vec3_t>(rotation_axis) *
           (rotation_order / static_cast<double>(subscript_value));
}

constexpr seitz_mat_t get_rotation(axis rotation_axis,
    int rotation_order)
{
    // make sure the rotation order is valid (only catches errors with
    // rotations for the principal axes)
    if (rotation_order < 1 || rotation_order > 6)
        return seitz_mat_t{};

    // "1-Fold" aka identity, doesn't matter what axis
    if (rotation_order == 1)
    {
        return {
            { 1,  0,  0,  0},
            { 0,  1,  0,  0},
            { 0,  0,  1,  0},
            { 0,  0,  0,  1}
        };
    }

    // we can use a switch statement since conversion of rotation_axis to
    // uint8_t is unique for all valid axes (and is constexpr)
    switch (static_cast<uint8_t>(rotation_axis))
    {
////////////////////////////////////////////////////////////////////////////////
// Rotation matrices for the principal axes:                                  //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Rotation Axis: a                                                           //
////////////////////////////////////////////////////////////////////////////////
    case axis::a():
        // five rotation matrices, which one we return depends on the
        // order of the rotation (the N in N-fold rotation)
        switch (rotation_order)
        {
        case 2:
            // 2-Fold
            return {
                { 1,  0,  0,  0},
                { 0, -1,  0,  0},
                { 0,  0, -1,  0},
                { 0,  0,  0,  1}
            };
        case 3:
            // 3-Fold
            return {
                { 1,  0,  0,  0},
                { 0,  0, -1,  0},
                { 0,  1, -1,  0},
                { 0,  0,  0,  1}
            };
        case 4:
            // 4-Fold
            return {
                { 1,  0,  0,  0},
                { 0,  0, -1,  0},
                { 0,  1,  0,  0},
                { 0,  0,  0,  1}
            };
        case 6:
            // 6-Fold
            return {
                { 1,  0,  0,  0},
                { 0,  1, -1,  0},
                { 0,  1,  0,  0},
                { 0,  0,  0,  1}
            };
        case 5:
            // Note: no 5-fold rotation
        default:
            return {};
        }

////////////////////////////////////////////////////////////////////////////////
// Rotation Axis: b                                                           //
////////////////////////////////////////////////////////////////////////////////
    case axis::b():
        // five rotation matrices, which one we return depends on the
        // order of the rotation (the N in N-fold rotation)
        switch (rotation_order)
        {
        case 2:
            // 2-Fold
            return {
                {-1,  0,  0,  0},
                { 0,  1,  0,  0},
                { 0,  0, -1,  0},
                { 0,  0,  0,  1}
            };
        case 3:
            // 3-Fold
            return {
                {-1,  0,  1,  0},
                { 0,  1,  0,  0},
                {-1,  0,  0,  0},
                { 0,  0,  0,  1}
            };
        case 4:
            // 4-Fold
            return {
                { 0,  0,  1,  0},
                { 0,  1,  0,  0},
                {-1,  0,  0,  0},
                { 0,  0,  0,  1}
            };
        case 6:
            // 6-Fold
            return {
                { 0,  0,  1,  0},
                { 0,  1,  0,  0},
                {-1,  0,  1,  0},
                { 0,  0,  0,  1}
            };
        case 5:
            // Note: no 5-fold rotation
        default:
            return {};
        }

////////////////////////////////////////////////////////////////////////////////
// Rotation Axis: c                                                           //
////////////////////////////////////////////////////////////////////////////////
    case axis::c():
        // five rotation matrices, which one we return depends on the
        // order of the rotation (the N in N-fold rotation)
        switch (rotation_order)
        {
        case 2:
            // 2-Fold
            return {
                {-1,  0,  0,  0},
                { 0, -1,  0,  0},
                { 0,  0,  1,  0},
                { 0,  0,  0,  1}
            };
        case 3:
            // 3-Fold
            return {
                { 0, -1,  0,  0},
                { 1, -1,  0,  0},
                { 0,  0,  1,  0},
                { 0,  0,  0,  1}
            };
        case 4:
            // 4-Fold
            return {
                { 0, -1,  0,  0},
                { 1,  0,  0,  0},
                { 0,  0,  1,  0},
                { 0,  0,  0,  1}
            };
        case 6:
            // 6-Fold
            return {
                { 1, -1,  0,  0},
                { 1,  0,  0,  0},
                { 0,  0,  1,  0},
                { 0,  0,  0,  1}
            };
        case 5:
            // Note: no 5-fold rotation
        default:
            return {};
        }

////////////////////////////////////////////////////////////////////////////////
// Rotation matrices for face-diagonal axes:                                  //
////////////////////////////////////////////////////////////////////////////////

        // Axis: b - c
    case axis::b() - axis::c():
        return {
            {-1,  0,  0,  0},
            { 0,  0, -1,  0},
            { 0, -1,  0,  0},
            { 0,  0,  0,  1}
        };

        // Axis: b + c
    case axis::b() + axis::c():
        return {
            {-1,  0,  0,  0},
            { 0,  0,  1,  0},
            { 0,  1,  0,  0},
            { 0,  0,  0,  1}
        };

        // Axis: a - c
    case axis::a() - axis::c():
        return {
            { 0,  0, -1,  0},
            { 0, -1,  0,  0},
            {-1,  0,  0,  0},
            { 0,  0,  0,  1}
        };

        // Axis: a + c
    case axis::a() + axis::c():
        return {
            { 0,  0,  1,  0},
            { 0, -1,  0,  0},
            { 1,  0,  0,  0},
            { 0,  0,  0,  1}
        };

        // Axis: a - b
    case axis::a() - axis::b():
        return {
            { 0, -1,  0,  0},
            {-1,  0,  0,  0},
            { 0,  0, -1,  0},
            { 0,  0,  0,  1}
        };

        // Axis: a + b
    case axis::a() + axis::b():
        return {
            { 0,  1,  0,  0},
            { 1,  0,  0,  0},
            { 0,  0, -1,  0},
            { 0,  0,  0,  1}
        };

////////////////////////////////////////////////////////////////////////////////
// Rotation matrix for the body-diagonal axis:                                //
////////////////////////////////////////////////////////////////////////////////

        // Axis: a + b + c
    case axis::a() + axis::b() + axis::c():
        return {
            { 0,  0,  1,  0},
            { 1,  0,  0,  0},
            { 0,  1,  0,  0},
            { 0,  0,  0,  1}
        };

    default:
        // invalid axis
        return {};
    }
}

constexpr seitz_mat_t get_rotation(axis rotation_axis,
    int rotation_order, bool improper)
{
    auto rotation = get_rotation(rotation_axis, rotation_order);
    if (!improper)
        return rotation;

    return -1 * rotation;
}


std::vector<seitz_mat_t> parse_hall_notation(std::string input)
{
    // For format information see:
    //  - http://cci.lbl.gov/sginfo/hall_symbols.html
    //  - [1] S.R. Hall; Space-Group Notation with an Explicit Origin;
    //        Acta Cryst. (1981). A37, 517-525
    //  - [2] International Tables Volume B 1994, Section 1.4. Symmetry in
    //        reciprocal space

    // Hall Symbols
    //
    // The explicit-origin space group notation proposed by Hall (1981) [1][2]
    // is based on the minimum number of symmetry operations, in the form of
    // Seitz matrices, needed to uniquely define a space group. The concise
    // unambiguous nature of this notation makes it well suited to handling
    // symmetry in computing and database applications.
    //
    // The notation has the general form:
    //     L [NAT]1 ... [NAT]p V
    // where:
    //  - L is the symbol specifying the lattice translational symmetry.
    //  - NAT identifies the 4x4 Seitz matrix of a symmetry element in the
    //      minimum set which defines the space-group symmetry.
    //  - p is the number of elements in the set,
    //  - V is a translation vector which shifts the origin of the generator
    //      matrices by fractions of the unit cell lengths a, b and c.

    // convert everything to lowercase to make things easier
    for (auto &c : input)
        c = static_cast<char>(std::tolower(c));

    // the result vector of seitz matrices
    std::vector<seitz_mat_t> matrices;
    // initial is just the identity matrix
    matrices.emplace_back({
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    });

    // The lattice symbol L is the first element (ignoring whitespace).
    // For centrosymmetric lattices the lattice symbols are preceded by a minus
    // sign `-', rotations are 1 and -1, and the total number of generator
    // matrices implied by each symbol is twice the number of implied lattice
    // translations.
    static const std::regex lattice_symbol_regex("^\\s*(-?[a-z])\\s+");

    std::smatch match;
    std::regex_search(input, match, lattice_symbol_regex);

    bool lat_error = (match.size() != 2);

    if (!lat_error)
    {
        // first sub-match contains L (ignores whitespace)
        std::string lattice_symbol = match[1];

        // remove the portion of 'input' which matched lattice_symbol_regex
        // to simplify later parsing
        input = match.suffix();

        // if centrosymmetric, lattice symbol is prefixed with '-'
        if (lattice_symbol[0] == '-')
        {
            matrices.emplace_back({
                {-1,  0,  0, 0},
                { 0, -1,  0, 0},
                { 0,  0, -1, 0},
                { 0,  0,  0, 1}
            });
        }

        // get the translations associated with the lattice symbol
        std::vector<vec3_t> lattice_trans = get_lattice_translations(
            lattice_symbol.back());

        if (lattice_trans.empty())
            lat_error = true;
        else
        {

        }
    }

    if (lat_error)
    {
        std::cerr << "Failed to parse lattice symbol in Hall notation input. "
                  << "Input string was '" << input << "'." << std::endl;

        return {};
    }

    // The origin shift translation vector V is the last element and we parse
    // it first before parsing the 4x4 Seitz matrix (NAT) elements to remove
    // some parsing ambiguities. Note: V looks like (0 0 -1)
    static const std::regex origin_shift_regex("\\(([^\\)]+)\\)\\s+$");

    std::regex_search(input, match, origin_shift_regex);

    std::string origin_shift;
    if (!match.empty())
    {
        // first sub-match contains the vector components
        origin_shift = match[1];

        // remove from input to simplify parsing
        input = match.prefix();
    }

    // The matrix symbol NAT is composed of three parts:
    static const std::regex seitz_regex(
        "(\\b-?[1-6])"        // N is the symbol denoting the n-fold order of
            //   the rotation matrix
            "([\"*'xyz]?)"        // A is a superscript symbol denoting
            //   the axis of rotation
            "([abcnuvwd1-6]*\\b)" // T is a subscript symbol denoting the
        //   translation vector (multiple are additive)
    );

    // find all matches
    auto seitz_matches = make_range(
        std::sregex_iterator(input.begin(), input.end(), seitz_regex),
        std::sregex_iterator()
    );

    // iterate through them
    int rotation_num   = 1,
        last_rot_order = 0;
    axis last_rot_axis = axis::none();

    for (auto &match : seitz_matches)
    {

    }
    return {};
}

struct symmetry_info_t
{
    int space_group_num;
    crystal_system sys;

    /// transformations used to generate all symmetry transformations
    std::vector<seitz_mat_t> base_transforms,
    /// all symmetry transformations
        full_transforms;


};

constexpr char simple_tolower(char c)
{
    return c < 'A' ? c :
           c > 'Z' ? c : (c - 'A') + 'a';
}

class hall_notation_parser_t
{
private:
    crystal_system sys = crystal_system::invalid;

    bool centrosymmetric = false;
    translational_symmetry trans_symm = translational_symmetry::invalid;

    int rotation_number = 1,
        rotation_order = 0;
    axis rotation_axis = axis::none();

    bool has_error = true;
    const char *current = nullptr;
public:
    hall_notation_parser_t() = default;

    hall_notation_parser_t(const char* str) :
        hall_notation_parser_t()
    {
        current = str;
        if (current != nullptr && *current != '\0')
            has_error = false;
    }

    operator bool() const {
        return !has_error;
    }

    bool parse()
    {
        if (has_error)
            return false;

        // first part is the lattice symbol
        if (!skip_whitespace())
            return false; // hit end of string before finding non-whitespace

        trans_symm = parse_lattice_symbol();
        if (trans_symm == translational_symmetry::invalid)
            return false;

        while (!has_error)
        {
            if (!skip_whitespace())
                return false;

            // char is at beginning of NAT
        }

        return !has_error;
    }

private:
    bool next()
    {
        return *++current;
    }

    char current_char() const
    {
        return simple_tolower(*current);
    }

    seitz_mat_t parse_matrix_symbol()
    {
        // The matrix symbol NAT is composed of three parts:
        // N is the symbol denoting the n-fold order of
        //   the rotation matrix
        // A is a superscript symbol denoting
        //   the axis of rotation
        // T is a subscript symbol denoting the
        //   translation vector (multiple are additive)

        // First is N:
        //   The 3x3 matrices for proper rotations along the three principal
        //   unit-cell directions. The matrices for improper rotations (-1, -2,
        //   -3, -4 and -6) are identical except that the signs are reversed.
        bool improper_rotation = false;
        if (*current == '-')
        {
            improper_rotation = true;
            current++;
        }

        int new_rotation_order = (*current - '0');
    }

    translational_symmetry parse_lattice_symbol()
    {
        // The lattice symbol L specifies one or more Seitz matrices which are
        // needed to generate the space-group symmetry elements.
        //
        // For noncentrosymmetric lattices the rotation matrices are for 1.
        // For centrosymmetric lattices the lattice symbols are preceded by a
        // minus sign `-', rotations are 1 and -1, and the total number of
        // generator matrices implied by each symbol is twice the number of
        // implied lattice translations.
        //
        if (*current++ == '-')
            centrosymmetric = true;

        if (!*current)
            return translational_symmetry::invalid;

        switch (*current++)
        {
            // P - Primitive
        case 'p': return translational_symmetry::primitive;
            // A - Base centered on A faces only
        case 'a': return translational_symmetry::base_centered_a;
            // B - Base centered on B faces only
        case 'b': return translational_symmetry::base_centered_b;
            // C - Base centered on C faces only
        case 'c': return translational_symmetry::base_centered_c;
            // I - Body centered (from the German "Innenzentriert")
        case 'i': return translational_symmetry::body_centered;
            // F - Face centered (from the German "Flächenzentriert")
        case 'f': return translational_symmetry::face_centered;
            // R - Rhombohedral
        case 'r': return translational_symmetry::rhombohedral;

        default:
            // exclude the 'unusual lattice symbols' S and T
        case 't':
        case 's':
            has_error = true;
            return translational_symmetry::invalid;
        }
    }

    bool skip_whitespace()
    {
        while (*current && (*current == ' ' || *current == '\t'))
            ++current;

        return *current;
    }
};



constexpr symmetry_info_t parse(const char* str, std::size_t len)
{
    return {};
}

template<std::size_t N>
constexpr symmetry_info_t parse(const char (&str)[N])
{
    return parse(str, N - 1);
}

}

#endif // SP2_ENABLE_TESTS
