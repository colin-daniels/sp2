#include "common/neighbor/utility_functions.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <stdexcept>

using namespace std;
using namespace sp2;

void fbc::array_rot(const int axis, double *input,
    double *output, const double theta)
{
    double ct = std::cos(theta),
        st = std::sin(theta);

    double R[3][3] = {{ct, 0, 0},
                      {0, ct, 0},
                      {0, 0, ct}};

    if (axis > 2 || axis < 0)
        cout << "Error, invalid axis (" << axis << ") passed to array_rot."
             << endl;

    R[(axis + 1) % 3][(axis + 2) % 3] = -st;
    R[(axis + 2) % 3][(axis + 1) % 3] = st;
    R[axis][axis] = 1;

    // apply the transformation
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            output[i * 3 + j] = 0;
            for (int k = 0; k < 3; ++k)
                output[i * 3 + j] += R[i][k] * input[k * 3 + j];
        }
    }
}

void rotate_lattice(double lattice[3][3], double R[3][3])
{
    double I0[3][3] = {};
    for (int i = 0; i < 3; ++i)
        I0[i][i] = 1;

    double temp[3][3] = {};
    copy_n(lattice[0], 9, temp[0]);

    // rotate to get into upper triangular form
    fbc::array_rot(2, I0[0], R[0], -atan2(temp[1][0], temp[0][0]));
    fbc::array_rot(2, temp[0], lattice[0], -atan2(temp[1][0], temp[0][0]));

    fbc::array_rot(1, R[0], I0[0], atan2(lattice[2][0], lattice[0][0]));
    fbc::array_rot(1, lattice[0], temp[0], atan2(lattice[2][0], lattice[0][0]));

    fbc::array_rot(0, I0[0], R[0], -atan2(temp[2][1], temp[1][1]));
    fbc::array_rot(0, temp[0], lattice[0], -atan2(temp[2][1], temp[1][1]));

    // prefer positive diagonals
    for (int i = 0; i < 3; ++i)
    {
        if (lattice[i][i] < 0)
        {
            for (int j = 0; j < 3; ++j)
            {
                lattice[j][i] = -lattice[j][i];
                R[j][i] = -R[j][i];
            }
        }
    }
}

void fbc::invert_3x3(const double input[3][3], double inverse[3][3])
{
    // 1 / det
    // 33 32 13 12 23 22
    // 23 22 33 32 13 12
    // 31 33 11 13 21 23
    // 21 23 31 33 11 13
    // 32 31 12 11 22 21
    // 22 21 32 31 12 11
    for (int i = 0; i < 3; ++i)
    {
        int jx1 = ((i + 2) % 3),
            jx2 = ((i + 1) % 3);

        for (int j = 0; j < 3; ++j)
        {
            int ix1 = ((j + 2) % 3),
                ix2 = ((j + 1) % 3);

            inverse[i][j] = input[ix1][jx1] * input[ix2][jx2]
                            - input[ix1][jx2] * input[ix2][jx1];
        }
    }

    double det = input[0][0] * inverse[0][0]
                 + input[0][1] * inverse[1][0]
                 + input[0][2] * inverse[2][0];

    if (det == 0)
    {
        cout << "Error, matrix passed to invert_3x3 is singular."
             << endl;
    }
    else
    {
        for (int i = 0; i < 9; ++i)
            inverse[i/3][i%3] /= det;
    }
}

double det_3x3(double input[3][3])
{
    double result = 0;
    for (int i = 0; i < 3; ++i)
        result += input[0][i] * (input[1][(i + 1) % 3] * input[2][(i + 2) % 3]
                                 - input[1][(i + 2) % 3] * input[2][(i + 1) % 3]);
    return result;
}

void fill_lattice(double lattice[3][3], int n_periodic)
{
    int n_total = n_periodic;

    // depending on the number of vectors, generate one, two, or
    // three temporary (orthogonal) lattice vectors
    if (n_total == 0)
    {
        lattice[0][0] = 1;
        lattice[1][1] = 1;
        lattice[2][2] = 1;
        n_total = 3;
    }

    if (n_total == 1)
    {
        int a = 0,
            id = 1;

        vector<int> nnz;
        for (int i = 0; i < 3; ++i)
            if (lattice[i][a] != 0)
                nnz.push_back(i);

        if (nnz.size() == 1)
            lattice[(nnz[0] + 1) % 3][id] = 1.0;
        else if (nnz.size() == 2)
        {
            lattice[nnz[0]][id] = -lattice[nnz[1]][a];
            lattice[nnz[1]][id] =  lattice[nnz[0]][a];
        }
        else
        {
            lattice[0][id] = -lattice[1][a];
            lattice[1][id] =  lattice[0][a];
        }

        n_total = 2;
    }

    if (n_total == 2)
    {
        int a  = 0,
            b  = 1,
            id = 2;

        // cross product
        lattice[0][id] = lattice[1][a] * lattice[2][b]
                         - lattice[2][a] * lattice[1][b];
        lattice[1][id] = lattice[2][a] * lattice[0][b]
                         - lattice[0][a] * lattice[2][b];
        lattice[2][id] = lattice[0][a] * lattice[1][b]
                         - lattice[1][a] * lattice[0][b];

        double len = sqrt(lattice[0][id] * lattice[0][id] +
                          lattice[1][id] * lattice[1][id] +
                          lattice[2][id] * lattice[2][id]);

        for (int i = 0; i < 3; ++i)
            lattice[i][id] /= len;
    }
}

int fbc::process_lattice(const double lattice[3][3],
    double mod_lattice[3][3], double transformation[3][3])
{
    // copy the original lattice into the modified one
    copy_n(lattice[0], 9, mod_lattice[0]);

    // get the lattice vector lengths
    double lattice_len[3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice_len[i] += mod_lattice[j][i] * mod_lattice[j][i];

    int n_periodic = 0;
    for (int i = 0; i < 3; ++i)
    {
        lattice_len[i] = sqrt(lattice_len[i]);
        if (lattice_len[i] > 0)
            ++n_periodic;
    }

    // put all the 'existing' (nonzero length) lattice vectors in the
    // beginning of the lattice matrix (the first columns)
    for (int i = 0; i < 2; ++i)
    {
        if (lattice_len[i] > 0)
            continue;
        // find a nonzero vector
        for (int j = i + 1; j < 3; ++j)
        {
            // swap it
            if (lattice_len[j] > 0)
            {
                swap(lattice_len[i], lattice_len[j]);
                for (int k = 0; k < 3; ++k)
                    swap(mod_lattice[k][i], mod_lattice[k][j]);
                break;
            }
        }
    }

    // generate orthogonal vectors if less than
    // three lattice vectors were defined
    fill_lattice(mod_lattice, n_periodic);

    // check linear independence
    if (det_3x3(mod_lattice) == 0)
        throw invalid_argument(
            "lattice vector input has linearly dependent vectors.");

    // rotate the lattice so that the first and second
    // lattice vectors are parallel to the xy-plane
    rotate_lattice(mod_lattice, transformation);

    // return number of periodic directions
    return n_periodic;
}
