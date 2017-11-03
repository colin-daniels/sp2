/************************************************************
* Nicholas Hackney & Colin Daniels           				*
* Rensellaer Polytechnic Institute          				*
* Fall 2017                                  				*
* Version 1.0                                				*
*                                            				*
* Implementation of the classical phosphorene  				*
* potential from Phys. Chem. Chem. Phys., 2016, 18, 23312   *
* Midtvedt & Croy                                           *
************************************************************/

#include "ephos.hpp"
#include <fstream>
#include "src/common/math/vec3_t.hpp"
#include "src/common/math/vec3_util.hpp"
#include "src/common/io/structure.hpp"
#include <src/common/math/numerical_diff.hpp>
#include <src/common/math/blas.hpp>
#include <src/common/util/random.hpp>
#include <vector>

using namespace sp2;

namespace {

inline std::vector<vec3_t> get_pristine_bond_array()
{
    constexpr double theta1 =  96.5 * (M_PI / 180), // same pucker horizontal
                     theta2 = 101.9 * (M_PI / 180), // different pucker vertical
         equilibrium_length =  2.22; // angstroms

    static const std::vector<vec3_t> pristine_bonds = []{
        vec3_t bonds[32] = {
            vec3_t{-std::cos(theta1 / 2.0), -std::sin(theta1 / 2.0), 0},
            vec3_t{-std::cos(theta1 / 2.0),  std::sin(theta1 / 2.0), 0},
            vec3_t{std::cos(theta1 / 2.0), -std::sin(theta1 / 2.0), 0},
            vec3_t{std::cos(theta1 / 2.0),  std::sin(theta1 / 2.0), 0},

            vec3_t{ std::cos(theta2)/std::cos(theta1/2.0), 0, std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},
         	vec3_t{ std::cos(theta2)/std::cos(theta1/2.0), 0, -std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},
         	vec3_t{ -std::cos(theta2)/std::cos(theta1/2.0), 0, std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},
         	vec3_t{ -std::cos(theta2)/std::cos(theta1/2.0), 0, -std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},

         	vec3_t{0, std::cos(theta2)/std::cos(theta1/2.0), std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},
         	vec3_t{0, std::cos(theta2)/std::cos(theta1/2.0), -std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},
         	vec3_t{0, -std::cos(theta2)/std::cos(theta1/2.0), std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},
         	vec3_t{0, -std::cos(theta2)/std::cos(theta1/2.0), -std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0)))},

         	vec3_t{-std::sin(theta1 / 2.0), -std::cos(theta1 / 2.0), 0},
            vec3_t{std::sin(theta1 / 2.0), -std::cos(theta1 / 2.0), 0},
            vec3_t{-std::sin(theta1 / 2.0), std::cos(theta1 / 2.0), 0},
            vec3_t{std::sin(theta1 / 2.0), std::cos(theta1 / 2.0), 0}
        };

        // reverse directions
        for (int i = 16; i < 32; ++i)
            bonds[i] = -bonds[i - 16];

        // scale to correct length
        for (auto &b : bonds)
            b *= equilibrium_length;

        std::vector<vec3_t> temp;
        for (int i = 0; i < 32; ++i)
            temp.push_back(bonds[i]);

        return temp;
    }();

    return pristine_bonds;
}

inline vec3_t get_pristine_vec(vec3_t bond)
{
    static const auto pristine_bonds = get_pristine_bond_array();

    int min_id = -1;
    double min_diff = std::numeric_limits<double>::max();

    for (int i = 0; i < pristine_bonds.size(); ++i)
    {
        double diff = (pristine_bonds[i] - bond).mag_sq();
        if (diff < min_diff)
        {
            min_id = i;
            min_diff = diff;
        }
    }
    assert (min_id!=-1);
    return pristine_bonds[min_id];
}

// std::vector<double> get_pristine()
// {
// 	structure_t output;
// 	io::read_structure("pristine.vasp", output);
// 	double bond_cutoff=3.0;
// 	fbc::bond_control_t bond_control(output.lattice, bond_cutoff, 0);
// 	bond_control.update(v3tod(output.positions));
// 	bond_control.lock_bonds();
// 	return bond_control.get_bond_deltas();
// }

void test_gradient(phos::phosphorene_sys_t &sys)
{
    constexpr double atom_offset_max = 0.15; // Angstroms
    static auto rng = []{
        util::rng_t temp;
        temp.seed_random();

        return temp;
    }();

    sys.update();
    const auto initial_pos = sys.get_position();

    auto get_rand_pos = [&rng, initial_pos]{
        std::vector<vec3_t> random_pos = dtov3(initial_pos);
        for (auto &p : random_pos)
            p += random_vec3(rng.get_gen()) * atom_offset_max;

        return v3tod(random_pos);
    };

    auto test_rand_dir = [&sys, &rng](std::vector<double> pos) {
        sys.set_position(pos);
        sys.update();

        std::vector<double> direction(pos.size(), 0);
        for (auto &v : direction)
            v = rng.rand(0.0, 1.0);

        // normalize random direction
        vnormalize(direction);

        // get analytical gradient in the random direction
        double analytical_grad = vdot(sys.get_gradient(), direction);

        auto objective_fn = sys.get_diff_fn();
        auto one_dim_func = [&](double a) {
            auto temp = pos;
            sp2::vaxpy(a, direction, temp);
            return objective_fn(temp).first;
        };

        // calculate numerical gradient
        double numerical_grad = util::central_difference<9, double>(
            one_dim_func, 0.0, 1e-5);

        return std::make_pair(
            numerical_grad, analytical_grad
        );
    };

    const int n_dir_trials = 100,
        n_offset_trials = 100;

    double max_rel_diff = std::numeric_limits<double>::lowest();
    double numerical_at_max = 0,
        analytical_at_max = 0;

    for (int i = 0; i < n_offset_trials; ++i)
    {
        // get random atom positions to test gradient at
        auto pos = get_rand_pos();
        for (int j = 0; j < n_dir_trials; ++j)
        {
            // test gradient in random direction at the given position
            auto grad = test_rand_dir(pos);
            auto relative = std::abs((grad.first - grad.second) / grad.first);

            if (relative > max_rel_diff)
            {
                max_rel_diff = relative;
                numerical_at_max = grad.first;
                analytical_at_max = grad.second;
            }

        }
    }

    std::cout << "Max relative difference: " << max_rel_diff << '\n'
              << "  Analytical gradient: " << analytical_at_max << '\n'
              << "  Numerical (actual): " << numerical_at_max << std::endl;

    exit(0);
}

vec3_t sin_grad(vec3_t r1, vec3_t r2)
{
	vec3_t grad;
	vec3_t cp= cross(r1,r2);
	double mcp= cp.mag(),
	       magr1=r1.mag(),
	       magr2=r2.mag(),
	       msr1=r1.mag_sq(),
	       msr2=r2.mag_sq();

	vec3_t num;
	num[0]=r1[2]-r2[2]+r2[1]-r1[1];
	num[1]=r2[2]-r1[2]+r1[0]-r2[0];
	num[2]=r1[1]-r2[1]+r2[0]-r1[0];

	grad=((num/mcp)*magr1*magr2-mcp*(r1*(magr2/magr1)+r2*(magr1/magr2)))/(msr1*msr2);

	return grad;
}

vec3_t cos_grad(vec3_t r1, vec3_t r2)
{
	vec3_t num1=(r1+r2)*r1.mag()*r2.mag();
	vec3_t num2 =dot(r1,r2)*(r1*(r2.mag()/r1.mag())+r2*(r1.mag()/r2.mag()));

	vec3_t grad=(num1-num2)/(r1.mag_sq()*r2.mag_sq());

	return grad;
}

} // anonymous namespace

phos::phosphorene_sys_t::phosphorene_sys_t(const sp2::structure_t &input) :
    structure(input),
    potential(0.0),
    forces(structure.types.size(), {0.0, 0.0, 0.0}),
    bond_control(structure.lattice, bond_cutoff, 0)
{
    // do initial update to calculate bonds
    bond_control.update(v3tod(structure.positions));
    // and immediately lock them since this potential is not a dynamic one
    bond_control.lock_bonds();
    bond_control.update(v3tod(structure.positions));
    for (vec3_t b : dtov3(bond_control.get_bond_deltas()))
    {
    	pristine_deltas.push_back(get_pristine_vec(b));
    }
    test_gradient(*this); //////////////////////////////////////////////////////************************//
    update();
}

diff_fn_t phos::phosphorene_sys_t::get_diff_fn()
{
    return [this](const auto &pos) {
        this->set_position(pos);
        this->update();

        return std::make_pair(this->get_value(), this->get_gradient());
    };
}

void phos::phosphorene_sys_t::set_position(const std::vector<double> &input)
{
    structure.positions = dtov3(input);
}

std::vector<double> phos::phosphorene_sys_t::get_gradient() const
{
    return v3tod(forces);
}

std::vector<double> phos::phosphorene_sys_t::get_position() const
{
    return v3tod(structure.positions);
}

double phos::phosphorene_sys_t::get_value() const
{
    return potential;
}

void phos::phosphorene_sys_t::update()
{
    // reset state
    potential = 0.0;  
    forces = std::vector<vec3_t>(structure.types.size(), {0, 0, 0});
    // update bonds from positions
    bond_control.update(v3tod(structure.positions));
	auto graph = bond_control.get_graph();
    auto bond_deltas = dtov3(bond_control.get_bond_deltas());

    // Define the constant parameters used in the calculations
	constexpr double Kr=11.17,
		KrP=10.3064,  // P denotes a ' (prime) on the variable
		Kt=1.18, // t is short for theta
		KtP=0.9259,
		Krrp=-0.6763, // p denotes a ' (prime) on the coordinate r
		KrrpP=1.2449,
		Krt=0.58,
		KrtP=1.932,
		KrtPP=0.797,
		d=2.22, // interatomic spacing [units: Angstroms]
		ds=4.9284;

	// double theta1 =  96.5 * (M_PI / 180);
	// double theta2 = 101.9 * (M_PI / 180);
	// vec3_t test= vec3_t(std::cos(theta2)/std::cos(theta1/2.0), 0, std::sqrt(1-(std::cos(theta2)/cos(theta1/2.0))*(std::cos(theta2)/cos(theta1/2.0))));
 //    std::cout<<test.x()<<" "<<test.y()<<" "<<test.z()<<std::endl;

    auto output_atom_into = [&](int id_a) {
        auto output_vec = [&](std::string name, vec3_t v) {
            std::cout << "      " << name << ": [ ";
            for (auto x : v)
                std::cout << x << ' ';
            std::cout << "]\n";
        };

        std::cout << "Atom ID: " << id_a << '\n'
                  << "  Type: " << sp2::enum_to_str(structure.types[id_a]) << '\n'
                  << "  Bonds: " << graph.degree(id_a) << std::endl;

        for (sp2::graph::ud_edge_t edge : graph.edges(id_a))
        {
            std::cout << "    ID " << edge.id << " "
                      << "[" << edge.a << "->" << edge.b << "]\n";

            output_vec("Components", bond_deltas[edge.id]);
            output_vec("Pristine", pristine_deltas[edge.id]);

            std::cout << "      Delta Mag: " << (pristine_deltas[edge.id] - bond_deltas[edge.id]).mag()
                                             << std::endl;
        }
    };
    static bool flag=false;
    assert(pristine_deltas.size() == bond_deltas.size());
    for (int id_a : graph.vertices()) // atom ids
	{
		if (not flag)
		{
			//output_atom_into(id_a);
		}
		auto type_a = structure.types[id_a];

		sp2::graph::ud_edge_t edge1{-1, -1, -1},
			      edge2{-1, -1, -1},
			      edge3{-1, -1, -1};

		assert(graph.degree(id_a)==3);
		for (sp2::graph::ud_edge_t edge : graph.edges(id_a))
		{
			if (type_a == structure.types[edge.b])
			{
				// same pucker
				if (edge1.id == -1)
				{
					edge1 = edge;
				}
				else
				{
					edge2 = edge;
				}
			}
			else
			{
				// different pucker
				edge3 = edge;
			}
		}

		// pristine
		vec3_t r12 = vec3_t(0, 0, 0),
			   r13 = vec3_t(0, 0, 0),
			   r14 = vec3_t(0, 0, 0);

		// deformed
		vec3_t r12_def = vec3_t(0, 0, 0),
			   r13_def = vec3_t(0, 0, 0),
			   r14_def = vec3_t(0, 0, 0);

		if (edge1.id != -1)
		{
			r12 = pristine_deltas[edge1.id];
			r12_def = bond_deltas[edge1.id];
		}

		if (edge2.id != -1)
		{
			r13 = pristine_deltas[edge2.id];
			r13_def = bond_deltas[edge2.id];
		}

		if (edge3.id != -1)
		{
			r14 = pristine_deltas[edge3.id];
			r14_def = bond_deltas[edge3.id];
		}

		double d12= r12.mag(),
		       d12_def= r12_def.mag();

		double d13= r13.mag(),
		       d13_def= r13_def.mag();

		double d14= r14.mag(),
		       d14_def= r14_def.mag();

		vec3_t dv1 = r12_def - r12,
			   dv2 = r13_def - r13,
			   dv3 = r14_def - r14;

		double drs1 = dv1.mag_sq()/ds,
			   drs2 = dv2.mag_sq()/ds,
			   drs3 = dv3.mag_sq()/ds;

		double dr1 = dv1.mag()/d,
			   dr2 = dv2.mag()/d,
			   dr3 = dv3.mag()/d;

		double cos_213 = 0.0,
		       cos_214 = 0.0,
		       cos_314 = 0.0;

		double cos_213_def = 0.0,
		       cos_214_def = 0.0,
			   cos_314_def = 0.0;

		double sin_213 = 0.0,
		       sin_214 = 0.0,
		       sin_314 = 0.0;

		if (d12!=0 and d13!=0)
		{
			cos_213 = dot(r12,r13)/(d12*d13);
			cos_213_def = dot(r12_def,r13_def)/(d12_def*d13_def);
			vec3_t numer213 = cross(r12,r13);
			sin_213 = numer213.mag()/(d12*d13);
		}

		if (d12!=0 and d14!=0)
		{
			cos_214 = dot(r12,r14)/(d12*d14);
			cos_214_def = dot(r12_def,r14_def)/(d12_def*d14_def);
			vec3_t numer214 = cross(r12,r14);
			sin_214 = numer214.mag()/(d12*d14);
		}

		if (d13!=0 and d14!=0)
		{
			cos_314 = dot(r13,r14)/(d13*d14);
			cos_314_def = dot(r13_def,r14_def)/(d13_def*d14_def);
			vec3_t numer314=cross(r13,r14);
			sin_314 = numer314.mag()/(d13*d14);
		}

		double d_theta213=0.0;
		if (sin_213!=0.0) {d_theta213= -1.0*(cos_213_def-cos_213)/sin_213;}

		double d_theta214=0.0;
		if (sin_214!=0.0) {d_theta214= -1.0*(cos_214_def-cos_214)/sin_214;}

		double d_theta314= 0.0;
		if (sin_314!=0.0) {d_theta314= -1.0*(cos_314_def-cos_314)/sin_314;}

		// Force Calculation

		///////////////// TERM 1 //////////////////

		vec3_t t1_a=-0.5*Kr*dv1;
		vec3_t t1_b=-0.5*Kr*dv2;

		// potential+=0.25*Kr*ds*(drs1+drs2);

		// forces[id_a]+=t1_a;
		// forces[edge1.b]-=t1_a;
		// forces[id_a]+=t1_b;
		// forces[edge2.b]-=t1_b;

		//////////////// TERM 2 ///////////////////
		vec3_t t2=-1.0*KrP*(dv3);

		// potential+=0.5*ds*KrP*drs3;

		// forces[id_a]+=t2;
		// forces[edge3.b]-=t2;

		////////// Calculate the Gradient of all the Angles ////////////////

		vec3_t grad213=vec3_t(0,0,0);
		vec3_t grad312=vec3_t(0,0,0);
		vec3_t grad214=vec3_t(0,0,0);
		vec3_t grad412=vec3_t(0,0,0);
		vec3_t grad314=vec3_t(0,0,0);
		vec3_t grad413=vec3_t(0,0,0);

		if (r12_def.mag()!=0 and r13_def.mag()!=0 and sin_213!=0)
		{
			grad213=(1.0/sin_213)*(r13_def/(r12_def.mag()*r13_def.mag())-cos_213_def*(r12_def/r12_def.mag_sq()));
			grad312=(1.0/sin_213)*(r12_def/(r12_def.mag()*r13_def.mag())-cos_213_def*(r13_def/r13_def.mag_sq()));
		}

		if (r12_def.mag()!=0 and r14_def.mag()!=0 and sin_214!=0)
		{	
			grad214=(1.0/sin_214)*(r14_def/(r12_def.mag()*r14_def.mag())-cos_214_def*(r12_def/r12_def.mag_sq()));
			grad412=(1.0/sin_214)*(r12_def/(r12_def.mag()*r14_def.mag())-cos_214_def*(r14_def/r14_def.mag_sq()));
		}

		if (r13_def.mag()!=0 and r14_def.mag()!=0 and sin_314!=0)
		{
			grad314=(1.0/sin_314)*(r14_def/(r13_def.mag()*r14_def.mag())-cos_314_def*(r13_def/r13_def.mag_sq()));
		}
		if (r13_def.mag()!=0 and r14_def.mag()!=0 and sin_314!=0)
		{
			grad413=(1.0/sin_314)*(r13_def/(r13_def.mag()*r14_def.mag())-cos_314_def*(r14_def/r14_def.mag_sq()));
		}

		////////// Term 3 /////////////////

		vec3_t t3_a=vec3_t(0,0,0);
		vec3_t t3_b=vec3_t(0,0,0);

		t3_a=1.0*Kt*ds*d_theta213*grad213;
		t3_b=1.0*Kt*ds*d_theta213*grad312;

		// potential+=0.5*Kt*ds*d_theta213*d_theta213;

		// forces[id_a]+=(t3_a+t3_b);
		// forces[edge1.b]-=t3_a;
		// forces[edge2.b]-=t3_b;

		/////////// Term 4 ///////////////

		vec3_t t4_a=vec3_t(0,0,0);
		vec3_t t4_b=vec3_t(0,0,0);
		vec3_t t4_c=vec3_t(0,0,0);
		vec3_t t4_d=vec3_t(0,0,0);

		t4_a=KtP*ds*d_theta214*grad214;
		t4_b=KtP*ds*d_theta214*grad412;

		t4_c=KtP*ds*d_theta314*grad314;
		t4_d=KtP*ds*d_theta314*grad413;

		// potential+=0.5*ds*KtP*(d_theta214*d_theta214+d_theta314*d_theta314);

		// forces[id_a]+=(t4_a+t4_b+t4_c+t4_d);
		// forces[edge1.b]-=(t4_a);
		// forces[edge2.b]-=(t4_c);
		// forces[edge3.b]-=(t4_b+t4_d);

		/////////////// TERM 5 ///////////////////////

		vec3_t t5_a = vec3_t(0.0,0.0,0.0);
		vec3_t t5_b = vec3_t(0.0,0.0,0.0);
		if (dr1!=0.0)
		{
			t5_a=-0.5*Krrp*(dr2*(dv1/dr1));
		}
		if (dr2 !=0.0)
		{
			t5_b=-0.5*Krrp*(dr1*(dv2/dr2));
		}

		// potential+=0.5*ds*Krrp*dr1*dr2;

		// forces[id_a]+=(t5_a+t5_b);
		// forces[edge1.b]-=t5_a;
		// forces[edge2.b]-=t5_b;

		////////////// TERM 6 /////////////////

		vec3_t t6_a = vec3_t(0.0,0.0,0.0);
		vec3_t t6_b = vec3_t(0.0,0.0,0.0);
		vec3_t t6_c = vec3_t(0.0,0.0,0.0);

		if (dr1!=0.0 and dr3!=0)
		{
			t6_a=-0.5*KrrpP*dr3*(dv1/dr1);
			t6_c+=-0.5*KrrpP*dr1*(dv3/dr3);
		}

		if (dr2!=0 and dr3!=0)
		{
			t6_b=-0.5*KrrpP*dr3*(dv2/dr2);
			t6_c+=-0.5*KrrpP*dr2*(dv3/dr3);
		}

		// potential+=0.5*ds*KrrpP*dr2*dr3;
		// potential+=0.5*ds*KrrpP*dr1*dr3;

		// forces[id_a]+=(t6_a+t6_b+t6_c);
		// forces[edge1.b]-=t6_a;
		// forces[edge2.b]-=t6_b;
		// forces[edge3.b]-=t6_c; 

		///////////// Term 7 /////////////////

		vec3_t t7_a = vec3_t(0.0,0.0,0.0);
		vec3_t t7_b = vec3_t(0.0,0.0,0.0);
		vec3_t t7_c = vec3_t(0.0,0.0,0.0);
		vec3_t t7_d = vec3_t(0.0,0.0,0.0);
		vec3_t t7_e = vec3_t(0.0,0.0,0.0);
		vec3_t t7_f = vec3_t(0.0,0.0,0.0);

		t7_b=0.5*ds*Krt*dr1*grad213;
		t7_e=0.5*ds*Krt*dr2*grad213;
		t7_c=0.5*ds*Krt*dr2*grad312;
		t7_f=0.5*ds*Krt*dr1*grad312;

		if (dr1!=0)
		{
			t7_a=-0.5*Krt*d_theta213*(dv1/dr1);
		}

		if(dr2!=0)
		{
			t7_d=-0.5*Krt*d_theta213*(dv2/dr2);
		}

		// potential+=0.5*ds*Krt*d_theta213*(dr1+dr2);
		
		// forces[id_a]+=(t7_a+t7_b+t7_c+t7_d+t7_e+t7_f);
		// forces[edge1.b]-=(t7_a+t7_b+t7_e);
		// forces[edge2.b]-=(t7_c+t7_d+t7_f);

		//////////// Term 8 ////////////////

		// The t8_a# terms refer to interactions between atoms 2-1-4

		vec3_t t8_aa = vec3_t(0.0,0.0,0.0);
		vec3_t t8_ab = vec3_t(0.0,0.0,0.0);
		vec3_t t8_ac = vec3_t(0.0,0.0,0.0);
		vec3_t t8_ad = vec3_t(0.0,0.0,0.0);
		vec3_t t8_ae = vec3_t(0.0,0.0,0.0);
		vec3_t t8_af = vec3_t(0.0,0.0,0.0);

		// the t8_b# terms refer to the interactions between atoms 3-1-4

		vec3_t t8_ba = vec3_t(0.0,0.0,0.0);
		vec3_t t8_bb = vec3_t(0.0,0.0,0.0);
		vec3_t t8_bc = vec3_t(0.0,0.0,0.0);
		vec3_t t8_bd = vec3_t(0.0,0.0,0.0);
		vec3_t t8_be = vec3_t(0.0,0.0,0.0);
		vec3_t t8_bf = vec3_t(0.0,0.0,0.0);

		t8_ab=0.5*ds*KrtP*dr1*grad214;
		t8_ae=0.5*ds*KrtPP*dr3*grad214;
		t8_ac=0.5*ds*KrtPP*dr3*grad412;
		t8_af=0.5*ds*KrtP*dr1*grad412;

		t8_bb=0.5*ds*KrtP*dr2*grad314;
		t8_be=0.5*ds*KrtPP*dr3*grad314;
		t8_bc=0.5*ds*KrtPP*dr3*grad413;
		t8_bf=0.5*ds*KrtP*dr2*grad413;

		if (dr1!=0)
		{
			t8_aa=-0.5*KrtP*d_theta214*(dv1/dr1);
		}

		if (dr2!=0)
		{
			t8_ba=-0.5*KrtP*d_theta314*(dv2/dr2);
		}

		if(dr3!=0)
		{
			t8_ad=-0.5*KrtPP*d_theta214*(dv3/dr3);
			t8_bd=-0.5*KrtPP*d_theta314*(dv3/dr3);
		}

		potential+=0.5*ds*((KrtP*dr1+KrtPP*dr3)*d_theta214+(KrtP*dr2+KrtPP*dr3)*d_theta314);
		
		forces[id_a]+=(t8_aa+t8_ab+t8_ac+t8_ad+t8_ae+t8_af);
		forces[id_a]+=(t8_ba+t8_bb+t8_bc+t8_bd+t8_be+t8_bf);
		forces[edge1.b]-=(t8_aa+t8_ab+t8_ae);
		forces[edge3.b]-=(t8_ac+t8_ad+t8_af);
		forces[edge2.b]-=(t8_ba+t8_bb+t8_be);
		forces[edge3.b]-=(t8_bc+t8_bd+t8_bf);
	}
	if(not flag)
	{
		for (int i=0;i<forces.size();i++)
		{
			// std::cout<<i<<" "<<forces[i].x()<<" "<<forces[i].y()<<" "<<forces[i].z()<<std::endl;
		}
		flag=true;
	}
}