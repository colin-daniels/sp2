#include <vector>
#include "vec3_t.hpp"
#include "ud_graph_t.hpp"

using namespace sp2;

struct system_state
{
	std::vector<int> atom_types; // indexed by atom/vertex ids

	sp2::graph::ud_graph_t graph;
	std::vector<vec3_t> bond_deltas; // indexed by bond/edge ids
	std::vector<vec3_t> pristine_deltas; // indexed by bond/edge ids
};

//vector<vec3_t> deltas;

//deltas.push_back(vec3_t(x, y, z));
//deltas.emplace_back(x, y, z);

vec3_t sin_grad (vec3_t r1, vec3_t r2)
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

vec3_t cos_grad (vec3_t r1, vec3_t r2)
{
	vec3_t num1=(r1+r2)*r1.mag()*r2.mag();
	vec3_t num2 =dot(r1,r2)*(r1*(r2.mag()/r1.mag())+r2*(r1.mag()/r2.mag()));

	vec3_t grad=(num1-num2)/(r1.mag_sq()*r2.mag_sq());

	return grad;
}

double calculate(system_state &state, std::vector<vec3_t> &forces)
{

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
		ds=d*d;

	//std::vector<vec3_t> forces(atom_types.size(), vec3_t(0, 0, 0));

	double sum=0.0;
	for (int id_a : state.graph.vertices()) // atom ids
	{
		int type_a = state.atom_types[id_a];

		// // if possible
		// for (ud_edge_t edge : graph.edges(atom_id))
		// {
		// 	int id_b = edge.b;
		// 	int type_b = atom_types[id_b];
		//
		// 	bool same_pucker = (type_a == type_b);
		//
		//
		// 	vec3_t delta = bond_deltas[edge.id];
		// 	const vec3_t pristine_delta = pristine_deltas[edge.id],
		// 				 deformation_delta = delta - pristine_delta;
		// }

		sp2::graph::ud_edge_t edge1{-1, -1, -1},
			      edge2{-1, -1, -1},
			      edge3{-1, -1, -1};

		for (sp2::graph::ud_edge_t edge : state.graph.edges(id_a))
		{
			if (type_a == state.atom_types[edge.b])
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

		// vec3_t example_vec;
		// vec3_t unit_vec = example_vec.unit_vector();
		// example_vec.normalize();
		// E = || r ||;
		// F = - r.unit_vector();
		//   = - r / r.mag();

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
			r12 = state.pristine_deltas[edge1.id];
			r12_def = state.bond_deltas[edge1.id];
			// edge1.b 
		}

		if (edge2.id != -1)
		{
			r13 = state.pristine_deltas[edge2.id];
			r13_def = state.bond_deltas[edge2.id];
		}

		if (edge3.id != -1)
		{
			r14 = state.pristine_deltas[edge3.id];
			r14_def = state.bond_deltas[edge3.id];
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

		double drs1 = dv1.mag_sq(),
			   drs2 = dv2.mag_sq(),
			   drs3 = dv3.mag_sq();

		double dr1 = dv1.mag(),
			   dr2 = dv2.mag(),
			   dr3 = dv3.mag();

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

		// Deformation Energy term calculations
		sum+= 0.5*Kr*(drs1+drs2); // Term 1
		sum+= 0.5*KrP*drs3; // Term2
		sum+= ds*Kt*d_theta213*d_theta213; // Term3
		sum+= ds*KtP*(d_theta214*d_theta214+d_theta314*d_theta314); // Term4
		sum+= Krrp*dr1*dr2; // Term5
		sum+= KrrpP*(dr1+dr2)*dr3; // Term6
		sum+= d*Krt*(dr1+dr2)*d_theta213; // Term7
		sum+= d*((KrtP*dr1+KrtPP*dr3)*d_theta214+(KrtP*dr2+KrtPP*dr3)*d_theta314); // Term8

		vec3_t t1=-1.0*Kr*(dv1+dv2);
		vec3_t t2=-1.0*KrP*dv3;

		// Term 3 of the Sum
		vec3_t grad_cos213_def = cos_grad(r12_def,r13_def);
		vec3_t grad_cos213 = cos_grad(r12,r13);
		vec3_t grad_sin213 = sin_grad(r12,r13);
		vec3_t d_dtheta213 = (-1.0*(grad_cos213_def-grad_cos213)*sin_213-(cos_213_def-cos_213)*grad_sin213)/(sin_213*sin_213);
		
		vec3_t t3=-2.0*ds*Kt*d_theta213*d_dtheta213;

		// Term 4 of the Sum
		vec3_t grad_cos214_def = cos_grad(r12_def,r14_def);
		vec3_t grad_cos214 = cos_grad(r12,r14);
		vec3_t grad_sin214 = sin_grad(r12,r14);
		vec3_t d_dtheta214 = (-1.0*(grad_cos214_def-grad_cos214)*sin_214-(cos_214_def-cos_214)*grad_sin214)/(sin_214*sin_214);

		vec3_t grad_cos314_def = cos_grad(r13_def,r14_def);
		vec3_t grad_cos314 = cos_grad(r13,r14);
		vec3_t grad_sin314 = sin_grad(r13,r14);
		vec3_t d_dtheta314 = (-1.0*(grad_cos314_def-grad_cos314)*sin_314-(cos_314_def-cos_314)*grad_sin314)/(sin_314*sin_314);
		
		vec3_t t4=-2.0*ds*KtP*(d_theta214*d_dtheta214+d_theta314*d_dtheta314);

		// Term 5 of the Sum
		vec3_t t5=-1.0*Krrp*(dv1*(dr2/dr1)+dv2*(dr1/dr2));

		// Term 6 of the Sum
		vec3_t t6=-1.0*(KrrpP*(dv1*(dr3/dr1)+dv3*(dr1/dr3))+KrrpP*(dv2*(dr3/dr2)+dv3*(dr2/dr3)));

		// Term 7 of the Sum
		vec3_t t7=-1.0*Krt*d*(d_theta213*((dv1/dr1)+(dv2/dr2))+(dr1+dr2)*d_dtheta213);

		// Term 8 of the Sum
		vec3_t t8_one=d*(d_theta214*(KrtP*(dv1/dr1)+KrtPP*(dv3/dr3))+(KrtP*dr1+KrtPP*dr3)*d_dtheta214);
		vec3_t t8_two=d*(d_theta314*(KrtP*(dv2/dr2)+KrtPP*(dv3/dr3))+(KrtP*dr2+KrtPP*dr3)*d_dtheta314);
		vec3_t t8=-1.0*(t8_one+t8_two);
		
	}
	return 0.5*sum;
}
