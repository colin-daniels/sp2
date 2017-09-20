#include <common/math/vec3_util.hpp>
#include <common/math/blas.hpp>
#include "common/minimize/minimize.hpp"

std::vector<double> sp2::minimize::fire(sp2::diff_fn_t grad_fn,
    double mass, std::vector<double> initial_position,
    const fire_settings_t &settings)
{
    std::size_t n = initial_position.size();

    double dt = settings.dt_initial,
        alpha = settings.alpha_initial;

    int counter = 0;

    // initialize velocity and force to zero
    std::vector<double> velocity(n, 0),
        force(n, 0);

    // make a copy of x to track previous values for verlet
    std::vector<double> x = initial_position;

    for (int k = 0; k < settings.max_iter; ++k)
    {
        // update position from velocity verlet
        for (std::size_t i = 0; i < n; ++i)
            x[i] += velocity[i] * dt + 0.5 * (force[i] / mass)  * (dt * dt);

        // get force from gradient at x(t + dt)
        auto force_prev = force;
        force = grad_fn(x).second;
        vscal(-1.0, force);


        // update velocity from velocity verlet
        for (std::size_t i = 0; i < n; ++i)
            velocity[i] += (force_prev[i] + force[i]) / (2 * mass) * dt;

        //check for convergence
        double force_mag = std::sqrt(vdot(force, force));

        // output intermediate something
        settings.output_fn(x);

        if (force_mag < settings.grad_tol)
            break;


        // Calculate P = F * v
        double P = vdot(force, velocity);
        double velocity_mag = std::sqrt(vdot(velocity, velocity));

        // set v = (1 - alpha) * v + alpha * (F / |F|) * |v|
        for (std::size_t i = 0; i < n; ++i)
            velocity[i] = (1 - alpha) * velocity[i]
                        + alpha * (force[i] / force_mag) * velocity_mag;

        std::cout <<   "i: "  << k
                  << "\ta: "  << alpha
                  << "\tdt: " << dt
                  << "\tv: "  << velocity_mag
                  << "\tf: "  << force_mag
                  << "\tcount: " << counter
                  << "\tP: "  << P
                  << std::endl;

        if (P <= 0)
        {
            // freeze the system if we're going uphill
            for (auto &v : velocity)
                v = 0;

            // decrease time step
            dt *= settings.dt_mult_decrease;

            // reset alpha and counter
            alpha = settings.alpha_initial;
            counter = 0;
        }
        else if (++counter > settings.N_min)
        {
            // increase time step and alpha
            dt = std::max(
                settings.dt_max,
                dt * settings.dt_mult_increase
            );

            alpha *= settings.alpha_mult;
        }
    }

    std::cout << "Done." << std::endl;
    return x;
}
