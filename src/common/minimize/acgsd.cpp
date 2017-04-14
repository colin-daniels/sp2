#include "common/minimize/minimize.hpp"
#include "common/math/numerical_diff.hpp"
#include "common/math/blas.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace sp2;

/// calculate beta (ACGSD)
double calc_beta_acgsd(const vector<double> &gradient,
    const vector<double> &delta_x, const vector<double> &delta_g)
{
    // inner products needed for beta
    double dg_dx = vdot(delta_g, delta_x), // dg.dx
        g_dg = vdot(gradient, delta_g),    // g.dg
        g_dx = vdot(gradient, delta_x);    // g.dx

    // beta = (g.dg / dg.dx) * (1.0 - (g.dx / dg.dx))
    return (g_dg / dg_dx) * (1.0 - g_dx / dg_dx);
}

/// whether the search should revert to steepest descent (ACGSD)
bool should_revert_acgsd(const vector<double> &gradient,
    const vector<double> &direction)
{
    double dots = sqrt(vdot(gradient, gradient) * vdot(direction, direction));
    // revert if: g.d > -0.001 * sqrt(d.d * g.g), basically just when the search
    // direction is almost parallel to (or pointing the same was as) the gradient
    return vdot(gradient, direction) > -1e-3 * dots;
}

std::vector<double> minimize::acgsd(scalar_fn_t objective_fn,
    vector_fn_t gradient_fn, std::vector<double> initial_position,
    const acgsd_settings_t &settings)
{
    return acgsd([&](const auto &v) {
        return make_pair(objective_fn(v), gradient_fn(v));
    }, initial_position, settings);
}

std::vector<double> minimize::acgsd(diff_fn_t objective_fn,
    std::vector<double> initial_position,
    const acgsd_settings_t &settings)
{
    if (settings.value_tolerance    <= 0 &&
        settings.grad_max_tolerance <= 0 &&
        settings.gradient_tolerance <= 0 &&
        settings.iteration_limit    <= 0)
        throw invalid_argument("all exit conditions disabled or set to 0.");

    // storage variables:
    int iter = 0,             // current iteration
        n_target_fail = 0;    // number of times the targeted minimization
    // exit condition has been failed in a row

    bool ls_failed = false;   // whether or not the linesearch has failed so far

    double alpha = 1.0,       // current step size
        ddot_last = 1.0,       // ||direction|| from the previous iteration
        value;                 // current objective function value

    vector<double> direction, // current search direction
        position,             // current position
        gradient;             // current gradient

    vector<double> delta_x,   // change in position between iterations
        delta_g;              // change in gradient between iterations

    // information output function for errors
    auto output_info = [&]() {
        cout << "Iteration: " << iter  << '\n'
             << "    Alpha: " << alpha << '\n'
             << "    Value: " << value << '\n'
             << "Grad Norm: " << sqrt(vdot(gradient, gradient))
             << endl;
    };


    do {
////////////////////////////////////////////////////////////////////////////////
// Calculate the search direction.                                            //
////////////////////////////////////////////////////////////////////////////////

        // whether or not we will use steepest descent
        bool use_steepest;

        if (iter == 0)
        {
            // initialize quantities
            position = initial_position;
            tie(value, gradient) = objective_fn(position);

            // if its the first iteration, we are forced to use steepest descent
            use_steepest = true;
        }
        else
        {
            // get beta (ACGSD)
            double beta = calc_beta_acgsd(gradient, delta_x, delta_g);

            // calculate the search direction = beta * dx - g
            direction = gradient;
            vscal(-1.0, direction);
            vaxpy(beta, delta_x, direction);

            // check if we should revert to using steepest descent
            use_steepest = should_revert_acgsd(gradient, direction);
        }

        // note: force us to use steepest descent if linesearch failed
        // last iteration
        if (use_steepest || ls_failed)
        {
            if (settings.output_level > 2)
                cout << "Using steepest descent." << endl;

            // descent direction is just the negative gradient
            direction = gradient;
            vscal(-1.0, direction);
        }

////////////////////////////////////////////////////////////////////////////////
// Record old data and calculate the next linesearch step size (alpha)        //
////////////////////////////////////////////////////////////////////////////////

        // record old data
        auto val_last = value,
            ddot = vdot(direction, direction);

        delta_x = position;
        delta_g = gradient;

        // call the output function if applicable
        if (settings.intermediate_output_interval > 0 &&
            iter % settings.intermediate_output_interval == 0)
            settings.output_fn(position);

        // adjust and check alpha for the next linesearch
        alpha = alpha * sqrt(ddot_last / ddot);

        // store old ddot
        ddot_last = ddot;

        // check alpha
        if (!std::isfinite(alpha))
        {
            if (settings.output_level > 0)
            {
                cout << "Error, non-finite alpha." << endl;
                output_info();
            }

            if (settings.except_on_fail)
                throw runtime_error("non-finite alpha");

            return position;
        }

////////////////////////////////////////////////////////////////////////////////
// Perform the linesearch and check for errors.                               //
////////////////////////////////////////////////////////////////////////////////

        // temporary to avoid re-calculation of the value post-linesearch
        auto ls_position = position;

        // 1D function that will be passed to the linsearch
        auto ls_function = [&](double alpha_in) {
            // set the new positon
            ls_position = position;
            vaxpy(alpha_in, direction, ls_position);

            // get the value/gradient (note: updated in acgsd by reference)
            tie(value, gradient) = objective_fn(ls_position);
            // slope by just dotting the direction with the gradient
            auto slope = vdot(gradient, direction);

            if (ls_failed && settings.output_level > 0)
            {
                cout << "LS: a: " << setprecision(14) << alpha_in
                     << "\tv: " << setprecision(14) << value
                     << "\ts: " << setprecision(14) << slope << endl;
            }
            return make_pair(value, slope);
        };

        // do the linesearch
        auto new_alpha = linesearch(ls_function, alpha);

        // get updated position/value/gradient
        vaxpy(new_alpha, direction, position);

        // only need to update value/grad if the last evaluation in
        // the linesearch was not the same as final alpha returned
        if (ls_position != position)
            tie(value, gradient) = objective_fn(position);

        // if the linesearch failed, note it and try one more iteration
        if (new_alpha == 0 && !ls_failed)
        {
            ls_failed = true;
            if (settings.output_level > 0)
            {
                cout << "Linesearch failure (first), switching "
                     << "to steepest descent." << endl;
                output_info();
            }
        }
        else if (new_alpha == 0)
        {
            if (settings.output_level > 0)
            {
                cout << "Linesearch failure (second), aborting."  << endl;
                output_info();

                auto test_dir = gradient;
                vscal(1.0 / static_cast<double>(sqrt(vdot(test_dir, test_dir))),
                    test_dir);

                double numerical_grad = util::central_difference<9, double>(
                    [&](double a) {
                        auto temp = position;
                        vaxpy(a, test_dir, temp);
                        return objective_fn(temp).first;
                    }, 0.0, 1e-3);

                cout << "Numerical gradient: "
                     << setprecision(14) << numerical_grad << endl;
            }

            if (settings.except_on_fail)
                throw runtime_error("linesearch failure");

            return position;
        }
        else
        {
            ls_failed = false;
            alpha = new_alpha;
        }

////////////////////////////////////////////////////////////////////////////////
// Update quantities following the linesearch.                                //
////////////////////////////////////////////////////////////////////////////////

        // finish calculating delta_x and delta_g
        vaxpy(-1.0, position, delta_x);
        vaxpy(-1.0, gradient, delta_g);
        vscal(-1.0, delta_x);
        vscal(-1.0, delta_g);

        // output current progress
        double grad_mag = sqrt(vdot(gradient, gradient)),
            grad_max = max_norm(gradient),
            delta_val = value - val_last;

        if (settings.output_level > 2)
            cout << " i: " << setw(6)  << iter << ' '
                 << " v: " << setw(18) << setprecision(14) << value     << ' '
                 << "dv: " << setw(13) << setprecision(7)  << delta_val << ' '
                 << " g: " << setw(13) << setprecision(7)  << grad_mag  << ' '
                 << endl;

        // targeted minimization
        if (settings.target_ratio_tol > 0 &&
            settings.target_exit_min > 0)
        {
            double delta_target = value - settings.target_value,
                delta_ratio = abs(delta_val / delta_target);

            n_target_fail += 1;
            if (delta_target < 0 || delta_ratio > settings.target_ratio_tol)
                n_target_fail = 0;
        }

        // exit conditions
        bool acgsd_exit = abs(delta_val) < settings.value_tolerance    ||
                                grad_mag < settings.gradient_tolerance ||
                      max_norm(gradient) < settings.grad_max_tolerance ||
                           n_target_fail > settings.target_exit_min    ||
            (settings.iteration_limit > 0 && iter >= settings.iteration_limit);

        if (acgsd_exit)
        {
            // output final info
            if (settings.output_level > 1)
                cout << "ACGSD Finished.\n"
                     << "Iteration: " << iter      << '\n'
                     << "    Value: " << value     << '\n'
                     << "Delta Val: " << delta_val << '\n'
                     << "Grad Norm: " << grad_mag  << '\n'
                     << " Grad Max: " << grad_max
                     << endl;

            break;
        }

        iter += 1;
    } while (true);

    return position;
}
