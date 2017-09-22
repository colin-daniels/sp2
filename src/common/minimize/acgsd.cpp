#include "common/minimize/minimize.hpp"
#include "common/math/numerical_diff.hpp"
#include "common/math/blas.hpp"
#include "settings.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include <optional>

using namespace std;
using namespace sp2;

// Helper types used in acgsd to group pieces of changing state.
namespace {

    // Helper type to collect position and output together, as they
    // are frequently used together.
    struct point_t {
        vector<double> position;
        double value;
        vector<double> gradient;

        point_t(vector<double> pos, double val, vector<double> grad)
            : position(move(pos)), value(val), gradient(move(grad))
        {}
    };

    // These are values at the "fencepost" that lies between each iteration.
    // They are all updated at once to assist in reasoning about the current
    // state of the algorithm.
    struct saved_t : public point_t
    {
        // this additionally has the fields from point_t (superclass)
        double alpha; // step size

        saved_t(double a, point_t point) : point_t(point), alpha(a) {}
        const point_t& point() const { return *this; }
    };

    // These are values that describe the events of the previous iteration.
    // They are invalid during the first iteration.
    struct last_t
    {
        vector<double> direction;   // direction searched (normalized)

        // NOTE: These next three are all zero when linesearch has failed.
        //       This can be a problem for d_value in particular.
        double d_value;             // change in value
        vector<double> d_position;  // change in position
        vector<double> d_gradient;  // change in gradient

        bool ls_failed;             // linesearch failed?

        last_t(vector<double> dir, double d_val, vector<double> d_pos,
            vector<double> d_grad, bool ls_fail)
            : direction(move(dir)), d_value(d_val), d_position(move(d_pos)),
            d_gradient(move(d_grad)), ls_failed(ls_fail)
        {}
    };
}

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

    auto compute_point = [&](const vector<double>& position) {
        double value;
        vector<double> gradient;
        std::tie(value, gradient) = objective_fn(position);
        return point_t{position, value, gradient};
    };

////////////////////////////////////////////////////////////////////////////////
// Loop start                                                                 //
////////////////////////////////////////////////////////////////////////////////

    // These are all updated only at the end of an iteration.

    // Describes data at the "fencepost" between iterations.
    saved_t saved(1.0, compute_point(initial_position));
    // Describes the previous iteration
    boost::optional<last_t> last;
    // Number of elapsed iterations
    int iterations = 0;

    while (true)
    {

////////////////////////////////////////////////////////////////////////////////
// Closures for things commonly done during an iteration                      //
////////////////////////////////////////////////////////////////////////////////

        // Compute at a position relative to saved.position
        auto compute_in_dir = [&](double alpha, vector<double> direction) {
            std::vector<double> position = saved.position;
            vaxpy(alpha, direction, position);

            return compute_point(position);
        };

        auto warning = [&](std::string msg, double alpha, const point_t &point,
            std::function<void(ostream&)> additional_output = {})
        {
            if (settings.output_level > 0) {
                cerr << msg << '\n'
                     << "Iterations: " << iterations << '\n'
                     << "     Alpha: " << alpha << '\n'
                     << "     Value: " << point.value << '\n'
                     << " Grad Norm: " << vmag(point.gradient) << '\n';
                if (additional_output)
                    additional_output(cerr);
            }
        };

        // use as 'return fatal(...);'
        // this will return or throw based on the 'except_on_fail' setting
        auto fatal = [&](std::string msg, double alpha, const point_t &point,
                std::function<void(ostream&)> additional_output = {})
        {
            warning("ACGSD Failed: " + msg, alpha, point, additional_output);
            if (settings.except_on_fail)
                throw runtime_error(msg);
            return point.position;
        };

////////////////////////////////////////////////////////////////////////////////
// Per-iteration output                                                       //
////////////////////////////////////////////////////////////////////////////////

        if (settings.output_level > 2)
        {
            double d_value = (last ? last->d_value : 0.0);
            double grad_mag = vmag(saved.gradient);
            cout << " i: " << setw(6) << iterations
                 << "  v: " << setw(18) << setprecision(14) << saved.value
                 << " dv: " << setw(13) << setprecision(7) << d_value
                 << "  g: " << setw(13) << setprecision(7) << grad_mag
                 << endl;
        }

        // call the output function if applicable
        if (settings.intermediate_output_interval > 0 &&
            iterations % settings.intermediate_output_interval == 0)
            settings.output_fn(saved.position);

////////////////////////////////////////////////////////////////////////////////
// Evaluate exit conditions                                                   //
////////////////////////////////////////////////////////////////////////////////

        { // scope
            double grad_mag = vmag(saved.gradient);
            double grad_max = max_norm(saved.gradient);

            bool done = false;
            done |= (grad_mag < settings.gradient_tolerance);
            done |= (grad_max < settings.grad_max_tolerance);
            if (last && !last->ls_failed)
                done |= (std::abs(last->d_value) < settings.value_tolerance);
            if (settings.iteration_limit > 0)
                done |= (iterations >= settings.iteration_limit);

            if (done)
            {
                // output final info
                if (settings.output_level > 1)
                {
                    double d_value = last ? last->d_value : 0.0;
                    cout << "ACGSD Finished.\n"
                         << "Iterations: " << iterations << '\n'
                         << "     Value: " << saved.value << '\n'
                         << " Delta Val: " << d_value << '\n'
                         << " Grad Norm: " << grad_mag << '\n'
                         << "  Grad Max: " << grad_max << '\n';
                }
                return saved.position;
            }
        } // scope

////////////////////////////////////////////////////////////////////////////////
// Calculate the search direction.                                            //
////////////////////////////////////////////////////////////////////////////////

        // whether or not we will use steepest descent
        // note: force us to use steepest descent if linesearch failed
        // last iteration
        //
        // an IIFE is used for complex control flow
        vector<double> direction = ([&] {

            // Consider the direction  'beta * dx - g'
            if (last && !last->ls_failed)
            {
                double beta = calc_beta_acgsd(
                    saved.gradient, last->d_position, last->d_gradient);

                auto direction = saved.gradient;
                vscal(-1.0, direction);
                vaxpy(beta, last->d_position, direction);

                // use this direction unless it is almost directly uphill
                if (!should_revert_acgsd(saved.gradient, direction))
                    return direction;
            }

            // Fallback to steepest descent:  '-g'
            if (settings.output_level > 2)
                cout << "Using steepest descent." << endl;

            auto direction = saved.gradient;
            vscal(-1.0, direction);
            return direction;
        })();

        // NOTE: The original source scaled alpha instead of normalizing
        //       direction, which seems to be a fruitless optimization
        //       that only serves to amplify the mental workload.
        if (!vnormalize(direction))
            return fatal("cannot normalize direction", saved.alpha, saved);

////////////////////////////////////////////////////////////////////////////////
// Perform the linesearch.                                                    //
////////////////////////////////////////////////////////////////////////////////

        // (these cache the best computation by linesearch)
        double ls_alpha = 0;
        auto ls_point = saved.point();

        // Linesearch along direction for a better point.
        // It is possible that no better point will be found, in which case
        //  the displacement returned will naturally be zero.
        auto next_alpha = linesearch(
            settings.linesearch_settings,
            saved.alpha,
            diff1d_fn_t([&](double alpha) {
                point_t point = compute_in_dir(alpha, direction);
                double slope = vdot(point.gradient, direction);

                // update cache, checking values to predict which
                //  point linesearch will prefer to use.
                // (future additions to linesearch may make this less reliable)
                if (point.value < ls_point.value) {
                    ls_alpha = alpha;
                    ls_point = point;
                }

                if (last && last->ls_failed && settings.output_level > 0)
                {
                    cout << "LS: a: " << setprecision(14) << alpha
                         << "\tv: " << setprecision(14) << point.value
                         << "\ts: " << setprecision(14) << slope << endl;
                }
                return make_pair(point.value, slope);
            }
        ));
        point_t next_point = (ls_alpha == next_alpha)
                             ? ls_point // extraneous computation avoided!
                             : compute_in_dir(next_alpha, direction);

        // if the linesearch failed, note it and try
        //  again next iteration with steepest descent
        bool ls_failed = (next_alpha == 0);
        if (ls_failed)
        {
            if (last && last->ls_failed)
            {
                return fatal(
                    "linesearch failure (second)", saved.alpha, saved.point(),
                    [&](ostream &out) {
                        auto test_dir = saved.gradient;
                        vnormalize(test_dir);

                        double numerical = util::central_difference<9, double>(
                            [&](double a) {
                                return compute_in_dir(a, test_dir).value;
                            }, 0.0, 1e-3);

                        out << "Numerical gradient: "
                            << setprecision(14) << numerical << endl;
                    }
                );
            }
            else
                warning("Linesearch failure, switching to steepest descent",
                    saved.alpha, saved.point());
        }

////////////////////////////////////////////////////////////////////////////////
// Update quantities following the linesearch.                                //
////////////////////////////////////////////////////////////////////////////////

        // NOTE: If linesearch failed, then next_alpha was zero and so
        //       next_point is saved.point(). Hence, all of these will be zero.
        auto d_value = next_point.value - saved.value;
        auto d_position = next_point.position;
        auto d_gradient = next_point.gradient;
        vaxpy(-1.0, saved.position, d_position);
        vaxpy(-1.0, saved.gradient, d_gradient);

        last = last_t(direction, d_value, d_position, d_gradient, ls_failed);
        saved = ls_failed ? saved
                          : saved_t{next_alpha, next_point};
        iterations++;
    }
    // unreachable
}
