#ifndef SP2_METROPOLIS_HPP
#define SP2_METROPOLIS_HPP

#include "settings.hpp"
#include "metropolis_enums.hpp"
#include "common/structure_t.hpp"
#include "common/python/bindings.hpp"
#include "common/python/environment.hpp"
#include "common/python/utility.hpp"

#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>

namespace sp2 {
namespace minimize {
namespace metropolis {

template<typename P>
using objective_fn_t = std::function<double(const P &)>;

template<typename P>
using mutate_fn_t = std::function<P(const P&)>;

template<typename P, typename D>
using generate_fn_t = std::function<D(const P&)>;

template<typename P, typename D>
using apply_fn_t = std::function<P(const P&, const D&)>;

template<typename P, typename D>
using applied_fn_t = std::function<void(const P&, const D&, double, double)>;

template<typename P, typename D>
using visit_fn_t = std::function<void(const P&, double, bool)>;

template<typename P, typename D>
using scale_fn_t = std::function<D(const P&, const D&, double)>;

template<typename P, typename D>
using is_repeatable_fn_t = std::function<bool(const P&, const D&)>;

/// Advanced callbacks for metropolis.
template<typename P, typename D>
struct callbacks_t
{
    /// Generates an abstract mutation.
    generate_fn_t<P, D> generate;
    /// Applies an abstract mutation to the position vector.
    apply_fn_t<P, D> apply;
    /// Called on each mutation (with the old structure)
    ///  after computing the objective.
    applied_fn_t<P, D> applied;
    /// Called on each structure after computing the objective.
    visit_fn_t<P, D> visit;
};

template<typename P, typename D>
P advanced(
    objective_fn_t<P> objective_fn,
    callbacks_t<P, D> callbacks,
    P initial_position,
    const metropolis_settings_t &settings)
{
    using namespace std;
    if (settings.iteration_limit <= 0 &&
        settings.improve_iteration_limit <= 0)
        throw invalid_argument("all exit conditions disabled or set to 0.");

    int mutations = 0,        // total number of mutations tested
        mutations_since = 0;  // number of mutations since last improvement

    // Initialize
    P position_best = initial_position;
    double value_best = objective_fn(position_best);

    callbacks.visit(position_best, value_best, true);

    if (settings.output_level > 1) {
        cout << " Initial Value: " << value_best << endl;
    }

    do {
        D mutation = callbacks.generate(position_best);

        auto position_new = callbacks.apply(position_best, mutation);
        mutations += 1;
        mutations_since += 1;
        double value_new = objective_fn(position_new);
        bool successful = value_new < value_best;

        callbacks.visit(position_best, value_best, true);
        callbacks.applied(position_best, mutation, value_best, value_new);
        if (successful) {
            if (settings.output_level > 1) {
                cout << "Improved Value: " << setprecision(13) << value_new
                     << "   after " << mutations_since << " mutations"
                     << endl;
            }

            // keep modified position and value
            value_best = value_new;
            position_best = std::move(position_new);
            mutations_since = 0;
        }

        // for now, simply reject all inferior values
        // (this is Metropolis at 0K)

        // exit conditions
        bool done = false;
        done = done || settings.improve_iteration_limit > 0 &&
                       mutations_since >= settings.improve_iteration_limit;

        done = done || settings.iteration_limit > 0 &&
                       mutations >= settings.iteration_limit;

        if (done)
        {
            // output final info
            if (settings.output_level > 1)
                cout << "Metropolis Finished.\n"
                     << "     Value: " << value_best << '\n'
                     << " Mutations: " << mutations << '\n'
                     << endl;

            break;
        }
    } while (true);

    return position_best;
}

// Simplified interface taking a single function that produces a new position.
template <typename P>
P basic(
    objective_fn_t<P> objective_fn,
    mutate_fn_t<P> mutation_fn,
    P initial_position,
    const metropolis_settings_t &settings)
{
    auto callbacks = callbacks_t<P, P> {
        // generate
        [&](const auto &pos) { return mutation_fn(pos); },
        // apply
        [&](const auto &pos, const auto &diff) { return diff; },
        // applied
        [&](const auto &pos, const auto &diff, double was, double now) { return; },
        // visit
        [&](const auto &pos, double value, bool better) { return; }
    };

    return advanced<P,P>(objective_fn,
        callbacks, initial_position, settings);
}

/*----------------------------------------------------------------------------*/
// scaling adapters
//
// This is a set of adapter utilities that extend a set of callbacks with the
//  ability to repeat mutations that appear to be working.

enum class scaling_action {
    SCALE_UP = 0,
    SCALE_DOWN = 1,
    REGENERATE = 2
};

constexpr static const int READY_TO_SCALE_UP = -1;
constexpr static const int READY_TO_GENERATE = 0;
struct scaling_state_t
{
    int downscales_max;

    // Models the algebraic data type
    //     | ReadyToScaleUp      // assigned the value -1
    //     | ReadyToGenerate     // assigned the value 0 for convenient logic
    //     | ReadyToScaleDown(turns_left: Nonzero<unsigned>)
    int value = 0;

    // create from config
    scaling_state_t(int max) : downscales_max(max)
    {}

    scaling_state_t(int max, int value)
        : downscales_max(max), value(value)
    {}

    // change state data while preserving config
    scaling_state_t derive(int value)
    { return {downscales_max, value}; }

    scaling_action recommended_action()
    {
        switch (value)
        {
        case READY_TO_GENERATE:
            return scaling_action::REGENERATE;
        case READY_TO_SCALE_UP:
            return scaling_action::SCALE_UP;
        default:
            return scaling_action::SCALE_DOWN;
        }
    }

    scaling_state_t succeed()
    { return this->derive(READY_TO_SCALE_UP); }

    scaling_state_t fail()
    {
        switch (value)
        {
        case READY_TO_GENERATE:
            return this->derive(READY_TO_GENERATE);
        case READY_TO_SCALE_UP:
            // Enter initial state of READY_TO_SCALE_DOWN
            return this->derive(downscales_max);
        default:
            if (value < 0)
                throw std::logic_error(":V");
            return this->derive(value - 1);
        }
    }
};

/// Adapter around callbacks_t that repeats mutations at different
/// strengths on success.
template<typename P, typename D>
struct scaling_control_t
{
    // config
    callbacks_t<P, D> cbs;
    is_repeatable_fn_t<P, D> is_repeatable;
    scale_fn_t<P, D> scale;

    metropolis_scaling_settings_t settings;

    // state
    scaling_state_t scaling_state;
    std::unique_ptr<D> prev;

    scaling_control_t(
        metropolis_scaling_settings_t settings,
        callbacks_t<P, D> cbs,
        is_repeatable_fn_t<P, D> is_repeatable,
        scale_fn_t<P, D> scale
    )
        : settings(settings), cbs(cbs), scale(scale),
        is_repeatable(is_repeatable),
        scaling_state{settings.downscale_max_attempts}
    {}

    callbacks_t<P, D> get_callbacks()
    {
        return {
            // generate
            [&](const auto &pos) {
                switch (scaling_state.recommended_action())
                {
                case scaling_action::SCALE_DOWN:
                    return scale(pos, *prev, settings.downscale_by);
                case scaling_action::SCALE_UP:
                    return scale(pos, *prev, settings.upscale_by);
                case scaling_action::REGENERATE:
                    return cbs.generate(pos);
                }
            },

            // apply
            cbs.apply,

            // applied
            [&](const auto &pos, const auto &mut, double was, double now) {
                cbs.applied(pos, mut, was, now);
                prev = std::make_unique<D>(mut);
                if (now < was && is_repeatable(pos, mut))
                    scaling_state = scaling_state.succeed();
                else
                    scaling_state = scaling_state.fail();
            },

            // visit
            cbs.visit,
        };
    }
};

/*----------------------------------------------------------------------------*/
// structural_metropolis
//
// This is a structure_t-aware wrapper around metropolis that allows the user
//  to provide mutations via a Python script.

#ifdef SP2_ENABLE_PYTHON

// NOTE: This is the core logic of 'structural_metropolis', which
// has been pulled out into a value-returning function to preempt any future
// bugs related to early returns and/or accidental modification of sys.
//
// This function may leave 'sys' in an arbitrary (but consistent) state on exit.
template<typename S>
structure_t _structural_metropolis(
    S &sys,
    structure_t initial_pos,
    python::py_opaque_t extra_kw,
    structural_metropolis_settings_t met_set)
{
    using namespace sp2::python;

    initial_pos.reduce_positions();

    size_t natom = sys.get_structure().positions.size(); // FIXME

    auto diff_fn = [&](const auto &pos) {
        sys.set_structure(pos);
        sys.update();
        return std::make_pair(sys.get_value(), sys.get_gradient());
    };

    auto value_fn = [&diff_fn, &met_set](const auto &pos) {
        return compute_objective(met_set.objective, diff_fn(pos));
    };

    auto get_param_pack = [&diff_fn, &extra_kw](const auto &pos) {
        auto carts = v3tod(pos.positions);
        auto force = diff_fn(pos).second;
        for (auto & x : force)
            x *= -1.0;

#warning FIXME move that out of there
        auto pp = sp2::python::structural_metropolis::make_param_pack(carts, pos.lattice, force);
        return merge_dictionaries(pp, extra_kw, merge_strategy::ERROR);
    };

    // Fair warning: This struct is immediately constructed.
    struct py_mutation_functions_t
    {
        bool advanced = false;
        // class_ isn't here; it's only used to resolve the others.
        py_opaque_t mutate;
        py_opaque_t generate;
        py_opaque_t apply;
        py_opaque_t applied;
        py_opaque_t visit;
        py_opaque_t is_repeatable;
        py_opaque_t scale;

        // This constructor is where we do all of our interpretation of the
        // provided script and config, so that all of the above objects should
        // be resolved by the end.
        py_mutation_functions_t(
            const structural_metropolis_settings_t &met_set,
            const py_opaque_t &initial_pp)
        {
            extend_sys_path(met_set.python_sys_path);

            auto module = import(met_set.python_module);
            auto &func_set = met_set.python_functions;

            py_opaque_t instance = [&] {
                if (func_set.class_.empty())
                    return module;
                else
                {
                    // return an instance of the class instead
                    return module.getattr(func_set.class_)
                                 .call({}, initial_pp);
                }
            }();

            // The logic for resolving many of these pieces of
            // configuration is spread all over the place.
            //
            // All except 'generate' and 'mutate' follow the following rules:
            // - If specified in config, the function must exist.
            // - If not specified, a default name is assumed.
            // - If not specified and no function is present with this
            //    default name, a default implementation is used.
            //
            // Logic for the default implementation is currently provided
            // by the python::structural_metropolis::call_XXX family of
            // functions, and is triggered by the function being NULL.
            //
            // This null comes from the default supplied to getattr().
            auto lookup = [&](auto &member, const auto &fs_member, const char* default_name) {
                if (fs_member.empty())
                    member = instance.getattr(default_name, {});
                else
                    member = instance.getattr(fs_member);
            };

            if (func_set.advanced) {

                generate = instance.getattr(func_set.generate);

                lookup(apply, func_set.apply, "apply");
                lookup(applied, func_set.applied, "applied");
                lookup(visit, func_set.visit, "visit");
                lookup(is_repeatable, func_set.is_repeatable, "is_repeatable");
                lookup(scale, func_set.scale, "scale");
            }
            else
            {
                // basic interface
                mutate = instance.getattr(func_set.mutate);
            }
        }
    } functions(met_set, get_param_pack(initial_pos));

    //----------------------------------------------------------

    // These wrap the python methods, containing all the conversion code
    //  and providing the default definitions for missing callbacks.

    auto call_mutate = [&](auto &pos) {
        auto kw = get_param_pack(pos);
        return functions.mutate
                        .call({}, kw)
                        .template parse_as<structural_mutation_t>();
    };

    auto call_apply = [&](auto &pos, auto &mutation) {
        if (!functions.apply)
            return mutation.template parse_as<structural_mutation_t>();

        auto args = py_opaque_t::tuple(mutation);
        auto kw = get_param_pack(pos);
        return functions.apply
                        .call(args, kw)
                        .template parse_as<structural_mutation_t>();
    };

    auto call_applied = [&](auto &pos, auto &mutation, double was, double now) {
        if (!functions.applied)
            return;

        auto values = py_opaque_t::tuple(was, now);
        auto args = py_opaque_t::tuple(mutation, values);
        auto kw = get_param_pack(pos);
        functions.applied.call(args, kw);
    };

    auto call_visit = [&](auto &pos, double value, bool better) {
        if (!functions.visit)
            return;

        auto args = py_opaque_t::tuple(value, better);
        auto kw = get_param_pack(pos);
        functions.visit.call(args, kw);
    };

    auto call_generate = [&](auto &pos) {
        auto kw = get_param_pack(pos);
        return functions.generate.call({}, kw);
    };

    auto call_is_repeatable = [&](auto &pos, auto &mutation) {
        if (!functions.is_repeatable)
            return false;

        auto args = py_opaque_t::tuple(mutation);
        auto kw = get_param_pack(pos);
        return functions.is_repeatable
                        .call(args, kw)
                        .template parse_as<bool>();
    };

    auto call_scale = [&](auto &pos, auto &mutation, double factor) {
        if (!functions.scale)
            return mutation;

        auto args = py_opaque_t::tuple(mutation, factor);
        auto kw = get_param_pack(pos);
        return functions.scale.call(args, kw);
    };

    //----------------------------------------------------------

    // These are the callbacks for metropolis.
    // - 'basic_xyz_fn' is callback 'xyz' when "advanced = false"
    // - 'advanced_xyz_fn' is callback 'xyz' when "advanced = true"

    auto basic_generate_fn = call_mutate;
    auto basic_apply_fn = [&](auto &in_pos, auto &mutation) {
        auto pos = in_pos;

        switch (mutation.type)
        {
        case structural_mutation_type::LATTICE:
        {
            const as_ndarray_t<double> &arr = mutation.data;
            const std::vector<double> &d = arr.data();
            if (arr.shape() != std::vector<size_t>{3, 3})
                throw std::runtime_error("script produced wrong shape lattice");

            auto new_lattice = mat3x3_t{{d[0], d[1], d[2]},
                                        {d[3], d[4], d[5]},
                                        {d[6], d[7], d[8]}};
            pos.rescale(new_lattice);

            break;
        }

        case structural_mutation_type::CART_COORDS:
        {
            const as_ndarray_t<double> &arr = mutation.data;
            if (arr.shape() != std::vector<size_t>{natom, 3})
                throw std::runtime_error("script produced wrong shape carts");

            pos.positions = dtov3(arr.data());
            pos.reduce_positions();

            break;
        }

        case structural_mutation_type::FRAC_COORDS:
            // TODO maybe
            throw std::runtime_error(
                "mutation of frac coords not yet implemented");

        default:
            throw std::runtime_error("unreachable");
        }

        return pos;
    };
    auto basic_applied_fn = [&](...) {};
    auto basic_visit_fn = [&](...) {};

    auto advanced_generate_fn = call_generate;
    auto advanced_apply_fn = [&](auto &in_pos, auto &mutation) {
        // convert the script's own mutation type to structural_mutation_t...
        auto cxx_mutation = call_apply(in_pos, mutation);
        // ...and delegate to the
        return basic_apply_fn(in_pos, cxx_mutation);
    };
    auto advanced_applied_fn = call_applied;
    auto advanced_visit_fn = call_visit;
    auto advanced_is_repeatable_fn = call_is_repeatable;
    auto advanced_scale_fn = call_scale;

    if (met_set.python_functions.advanced)
    {
        using pos_t = structure_t;
        using diff_t = py_opaque_t;
        minimize::metropolis::callbacks_t<pos_t, diff_t> raw_callbacks{
            advanced_generate_fn,
            advanced_apply_fn,
            advanced_applied_fn,
            advanced_visit_fn,
        };

        auto scaling_control =
            minimize::metropolis::scaling_control_t<pos_t, diff_t>{
                met_set.scaling_settings,
                raw_callbacks,
                advanced_is_repeatable_fn,
                advanced_scale_fn,
            };
        auto callbacks = scaling_control.get_callbacks();

        return minimize::metropolis::advanced<pos_t, diff_t>
            (value_fn, callbacks, initial_pos, met_set.settings);
    }
    else
    {
        // Simple script interface.
        // Ironically, this still uses the advanced metropolis interface.
        // (perhaps we don't need metropolis::basic?)
        using pos_t = structure_t;
        using diff_t = structural_mutation_t;
        minimize::metropolis::callbacks_t<pos_t, diff_t> callbacks{
            basic_generate_fn,
            basic_apply_fn,
            [&](...) { },
        };
        return minimize::metropolis::advanced<pos_t, diff_t>
            (value_fn, callbacks, initial_pos, met_set.settings);
    }
};

/// Perform structure-aware metropolis minimization, leaving the optimal
/// structure in 'sys' on exit.
template<typename S>
void structural(S &sys,
    python::py_opaque_t extra_kw,
    structural_metropolis_settings_t met_set)
{
    auto initial = sys.get_structure();
    auto final = _structural_metropolis(sys, initial, extra_kw, met_set);
    sys.set_structure(final);
    sys.update();
}

#endif // SP2_ENABLE_PYTHON

} // namespace metropolis
} // namespace minimize
} // namespace sp2
#endif // SP2_METROPOLIS_HPP
