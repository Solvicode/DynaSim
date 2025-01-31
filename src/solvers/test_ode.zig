// test_ode.zig
const std = @import("std");
const ode = @import("ode.zig");

const solvers: [3]ode.IntegrationMethod = .{
    .euler,
    .huen,
    .rk4,
};

// static vector field that always returns [0.0, 0.0]
fn staticField(x: @Vector(2, f64), u: ?@Vector(0.0, f64), t: ?f64) @Vector(2, f64) {
    _ = x;
    _ = u;
    _ = t;
    return .{ 0.0, 0.0 };
}

fn secondOrderSystem(x: @Vector(2, f64), u: ?@Vector(2, f64), t: ?f64) @Vector(2, f64) {
    _ = t; // no time dependence
    const input = u orelse @Vector(2, f64){ 0, 0 };
    return .{
        x[1],
        -x[0] - 2 * x[1] + input[1],
    };
}

test "static field preserves state" {

    // initial conditions
    const xinit = @Vector(2, f64){ 0.0, 0.0 };
    const t0 = 0.0;
    const dt = 0.1;

    // create initial ODE state
    var initial_ode = ode.ODE(2, 0).init(xinit, null, t0, staticField);
    const initial_x = initial_ode.x;

    inline for (solvers) |solver| {
        // take a step
        const next_ode = try initial_ode.step(null, dt, solver, null, null);

        // check that original state is unchanged
        try std.testing.expectEqual(initial_ode.x, initial_x);

        // for a static field [0,0], after time dt the position should be [0,0]
        try std.testing.expectEqual(next_ode.x, initial_x);

        // test multiple steps to ensure consistency
        const next_next_ode = try next_ode.step(null, dt, solver, null, null);
        try std.testing.expectEqual(next_next_ode.x, initial_x);
    }
}

test "zero time delta gives error" {
    const xinit = @Vector(2, f64){ 0.0, 0.0 };
    const t0 = 0.0;
    const dt = 0.0;

    var initial_ode = ode.ODE(2, 0).init(xinit, null, t0, staticField);
    inline for (solvers) |solver| {
        const result = initial_ode.step(null, dt, solver, null, null);
        try std.testing.expectError(ode.SolverError.InvalidTimeDelta, result);
    }
}

test "final value is as expected" {
    // initial conditions
    const xinit = @Vector(2, f64){ 0.0, 0.0 };
    const ustep = @Vector(2, f64){ 0.0, 1.0 };
    const t0 = 0.0;
    const dt = 0.001;
    const initial_ode = ode.ODE(2, 2).init(xinit, ustep, t0, secondOrderSystem);
    inline for (solvers) |solver| {
        var current_ode = initial_ode;
        var t: f64 = t0;

        var i: usize = 0;
        while (i < 10000) : (i += 1) {
            current_ode = try current_ode.step(ustep, dt, solver, null, null);
            t += dt;
        }
        try std.testing.expectApproxEqAbs(1.0, current_ode.x[0], 1e-3);
        try std.testing.expectApproxEqAbs(0.0, current_ode.x[1], 1e-3);
    }
}
