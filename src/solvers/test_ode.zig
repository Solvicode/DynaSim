// test_ode.zig
const std = @import("std");
const ode = @import("ode.zig");

// static vector field that always returns [0.0, 0.0]
fn staticField(x: [2]f64, t: f64) [2]f64 {
    _ = x;
    _ = t;
    return .{ 0.0, 0.0 };
}

test "static field preserves state" {

    // test setup
    const solvers: [2]ode.IntegrationMethod = .{ .euler, .huen };

    // initial conditions
    const x0 = [2]f64{ 0.0, 0.0 };
    const t0 = 0.0;
    const dt = 0.1;

    // create initial ODE state
    var initial_ode = ode.ODE(2).init(x0, t0, staticField);
    const initial_x = initial_ode.x;

    inline for (solvers) |solver| {
        // take a step
        const next_ode = try initial_ode.step(dt, solver);

        // check that original state is unchanged
        try std.testing.expectEqual(initial_ode.x, initial_x);

        // for a static field [0,0], after time dt the position should be [0,0]
        try std.testing.expectEqual(next_ode.x, initial_x);

        // test multiple steps to ensure consistency
        const next_next_ode = try next_ode.step(dt, solver);
        try std.testing.expectEqual(next_next_ode.x, initial_x);
    }
}

test "zero time delta gives error" {}
