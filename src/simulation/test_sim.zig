const std = @import("std");
const core_sim = @import("core_sim.zig");
const ODE = @import("../simulation/ode.zig").ODE;
const IntegrationMethod = @import("../simulation/ode.zig").IntegrationMethod;

const StateDimension = 2;
const InputDimension = 1;

fn simple_derivative(x: @Vector(StateDimension, f64)) @Vector(StateDimension, f64) {
    // Exampe derivative: dx/dt = -x
    var result: @Vector(StateDimension, f64) = undefined;
    result[0] = -x[0];
    result[1] = -x[1];
    return result;
}

pub fn test_simulation() void {
    const initial_dt = 0.1;
    const end_time = 2.0;

    // using Euler method
    std.debug.print("Running simulation with Euler method...\n", .{});
    _ = core_sim.simulate(StateDimension, InputDimension, initial_dt, end_time, simple_derivative, IntegrationMethod.euler);

    // using Runge-Kutta method
    std.debug.print("Running simulation with Runge-Kutta method...\n", .{});
    _ = core_sim.simulate(StateDimension, InputDimension, initial_dt, end_time, simple_derivative, IntegrationMethod.rk4);
}
