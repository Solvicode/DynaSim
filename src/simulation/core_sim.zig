const std = @import("std");
const ODE = @import("../solvers/ode.zig.").ODE;
const IntegrationMethod = @import("../solvers/ode.zig.").IntegrationMethod;
const SolverError = @import("../solvers/ode.zig.").SolverError;

pub fn simulate(StateDimension: usize, InputDimension: usize, initial_dt: f64, end_time: f64, derivative: fn (@Vector(StateDimension, f64)) @Vector(StateDimension, f64), method: IntegrationMethod) !void {
    // Create initial state (assuming 0s for now)
    var state = ODE(StateDimension, InputDimension).init(
        @Vector([]f64, StateDimension),
        null,
        initial_dt,
        0.0,
        derivative,
    );

    var current_time = 0.0;
    while (current_time < end_time) {
        state = try state.step(null, method, null, null, null);

        std.debug.print("Time: {}, State: {}\n", .{ state.t, state.x });

        current_time += state.dt;
    }
}
