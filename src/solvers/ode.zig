pub const SolverError = error{
    InvalidTimeDelta,
    ShrunkTimeDelta,
};

pub const IntegrationMethod = enum {
    euler,
    huen,
    rk4,
};

pub fn ODE(comptime StateDimension: usize, comptime InputDimension: usize) type {
    const DerivativeFn = *const fn (
        x: @Vector(StateDimension, f64),
        u: ?@Vector(InputDimension, f64),
        t: ?f64,
    ) @Vector(StateDimension, f64);

    return struct {
        x: @Vector(StateDimension, f64),
        u: ?@Vector(InputDimension, f64),
        dt: f64, // timestep is stored in the ODE func because solvers will modify it
        t: f64,
        derivative: DerivativeFn,

        pub fn init(
            x: @Vector(StateDimension, f64),
            u: ?@Vector(InputDimension, f64),
            dt: f64,
            t: f64,
            derivative: DerivativeFn,
        ) @This() {
            return .{
                .x = x,
                .u = u,
                .dt = dt,
                .t = t,
                .derivative = derivative,
            };
        }

        pub fn step(
            self: @This(),
            u: ?@Vector(InputDimension, f64),
            method: IntegrationMethod,
            adaptive: ?bool,
            stateStolerance: ?f64,
            timeTolerance: ?f64,
        ) !@This() {
            if (self.dt <= 0) {
                return error.InvalidTimeDelta;
            }
            // return a new object
            const new_state = @This(){
                .x = self.x,
                .u = u,
                .dt = self.dt,
                .t = self.t,
                .derivative = self.derivative,
            };

            switch (method) {
                .euler => return Euler(StateDimension, InputDimension, new_state),
                .huen => return Huen(StateDimension, InputDimension, new_state),
                .rk4 => return RungeKutta4(
                    StateDimension,
                    InputDimension,
                    new_state,
                    adaptive,
                    stateStolerance,
                    timeTolerance,
                ),
            }
        }
    };
}

pub fn Euler(
    comptime StateDimension: usize,
    comptime InputDimension: usize,
    ode: ODE(StateDimension, InputDimension),
) !ODE(StateDimension, InputDimension) {
    const dt_vec: @Vector(StateDimension, f64) = @splat(ode.dt);
    const dx = ode.derivative(ode.x, ode.u, ode.t);
    const x_next = ode.x + dt_vec * dx;

    return .{
        .x = x_next,
        .u = ode.u,
        .dt = ode.dt,
        .t = ode.t + ode.dt,
        .derivative = ode.derivative,
    };
}

pub fn Huen(
    comptime StateDimension: usize,
    comptime InputDimension: usize,
    ode: ODE(StateDimension, InputDimension),
) !ODE(StateDimension, InputDimension) {
    const Vec = @Vector(StateDimension, f64);
    const dt_vec: Vec = @splat(ode.dt);
    const half_dt_vec: Vec = @splat(ode.dt / 2.0);

    const k1 = ode.derivative(ode.x, ode.u, ode.t);
    const x_intermediate = ode.x + dt_vec * k1;
    const t_intermediate = ode.t + ode.dt;
    const k2 = ode.derivative(x_intermediate, ode.u, t_intermediate);
    const x_next = ode.x + half_dt_vec * (k1 + k2);

    return .{
        .x = x_next,
        .u = ode.u,
        .dt = ode.dt,
        .t = ode.t + ode.dt,
        .derivative = ode.derivative,
    };
}

pub fn RungeKutta4(
    comptime StateDimension: usize,
    comptime InputDimension: usize,
    ode: ODE(StateDimension, InputDimension),
    adaptive: ?bool,
    stateTolerance: ?f64,
    timeTolerance: ?f64,
) !ODE(StateDimension, InputDimension) {
    const Vec = @Vector(StateDimension, f64);
    var current_dt = ode.dt;
    var current_ode = ode;
    const stateTol = stateTolerance orelse 1e-6;
    const timeTol = timeTolerance orelse 1e-8;

    while (true) {
        const dt_vec: Vec = @splat(current_dt);
        const dt_half: f64 = current_dt / 2.0;
        const dt_sixth: Vec = @splat(current_dt / 6.0);

        const half: @Vector(StateDimension, f64) = @splat(0.5);
        const double: @Vector(StateDimension, f64) = @splat(2.0);

        const k1 = current_ode.derivative(current_ode.x, current_ode.u, current_ode.t);

        const x2 = current_ode.x + dt_vec * half * k1;
        const k2 = current_ode.derivative(x2, current_ode.u, current_ode.t + dt_half);

        const x3 = current_ode.x + dt_vec * half * k2;
        const k3 = current_ode.derivative(x3, current_ode.u, current_ode.t + dt_half);

        const x4 = current_ode.x + dt_vec * k3;
        const k4 = current_ode.derivative(x4, current_ode.u, current_ode.t + current_dt);

        const x_next = current_ode.x + dt_sixth * (k1 + double * k2 + double * k3 + k4);

        if (adaptive orelse false) {
            // calculate with half step size
            const half_dt_vec: @Vector(StateDimension, f64) = @splat(dt_half);
            const half_dt_sixth: @Vector(StateDimension, f64) = @splat(dt_half / 6.0);

            // first half step
            const k1_h1 = current_ode.derivative(current_ode.x, current_ode.u, current_ode.t);
            const x2_h1 = current_ode.x + half_dt_vec * half * k1_h1;
            const k2_h1 = current_ode.derivative(x2_h1, current_ode.u, current_ode.t + dt_half / 2.0);
            const x3_h1 = current_ode.x + half_dt_vec * half * k2_h1;
            const k3_h1 = current_ode.derivative(x3_h1, current_ode.u, current_ode.t + dt_half / 2.0);
            const x4_h1 = current_ode.x + half_dt_vec * k3_h1;
            const k4_h1 = current_ode.derivative(x4_h1, current_ode.u, current_ode.t + dt_half);
            const x_mid = current_ode.x + half_dt_sixth * (k1_h1 + double * k2_h1 + double * k3_h1 + k4_h1);

            // second half step
            const k1_h2 = current_ode.derivative(x_mid, current_ode.u, current_ode.t + dt_half);
            const x2_h2 = x_mid + half_dt_vec * half * k1_h2;
            const k2_h2 = current_ode.derivative(x2_h2, current_ode.u, current_ode.t + dt_half * 1.5);
            const x3_h2 = x_mid + half_dt_vec * half * k2_h2;
            const k3_h2 = current_ode.derivative(x3_h2, current_ode.u, current_ode.t + dt_half * 1.5);
            const x4_h2 = x_mid + half_dt_vec * k3_h2;
            const k4_h2 = current_ode.derivative(x4_h2, current_ode.u, current_ode.t + current_dt);
            const x_next_half = x_mid + half_dt_sixth * (k1_h2 + double * k2_h2 + double * k3_h2 + k4_h2);

            // estimate the error
            const error_vec = @abs(x_next - x_next_half);
            const max_error = @reduce(.Max, error_vec);

            if (max_error > stateTol) {
                // if error is too large decrease the step size
                current_dt = current_dt / 2.0;
                if (current_dt < timeTol) {
                    return error.ShrunkTimeDelta;
                }
                continue;
            } else if (max_error < stateTol / 10.0) {
                // if error is very small increase the step size
                current_dt = current_dt * 2.0;
                continue;
            }
        }

        return .{
            .x = x_next,
            .u = current_ode.u,
            .dt = current_dt,
            .t = current_ode.t + current_dt,
            .derivative = current_ode.derivative,
        };
    }
}
