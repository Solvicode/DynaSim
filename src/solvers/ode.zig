pub const SolverError = error{
    InvalidTimeDelta,
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
        t: f64,
        derivative: DerivativeFn,

        pub fn init(
            x: @Vector(StateDimension, f64),
            u: ?@Vector(InputDimension, f64),
            t: f64,
            derivative: DerivativeFn,
        ) @This() {
            return .{
                .x = x,
                .u = u,
                .t = t,
                .derivative = derivative,
            };
        }

        pub fn step(
            self: @This(),
            u: ?@Vector(InputDimension, f64),
            dt: f64,
            method: IntegrationMethod,
            adaptive: ?bool,
            tolerance: ?f64,
        ) !@This() {
            if (dt <= 0) return SolverError.InvalidTimeDelta;

            // return a new object
            const new_state = @This(){
                .x = self.x,
                .u = u,
                .t = self.t,
                .derivative = self.derivative,
            };

            switch (method) {
                .euler => return Euler(StateDimension, InputDimension, new_state, dt),
                .huen => return Huen(StateDimension, InputDimension, new_state, dt),
                .rk4 => return RungeKutta4(
                    StateDimension,
                    InputDimension,
                    new_state,
                    dt,
                    adaptive,
                    tolerance,
                ),
            }
        }
    };
}

pub fn Euler(
    comptime StateDimension: usize,
    comptime InputDimension: usize,
    ode: ODE(StateDimension, InputDimension),
    dt: f64,
) !ODE(StateDimension, InputDimension) {
    const dt_vec: @Vector(StateDimension, f64) = @splat(dt);
    const dx = ode.derivative(ode.x, ode.u, ode.t);
    const x_next = ode.x + dt_vec * dx;

    return .{
        .x = x_next,
        .u = ode.u,
        .t = ode.t + dt,
        .derivative = ode.derivative,
    };
}

pub fn Huen(
    comptime StateDimension: usize,
    comptime InputDimension: usize,
    ode: ODE(StateDimension, InputDimension),
    dt: f64,
) !ODE(StateDimension, InputDimension) {
    const Vec = @Vector(StateDimension, f64);
    const dt_vec: Vec = @splat(dt);
    const half_dt_vec: Vec = @splat(dt / 2.0);

    const k1 = ode.derivative(ode.x, ode.u, ode.t);
    const x_intermediate = ode.x + dt_vec * k1;
    const t_intermediate = ode.t + dt;
    const k2 = ode.derivative(x_intermediate, ode.u, t_intermediate);
    const x_next = ode.x + half_dt_vec * (k1 + k2);

    return .{
        .x = x_next,
        .u = ode.u,
        .t = ode.t + dt,
        .derivative = ode.derivative,
    };
}

pub fn RungeKutta4(
    comptime StateDimension: usize,
    comptime InputDimension: usize,
    ode: ODE(StateDimension, InputDimension),
    dt: f64,
    adaptive: ?bool,
    tolerance: ?f64,
) !ODE(StateDimension, InputDimension) {
    const Vec = @Vector(StateDimension, f64);
    var current_dt = dt;
    var current_ode = ode;
    const tol = tolerance orelse 1e-6;

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
            // perform calculation again with half step size
            const half_dt = current_dt / 2.0;
            const half_step1 = try RungeKutta4(StateDimension, InputDimension, current_ode, half_dt, false, null);
            const half_step2 = try RungeKutta4(StateDimension, InputDimension, half_step1, half_dt, false, null);

            // estimate the error
            const error_vec = @abs(x_next - half_step2.x);
            const max_error = @reduce(.Max, error_vec);

            if (max_error > tol) {
                // if error is too large decrease the step size
                current_dt = current_dt / 2.0;
                continue;
            } else if (max_error < tol / 10.0) {
                // if error is very small increase the step size
                current_dt = current_dt * 2.0;
                continue;
            }
        }

        return .{
            .x = x_next,
            .u = current_ode.u,
            .t = current_ode.t + current_dt,
            .derivative = current_ode.derivative,
        };
    }
}
