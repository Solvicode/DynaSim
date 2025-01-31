pub const SolverError = error{
    InvalidTimeDelta,
};

pub const IntegrationMethod = enum {
    euler,
    huen,
};

pub fn ODE(comptime StateDimension: usize, comptime InputDimension: usize) type {
    const DerivativeFn = *const fn (x: @Vector(StateDimension, f64), u: ?@Vector(InputDimension, f64), t: ?f64) @Vector(StateDimension, f64);

    return struct {
        x: @Vector(StateDimension, f64),
        u: ?@Vector(InputDimension, f64),
        t: ?f64,
        derivative: DerivativeFn,

        pub fn init(
            x: @Vector(StateDimension, f64),
            u: ?@Vector(InputDimension, f64),
            t: ?f64,
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
        ) !@This() {
            if (dt <= 0) return SolverError.InvalidTimeDelta;

            // assign the input
            self.u = u;

            switch (method) {
                .euler => return Euler(StateDimension, self, dt),
                .huen => return Huen(StateDimension, self, dt),
            }
        }
    };
}

pub fn Euler(comptime N: usize, ode: ODE(N), dt: f64) !ODE(N) {
    const dt_vec: @Vector(N, f64) = @splat(dt);
    const dx = ode.derivative(ode.x, ode.u, ode.t);
    const x_next = ode.x + dt_vec * dx;

    return .{
        .x = x_next,
        .t = ode.t + dt,
        .derivative = ode.derivative,
    };
}

pub fn Huen(comptime N: usize, ode: ODE(N), dt: f64) !ODE(N) {
    const Vec = @Vector(N, f64);
    const dt_vec: Vec = @splat(dt);
    const half_dt_vec: Vec = @splat(dt / 2.0);

    const k1 = ode.derivative(ode.x, ode.u, ode.t);

    const x_intermediate = ode.x + dt_vec * k1;

    const k2 = ode.derivative(x_intermediate, ode.u, ode.t + dt);

    const x_next = ode.x + half_dt_vec * (k1 + k2);

    return .{
        .x = x_next,
        .t = ode.t + dt,
        .derivative = ode.derivative,
    };
}
