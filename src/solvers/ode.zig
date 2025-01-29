const SolverError = error{
    InvalidTimeDelta,
};

pub const IntegrationMethod = enum {
    euler,
    huen,
};
pub fn ODE(comptime N: usize) type {
    const DerivativeFn = *const fn (x: [N]f64, t: f64) [N]f64;

    return struct {
        x: [N]f64,
        t: f64,
        derivative: DerivativeFn,

        pub fn init(
            x: [N]f64,
            t: f64,
            derivative: DerivativeFn,
        ) @This() {
            return .{
                .x = x,
                .t = t,
                .derivative = derivative,
            };
        }

        pub fn step(
            self: @This(),
            dt: f64,
            method: IntegrationMethod,
        ) !@This() {
            if (dt <= 0) return SolverError.InvalidTimeDelta;

            switch (method) {
                .euler => return Euler(N, self, dt),
                .huen => return Huen(N, self, dt),
            }
        }
    };
}

pub fn Euler(comptime N: usize, ode: ODE(N), dt: f64) !ODE(N) {
    const dx = ode.derivative(ode.x, ode.t);

    const Vec = @Vector(N, f64);
    const dx_vec = @as(Vec, dx);
    const x_vec = @as(Vec, ode.x);
    const dt_vec: Vec = @splat(dt);

    const x_next_vec = x_vec + dt_vec * dx_vec;
    const x_next = @as([N]f64, x_next_vec);

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

    const k1 = ode.derivative(ode.x, ode.t);
    const k1_vec = @as(Vec, k1);
    const x_vec = @as(Vec, ode.x);

    const x_intermediate_vec = x_vec + dt_vec * k1_vec;
    const x_intermediate = @as([N]f64, x_intermediate_vec);

    const k2 = ode.derivative(x_intermediate, ode.t + dt);
    const k2_vec = @as(Vec, k2);

    const x_next_vec = x_vec + half_dt_vec * (k1_vec + k2_vec);
    const x_next = @as([N]f64, x_next_vec);

    return .{
        .x = x_next,
        .t = ode.t + dt,
        .derivative = ode.derivative,
    };
}
