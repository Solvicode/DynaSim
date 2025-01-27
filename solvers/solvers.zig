pub fn ODE(N: usize) type {
    return struct {
        x: @Vector(N, f64),
        t: f64,
        derivative: *const fn (x: @Vector(N, f64), t: f64) @Vector(N, f64),

        pub fn init(
            x: @Vector(N, f64),
            t: f64,
            derivative: *const fn (x: @Vector(N, f64), t: f64) @Vector(N, f64),
        ) @This() {
            return .{
                .x = x,
                .t = t,
                .derivative = derivative,
            };
        }
    };
}

pub fn EulerStep(ode: ODE, dt: f64, t: f64) ODE {
    const x_next = ode.x + dt * ode.derivative(ode.x, t);
    return .{
        .x = x_next,
        .t = t,
        .derivative = ode.derivative,
    };
}

pub fn HuenStep(ode: ODE, dt: f64, t: f64) ode {
    const k1 = ode.derivative(ode.x, t);
    const k2 = ode.derivative(ode.x + dt * k1, t + dt);
    const x_next = ode.x + (dt / 2) * (k1 + k2);

    return .{
        .x = x_next,
        .t = t,
        .derivative = ode.derivative,
    };
}
