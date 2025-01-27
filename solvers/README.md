[[toc]]

With all these integration methods, the following is true:

$$
\frac{dx}{dt} = f(x, t) \\
t_{k+1} = t_k + \Delta t
$$

## Eulers Method

Eulers method performs a step update based on the raw evaluation of the function, over a finite step size.

As the step size approaches zero, Eulers method approaches the true value as the step approaches zero.

### Method

$$
x_{k+1}= x_{k-1} + \Delta t \cdot f(x_k, t_k)
$$

The function, $f$ must be able to be evaluated at time $t_k$

## Huen's Method

Heun's method is an extension to Euler's, where Euler's is treated as an intermediate step. Instead, however,
a second step is performed where the current state is updated with the step-size multipled average of the
derivative function evaluation at $t_k$ and derivative function evaluation at $t_{k+1}$, but with the
estimated state given by Euler's method.

$$
\tilde{x}_{k+1}= x_{k-1} + \Delta t \cdot f(x_k, t_k) \\
\tilde{x}_{k+1} = x_k + \frac{\Delta t}{2} \cdot (f(t_k, x_k) + f(t_{k+1}, \tilde{x}_{k+1}))
$$
