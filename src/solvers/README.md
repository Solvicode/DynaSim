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
x_{k+1}= x_{k} + \Delta t \cdot f(x_k, t_k)
$$

The function, $f$ must be able to be evaluated at time $t_k$

## Huen's Method

Heun's method is an extension to Euler's, where Euler's is treated as an intermediate step. Instead, however,
a second step is performed where the current state is updated with the step-size multipled average of the
derivative function evaluation at $t_k$ and derivative function evaluation at $t_{k+1}$, but with the
estimated state given by Euler's method.

$$
\tilde{x}_{k+1}= x_k + \Delta t \cdot f(x_k, t_k) \\
x_{k+1} = x_k + \frac{\Delta t}{2} \cdot (f(x_k, t_k) + f(\tilde{x}_{k+1}, t_{k+1}))
$$

When expanded it's clear that Heun's method is equivalent to a 2 stage second order Runge-Kutta.

$$
\begin{align*}
\text{Let;} \\
& & k_1     &= f(x_k, t_k) \\
& & k_2     &= f(x_k + \Delta t \cdot k_1, t_{k+1}) \\
\text{Then;} \\
& & x_{k+1} &= x_k + \frac{\Delta t}{2}(k_1 + k_2)
\end{align*}
$$

## Runge-Kutta

The 4th order Runge-Kutta method flows well off of Huen's method, and the pattern starts to become clear.

$$
\begin{align*}
\text{Let;} & & & \\
& & k_1 &= f(t_k, x_k)  \\
& & k_2 &= f(t_k + \frac{\Delta t}{2} ,  \frac{\Delta t \cdot k_1}{2})  \\
& & k_3 &= f(t_k + \frac{\Delta t}{2} ,  \frac{\Delta t \cdot k_2}{2})  \\
& & k_4 &= f(t_k + \Delta t ,  \Delta t \cdot k_3 )  \\
\text{Then;} \\
& & x_{k+1} &= x_k + \frac{1}{6} \left( k_1 + 2k_2 + 2k_3 + k4 \right)
\end{align*}
$$

Runge-Kutta 4 is considered a best in class solver, and is generally the first attempted solver for complex
ODE problems. Extensions of this method exist including explicit methods where the method is generalised into
the following form:

$$
\begin{align*}
& & x_{k+1} &= x_k + \Delta t \sum_{i=1}^n b_i k_i \\
& & k_1 &= f(t_k, x_k) \\
& & k_2 &= f(t_k + c_2 \Delta t, x_k + (a_{21}k_1) \Delta t) \\
& & k_3 &= f(t_k + c_3 \Delta t, x_k + (a_{31}k_1 + a_{32}k_2) \Delta t) \\
& & k_s &= f(t_k + c_s \Delta t, x_k + (a_{s1}k_1 + \dots + a_{s,s-1}k_{s-1})\Delta t) \\
\end{align*}
$$

The generalised coefficients, $a_n, b_i, c_s$ are 'free' and can be tuned in any manner of ways to form
a specific solving regime. They can be well summarised into the following form:

$$
\begin{array}{c|cccc}
0 & & & & \\
c_2 & a_{21} & & & \\
c_3 & a_{31} & a_{32} & & \\
\vdots & \vdots & \vdots & \ddots & \\
c_s & a_{s1} & a_{s2} & \cdots & a_{s,s-1} \\
\hline
& b_1 & b_2 & \cdots & b_s
\end{array}
$$
