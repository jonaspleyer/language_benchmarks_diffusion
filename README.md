# Language benchmark to solve the diffusion equation
## Mathmatical Setting
We benchmark different languages with reasonable optimizations enabled (such as -O3, target=native, ...) in their ability to numerically solve the 2D diffusion equation 
$$
\partial_t u = F(u,v) + D_u\Delta u\\
\partial_t v = G(u,v) + D_v\Delta v
$$
with source terms given by a schnakenberg substrate depletion model.
$$\begin{aligned}
F(u,v) &= &&k_1 - k_2u &+ k_3u^2v\\
G(u,v) &= &&k_2v &- k_3u^2v
\end{aligned}$$
To model this we use a finite difference approach
$$\begin{aligned}
\frac{u^{n+1}_{i,j} - u^n_{i,j}}{\Delta t} &= &&F(u^n_{i,j},v^n_{i,j}) &&+ D_u\Big(\frac{u^n_{i+1,j} - 2u^n_{i,j} + u^n_{i-1,j}}{\Delta x^2} + \frac{u^n_{i,j+1} - 2u^n_{i,j} + u^n_{i,j-1}}{\Delta y^2}\Big)\\
\frac{v^{n+1}_{i,j} - v^n_{i,j}}{\Delta t} &= &&G(u^n_{i,j},v^n_{i,j}) &&+ D_v\Big(\frac{v^n_{i+1,j} - 2v^n_{i,j} + v^n_{i-1,j}}{\Delta x^2} + \frac{v^n_{i,j+1} - 2v^n_{i,j} + v^n_{i,j-1}}{\Delta y^2}\Big)
\end{aligned}$$
such that the iterative calculation we want to perform is:
$$\begin{aligned}
u^{n+1}_{i,j} &= &\Delta tF(u^n_{i,j},v^n_{i,j}) &+ \frac{\Delta t D_u}{\Delta x^2}(u^n_{i+1,j} - 2u^n_{i,j} + u^n_{i-1,j}) &+ \frac{\Delta t D_u}{\Delta y^2}(u^n_{i,j+1} - 2u^n_{i,j} + u^n_{i,j-1}) &+ u^n_{i,j}\\
v^{n+1}_{i,j} &= &\Delta tG(u^n_{i,j},v^n_{i,j}) &+ \frac{\Delta t D_v}{\Delta x^2}(v^n_{i+1,j} - 2v^n_{i,j} + v^n_{i-1,j}) &+ \frac{\Delta t D_v}{\Delta x^2}(v^n_{i,j+1} - 2v^n_{i,j} + v^n_{i,j-1}) &+ v^n_{i,j}
% TODO make this formula correct
\end{aligned}$$
We solve this equation with Neumann-Boundary conditions $\partial_x u,\partial_x v=0$ and initialise with a random white noise.
A complete list of parameters can be found ...

## Languages
Current
Future Plans
- Bash
- C
- C++
- Fortran
- Haskell
- Julia
- LLVM assembly
- Python
- Rust
- x86-64 assembly

We also plan to provide a PKGBUILD file to the Arch-User-Repository
## Overview Results