# Quantum Mechanics

## The Schrodinger Equation
Quantum theory essentially represents the theory of energy and matter on the smallest scales of nature.  It has foundations in the theory of 
waves, however the ideas concerned can be used to develop all the modern theories of Chemistry, Physics and Structural Biology.

### Basic Wave Theory}
The wave equation, in 1D, in classical physics is given by:
```{math}
\frac{\partial^2 \Psi}{\partial t^2} = c^2 \frac{\partial^2 \Psi}{\partial x^2}
```
where $\Psi$ represents some displacement from equilibrium. If we want to extend this to 2D or 3D we can use the Laplacian 
$\nabla^2 = \nabla\cdot \nabla$:
```{math}
\frac{\partial^2 \Psi}{\partial t^2} = c^2 \nabla^2 \Psi
```
Finding solutions to PDEs is typically tricky, however through analysis of the PDE, we find they have the form:
```{math}
\Psi(x,\,t) = A\,\exp(x - ct) + B\,\exp(x + ct)
```
where $A,\,B$ are constants.  We can show that the $A$ solutions represent travelling waves moving left to right (forwards) and the 
$B$ solutions represent travelling waves moving right to left (backwards).  We can find solution of this PDE that are of a separable form:
```{math}
\Psi(x,\,t) = f(x)\,g(t)
```
substituting this in gives us:
```{math}
f\,\frac{\mathrm{d}^2 g}{\mathrm{d} t^2} = c^2 g\,\frac{\mathrm{d}  f}{\mathrm{d}^2 x^2}
```
where we find the partial derivatives are are not just regular one variable derivatives.  Rearranging gives:
```{math}
\frac{1}{c^2\,g}\,\frac{\mathrm{d}^2 g}{\mathrm{d} t^2} = \frac{1}{f}\,\frac{\mathrm{d}  f}{\mathrm{d}^2 x^2}
```
Looking at this equation, since we separated out the $x$ and $t$ dependence's into functions $f$ and $g$ respectively, then although the 
LHS here appears to <b>only</b> depend on $t$, since it equals the RHS and that <b>does not</b> depend on $t$, but the LHS cannot depend 
on $t$ either.  Likewise looking at the RHS, it appears to <b>only</b> depend on $x$, but since it equals the LHS and that <b>does not</b> 
depend on $x$, the RHS cannot depend on $x$ either.  Therefore the only thing both sides can depend on is a constant, which we will 
denote as $-k^2$:
```{math}
\frac{1}{c^2\,g}\,\frac{\mathrm{d}^2 g}{\mathrm{d} t^2} = -k^2 &\Longrightarrow& \frac{\mathrm{d}^2 g}{\mathrm{d} t^2} = -k^2\,c^2\,g\\ 
\frac{1}{f}\,\frac{\mathrm{d}^2 f}{\mathrm{d} x^2} = -k^2 &\Longrightarrow& \frac{\mathrm{d}^2 f}{\mathrm{d} x^2} = -k^2\,f
```
We can solve these resulting ODEs:
```{math}
f(x) &=& C\,\exp(i\,k\,x) + D\,\exp(-i\,k\,x) \\
g(t) &=& E\,\exp(i\,k\,c\,t) + F\,\exp(-i\,k\,c\,t) 
```
where $C,\,D,\,E,\,F$ are constants.  If $x,\,t$ are variables that carry units $[\text{length}]$, $[\text{times}]$ respectively, then the arguments 
of the functions should be dimensionless, so lets call $k$ our wavevector and $kc = \omega$ our angular frequency.  We notice that these 
solutions are periodic:
```{math}
\omega t \rightarrow \omega t + 2\pi &\Longrightarrow& t \rightarrow  t + \frac{2\pi}{\omega} = t + T\\
kx \rightarrow kx + 2\pi &\Longrightarrow& x \rightarrow x + \frac{2\pi}{k} = x + \lambda
```
therefore we find the regular wave quantities $T = \frac{2\pi}{\omega}$ and $\lambda = \frac{2\pi}{k}$.  Lets take a trial wave solution:
```{math}
\Psi = A \,\exp(i(kx - \omega t))
```
and substitute into the wave equation:
```{math}
\frac{\partial^2 \Psi}{\partial t^2} = c^2 \frac{\partial^2 \Psi}{\partial x^2} \Rightarrow -\omega^2\,\Psi = -c^2 \,k^2\,\Psi
```
and therefore 
```{math}
\omega = \pm k\,c
```
where each sign here represents the left or right moving waves we discussed before.  We call such a formula that links together 
the angular frequency with wavevector a <b>dispersion relation</b>.  In higher dimensions, we find that we have wave solutions of the form:
```{math}
\Psi({\bf x},\,t) = A\,\exp(i({\bf k}\cdot {\bf x} - \omega t))
```
and the dispersion relation depends on $k = |{\bf k}|$.  So far we have used the wave equation to construct the dispersion relation for 
travelling waves, but we can also try to reverse engineering this process to construct the wave equation from the form of the dispersion relation. 

### Quantum Matter Waves
Lets start with some basic quantum physics, for a particle moving with momentum ${\bf p}$, its magnitude $p = |{\bf p}|$ is related to the 
de Broglie wavelength $\lambda$ by:
```{math}
\lambda = \frac{h}{p} 
```
which we can rewrite in terms of wavevector $k = |{\bf k}| = \frac{2\pi}{\lambda}$
```{math}
:label: debrogliep
p = \hbar k
```
where $\hbar = \frac{h}{2\pi}$ is the reduced Planck's constant.
Likewise the Planck equation for energy $E$ and frequency $f$, can be rewritten using the angular frequency $\omega = 2 \pi f$:
```{math}
:label: planckE
E = hf = \hbar \omega
```
For a free particle, moving with velocity ${\bf v}$, the only contribution to the energy $E$ is the kinetic energy 
$\frac{1}{2}mv^2$ where $v^2 = {\bf v}\cdot{\bf v}$ and since ${\bf p} = m{\bf v}$, this means we can rewrite 
```{math}
E = \frac{p^2}{2m}
```
where $p^2 = {\bf p}\cdot{\bf p}$.  Substituting in the quantum equations {eq}`debrogliep`-{eq}`planckE`Sxx	 here gives:
```{math}
\hbar \omega = \frac{\hbar^2 k^2}{2m}
```
producing a radically different dispersion relation to those from travelling waves.  But can we construct a wave equation that 
similarly produces this? Well if we have the complex form of wave as
```{math}
\psi(x,t) = \psi_0 \exp(i(kx-\omega t)) 
```
then we can see that in order to get a $k^2$ here we need a $\frac{\partial^2}{\partial x^2}$ term in the wave equation and for a 
$\omega$ term we need a $\frac{\partial}{\partial t}$ term in the wave equation.  

\noindent These will result in an equation of the form:
```{math}
A \frac{\partial \psi}{\partial t} = B\frac{\partial^2 \psi}{\partial x^2}
```
where $A = i \hbar$ because $\frac{\partial \psi}{\partial t} = -i\omega \psi$ and $B = -\frac{\hbar^2}{2m}$ because 
$\frac{\partial^2 \psi}{\partial x^2} = -k^2 \psi$ producing a one dimensional free particle equation of the form:
```{math}
i \hbar\frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \psi}{\partial x^2}
```
If we want to include a potential energy $V = V(x,t)$ in the system, then we first rewrite the energy as KE + PE:
```{math}
E = \frac{p^2}{2m} + V
```
which results in a modified dispersion relation
```{math}
\hbar \omega = \frac{\hbar^2 k^2}{2m} + V
```
and so a modified wave equation of the form:
```{math}
i \hbar\frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \psi}{\partial x^2} + V\psi
```
which is known as the <b>Time Dependent Schr\"odinger Equation (TDSE)</b>.  Sometimes this is written in the form of an 
eigenvalue problem:
```{math}
 -\frac{\hbar^2}{2m}\frac{\partial^2 \psi}{\partial x^2} + V\psi= E\psi
```
which is known as the <b>Time Independent Schrodinger Equation (TISE)</b> where $E$ are the energy eigenvalues. 
