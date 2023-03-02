# Computational Methods

Solving the Schrodinger equation for different quantum systems, we essentially add different choices of potential $V({\bf x},\,t)$.  In each case 
the Hamilonian operator:
```{math}
\hat{H} = \hat{T} + \hat{V} = -\frac{\hbar^2}{2m}\nabla^2 + V({\bf x},\,t)
```
will vary depending on the choice of potential.


## Bohr Model
This one is a bit of a cheat, because solving it does not involve the Schroodinger equation, however it does allow us to address 
some conceptual ideas in quantum theory.  To do this properly, we need to consider the hydrogen atom properly, which we do in a later section.  The 
basic idea of the Bohr model, starting with Hydrogen, is to assume that electrons are in circular motion around the nucleus:
```{math}
F_{centripetal} = \frac{m_e v^2}{r}
```
which is equal to the Coulombic attraction between the electron $-e$ and proton $+e$:
```{math}
F_{Coulomb} = \frac{e^2}{4\pi \,\epsilon_0\,r^2}
```
Thus we find the velocity of these bound electrons is found by:
```{math}
\frac{m_ev^2}{r} = \frac{e^2}{4\pi \,\epsilon_0\,r^2} \Rightarrow v^2 = \frac{e^2}{4\pi\,m_e\,\epsilon_0\,r}
```
The kinetic energy of these electrons is given by:
```{math}
E_K = -\frac{1}{2}m_e\,v^2 = -\frac{e^2}{8\pi\,\epsilon_0\,r}
```
Where the $-$ sign is introduced since these electrons are bound and we would need to give $E = \frac{1}{2}m_e\,v^2$ worth of 
energy for the electron to be free with net zero energy.  Notice that since $E_K \propto \frac{1}{r}$, this formally happens as 
$r \rightarrow \infty$.  Bohr's key insight was that the electrons do not orbit the atomic nucleus akin to planets around a star, 
but instead behave as electron matter waves.  Electrons have a de Broglie wavelength:
```{math}
\lambda = \frac{h}{mv}
```
then the key idea here is that different orbits around the nucleus would correspond to different wavelengths, with <b>only</b> complete 
wavelengths are permitted.  The circumference of each orbit would be:
```{math}
2\pi\,r = n\,\lambda,\quad n \in \mathbb{N}
```
this means the angular momentum $L = m\,v\,r$ is given by:
```{math}
m_e\,v\,r = \frac{n\,h}{2\pi} = n\,\hbar
```
Substituting in the expression for $v$ allows us to write the allowed radius of orbits:
```{math}
r_n = \frac{4\pi\,\epsilon_0\,\hbar^2}{e^2\,m_e}\,n^2
```
and hence the energy of the electrons is found to be:
```{math}
E_n = -\frac{e^4\,m_e}{2(4\pi\,\epsilon_0)^2\,\hbar^2}\frac{1}{n^2} \simeq -\frac{13.6}{n^2} \text{eV}
```
where the coefficients are collected together as the Rydeburg energy $R_E \simeq 13.6 eV$.  Recall that electron volts, eV,  
represent the amount of energy given to a particle with charge of $e = 1.6\times 10^{-19}$ C travelling across a potential 
difference of 1 V.  This expression suggests electrons can move between energy levels, e.g. the energy emitted moving from 
the $n_1$ to $n_2$ energy level would be:
```{math}
\Delta E_{n_1\rightarrow n_2} = 13.6\left(\frac{1}{{n_1}^2}-\frac{1}{{n_2}^2}\right)
```
The typical size of an atom can be represented by the ground state radius, $r_{n=1}$, which is denoted by $a_0$ and known as the Bohr radius:
```{math}
a_0 = \frac{4\pi\,\epsilon_0\,\hbar^2}{e^2\,m_e}
```
We stress that this result is <b>semi-classical</b>, it has its limitations and does not correctly account for additional quantum 
effects.  One obvious issue is how to generalise this expression to larger atoms, initially the most obvious thing to do is adjust 
the size of the atomic nuclear:
```{math}
F_{Coulomb} = \frac{e^2}{4\pi\,\epsilon_o\,r^2} \rightarrow \frac{Z\,e^2}{4\pi\,\epsilon_o\,r^2}
```
where $Z$ is the atomic number of the element in question.  The biggest issue with this approximation is that it leaves out the $e-e$ 
interactions within each energy level, something which is not typically present in atomic Hydrogen.



## Infinite Square Well
One of the simplest quantum systems we can investigate is the infinite square well, where the potential energy of the system traps the wave 
function within a well, as shown in {numref}`InfiniteSquareWell`.
```{figure} ./figures/InfiniteSquareWell.png
---
name: InfiniteSquareWell
---
The infinite square well system
``` 

To solve the Schrodinger equation for this system:
```{math}
-\frac{\hbar^2}{2m}\frac{\partial^2 \Psi}{\partial x^2} + V\,\Psi = i \hbar \frac{\partial \Psi}{\partial t}
```
we can start by using a separable ansatz:

```{math}
\Psi(x,\,t) = f(x)\,g(t)
```

and therefore drop the partial derivatives and set the problem up as two ODE problems with eigenvalues $E$:

```{math}
-\frac{\hbar^2}{2m}\frac{\mathrm{d}^2 f}{\mathrm{d} x^2} + V\,f &=&\, E\,f \\
i \hbar \frac{\mathrm{d} g}{\mathrm{d} t} &=&\, E\,g
```
Inside the well, $V = 0$ which simplifies the $x$ ODE further:
```{math}
\frac{\mathrm{d}^2 f}{\mathrm{d} x^2} + \frac{2mE}{\hbar^2}\,f = 0
```
and therefore has solutions:
```{math}
f(x) = A\,\sin\left(\frac{\sqrt{2mE}}{\hbar}x\right) + B\,\cos\left(\frac{\sqrt{2mE}}{\hbar}x\right)
```
At the boundaries of the well, $x = 0,\,L$, we find an infinite potential $V$ which is usually problematic for a physical system, 
however if we make $\Psi(0,t) = \Psi(L,t) = 0$, then this can be fixed and this gives us the relevant boundary conditions:
```{math}
f(0) &=& A(0) + B(1) = 0 \Rightarrow B = 0\\
f(L) &=& A\,\sin\left(\frac{\sqrt{2mE}}{\hbar}L\right) + B\,\cos\left(\frac{\sqrt{2mE}}{\hbar}L\right) = 0
```
and hence
```{math}
\sin\left(\frac{\sqrt{2mE}}{\hbar}L\right) = 0  \Rightarrow \frac{\sqrt{2mE}}{\hbar}L = n\pi,\,n \in \mathbb{N}
```
which means that there are energies $E$ of the form:
```{math}
E_n = \frac{\hbar^2\,\pi^2\,n^2}{2m\,L^2},\,n \in \mathbb{N}
```
and the expressions for the different wave function solutions are given by:
```{math}
\frac{\sqrt{2mE}}{\hbar} = \frac{n\pi}{L} \Longrightarrow f(x) = A\,\sin\left(\frac{n\,\pi\,x}{L}\right),\,n \in \mathbb{N}
```
and we can plot these as shown in {numref}`InfiniteSquareSolns`
```{figure} ./figures/WaveFunctionSolns.png
---
name: InfiniteSquareSolns
---
Wave functions for the infinite square well for the first three energy levels.
```
We see that this system is simply just the standing wave system with the string fixed at both ends!  Rearranging we see that 
$L = n \lambda/2, \,n \in \mathbb{N}$.  The temporal ODE is of the form:
```{math}
\frac{\mathrm{d} g}{\mathrm{d} t} + \frac{i\,E}{\hbar}\,g = 0
```
which has solution:
```{math}
g = C\,\exp\left(-i\,\frac{E}{\hbar}\,t \right)
```
where we found the form of $E_n$ from the spatial solutions, therefore the full wave function solutions are:
```{math}
\Psi(x,\,t) &=& A'\,\exp\left(-i\,\frac{E_n}{\hbar}\,t \right)\sin\left(\frac{n\,\pi\,x}{L}\right)\, 0\leq x \leq L\\ 
\Psi(x,\,t) &=& \text{otherwise}\\
E_n &=& \frac{\hbar^2\,\pi^2\,n^2}{2m\,L^2},\,n \in \mathbb{N}
```
We can find $A'$ from the necessity that wave functions are normalised, since 
```{math}
\int_{-\infty}^{\infty}\psi\,\psi^*\,\mathrm{d} x = 1
```
to make the function a probability distribution, that would fix:
```{math}
(A')^2\int_0^L  \sin^2\left(\frac{n\,\pi\,x}{L}\right)\,\mathrm{d} x = 1 
```
Making use of the double angle formula $\cos(2A) = 1 -2\sin^2(A)$ to make the integration doable, we find:
```{math}
\frac{(A')^2}{2}\left[x - \frac{L}{2n\,\pi}\cos\left(\frac{2 n\,\pi\,x}{L}\right) \right]_0^L = 
\frac{(A')^2}{2}L = 1 \Rightarrow A' = \sqrt{\frac{2}{L}}
```
and therefore the wavefunction over all space is given by:
```{math}
\Psi(x) = 
\left\{ 
\begin{array}{ll}
\sqrt{\frac{2}{L}}\,\exp\left(-\frac{i\,E_n}{\hbar}\,t \right)\sin\left(\frac{n\,\pi\,x}{L}\right) & 0 \leq x \leq L \\
0 & x < 0,\, x > L
\end{array}
\right.
```
## Finite Square Well
If the system now has a finite square well, $V = U_0$, then <em>is</em> possible to have wavefunction solutions outside of 
$0 \leq x \leq L$, lets consider the Schrodinger equation again:
```{math}
-\frac{\hbar^2}{2m}\frac{\mathrm{d}^2 f}{\mathrm{d} x^2} + V\,f = E\,f \longrightarrow \frac{\mathrm{d}^2 f}{\mathrm{d} x^2} + 
fix\frac{2m(E-U_0)}{\hbar^2} f = 0
```
which for lower energy solutions $E < U_0$ means that:

```{math}
\frac{\mathrm{d}^2 f}{\mathrm{d} x^2} - \frac{2m|U_0-E|}{\hbar^2} f = 0
```

which can be solved to find:

```{math}
f(x) = A\exp\left(-\sqrt{\frac{2m|U_0-E|}{\hbar}}x\right) + B\exp\left(\sqrt{\frac{2m|U_0-E|}{\hbar}}x\right)
```

since we do not find wave functions growing as $x \rightarrow \infty$, $B = 0$, meaning that 
```{math}
f(x) = A\exp\left(-\sqrt{\frac{2m|U_0-E|}{\hbar}}x\right)
```
which are exponentially decaying solutions - this is quantum tunnelling!  We can see a pictorial representation of these 
finite well solutions alongside the infinite well solutions in {numref}`FiniteSquareWellSolns`.  To find exact expressions 
for the wave functions in the finite well case is not so easy however, as $\Psi(x,\,t)$ is expected to be continuous.

```{figure} ./figures/FiniteWellSolns1.png
---
name: FiniteSquareWellSolns
---
Infinite (<b>left pane</b>) and finite (<b>right pane</b>) square well wave functions.
```

## Quantum Harmonic Oscillator
This system is the quantum analogue of the classical simple harmonic oscillator, the Hamiltonian is given by:
```{math}
\hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + \frac{1}{2}m\omega^2\,x^2 
```
We can solve this system to find wave-functions of the form
```{math}
\psi _{n}(x) = {\frac {1}{\sqrt {2^{n}\,n!}}}\cdot \left({\frac {m\omega }{\pi \hbar }}\right)^{1/4}\cdot 
e^{-{\frac {m\omega x^{2}}{2\hbar }}}\cdot H_{n}{\Bigl (}{\sqrt {\frac {m\omega }{\hbar }}}x{\Bigr )},\qquad n=0,1,2,\ldots
```
Where the functions $H_n$ are known as Hermite polynomials, which satisfy:

```{math}
H_{n}(z)=(-1)^{n}~e^{z^{2}}{\frac {d^{n}}{dz^{n}}}\left(e^{-z^{2}}\right)
```
The corresponding energy levels for this system are:
```{math}
E_{n}=\hbar \omega {\bigl (}n+{\tfrac {1}{2}}{\bigr )}
```

## Revisiting the Hydrogen Atom

We can put in a Coulomb law's potential into the TISE as a model for the Hydrogen atom, this gives:
```{math}
\left(-\frac{\hbar^2}{2\mu}\nabla^2 - \frac{e^2}{4\pi\,\epsilon_0\,r} \right)\psi(r,\,\theta,\,\phi) = E \,\psi(r,\,\theta,\,\varphi)
```
where $\mu$ is the {\it reduced mass} of the system, which takes into account the slight influence of the electron mass $m_e$ on the proton mass $M_p$, the centre of mass of the system is shifted slightly:
```{math}
\mu = \frac{m_e\,M_p}{m_e + M_p} = m_e\left(1 + \frac{m_e}{M_p}\right)^{-1} \simeq m_e - \dots 
```
If we expand the Laplacian $\nabla^2$ in spherical polar coordinates, this produces a PDE:
```{math}
-\frac{\hbar^2}{2\mu}\frac{1}{r^2}\left[ \frac{\partial}{\partial r}\left( r^2\,\frac{\partial\psi}{\partial r}\right) + \frac{1}{\sin\theta}\frac{\partial}{\partial \theta}\left(\sin\theta\, \frac{\partial \psi}{\partial \theta} \right) + \frac{1}{\sin^2\theta}\frac{\partial^2 \psi}{\partial \varphi^2}\right] - \frac{e^2}{4\pi\,\epsilon_0\,r} \psi = E \,\psi
```
which turns out to be separable using the ansatz:
```{math}
\psi(r,\,\theta,\varphi) = R(r)\,\,\Theta(\theta)\,\,\Phi(\varphi)
```
which reduces the PDE into three ODE's:
```{math}
\frac{\mathrm{d} }{\mathrm{d} r}\left( r^2\frac{\mathrm{d} R}{\mathrm{d} r}\right) + \frac{2\mu\,r^2}{\hbar^2}\left(E + \frac{Z\,e^2}{4 \pi \epsilon_0 \,r} \right)R - A\,R &=& 0\\
\frac{\sin \theta}{\Theta}\frac{\mathrm{d} }{\mathrm{d} \theta}\left( \sin\theta\frac{\mathrm{d} \Theta}{\mathrm{d} \theta} \right) + A\sin^2 \theta - B &=& 0 \\
\frac{1}{\Phi}\frac{\mathrm{d}^2 \Phi}{\mathrm{d} \varphi} + B &=& 0
```
where $A,\,B$ are the separation constants for this system.  This system can be solved using a combination of integrating factor for the 
$R(r)$ functions and <b>Spherical Harmonics</b> for $\Theta(\theta),\,\Phi(\varphi)$, which are complete functions that can be used to 
decompose the surface of a sphere.  The overall solutions are:
```{math}
\psi _{n\ell m}(r,\theta ,\varphi ) = {\sqrt {{\left({\frac {2}{n\,{a_0}^{*}}}\right)}^{3}\,\,{\frac {(n-\ell -1)!}{2n(n+\ell )!}}}}\,
e^{-\rho /2}\rho ^{\ell }\,{L_{n-\ell -1}}^{2\ell +1}(\rho )\,Y_{\ell }^{m}(\theta ,\varphi )
```
Where the radial coordinate $r$ is best represented by a dimensionless parameter $\rho$:
```{math}
\rho = \frac{2r}{n{a_0}^{*}}
```
and the parameter ${a_0}^{*}$ is known as the reduced Bohr radius:
```{math}
{a_0}^{*} = \frac{4\pi \epsilon_0\,\hbar^2}{\mu\,e^2}
```
taking into account the reduced mass $\mu$.  
The function ${L_{n-\ell -1}}^{2\ell +1}(\rho )$ is a generalized Laguerre polynomial of degree $ n-\ell -1$, and the function 
$Y_{\ell}^{m}(\theta,\,\varphi)$ is a spherical harmonic function of degree $\ell$ and order $m$.  The parameters $\ell,\,m,\,n$ 
are known as quantum numbers and can take the following values:

- $ n = 1,2,3,\ldots $ this is known as the <b>principal</b> quantum number.

- $ \ell = 0,1,2,\ldots ,n-1 $ this is known as the <b>azimuthal</b> quantum number.
 
- $ m = -\ell ,\ldots ,\ell $ this is known as the <b>magnetic</b> quantum number.

Additionally, the wave-functions are both normalized and orthogonal:
```{math}
\int _{0}^{\infty }r^{2}\,\mathrm{d} r\int _{0}^{\pi }\sin \theta \,\mathrm{d} \theta \int _{0}^{2\pi }\mathrm{d} \varphi 
\,\,\psi _{n\ell m}^{*}(r,\,\theta,\,\varphi)\,\,\psi_{n'\ell 'm'}(r,\theta ,\varphi )= \delta_{n\,n'}\delta_{\ell \,\ell '}\delta_{m\,m'}
```

## A note on units
It is usually easier to use <b>atomic units</b> when considering the computing systems, this saves memory storing excessively large or 
small values and only converting values at the end,  In Hartree atomic units, we set $\hbar = m_e = a_0 = e = 1$ and energies measured in
 atomic units can be converted into electron volts using the Hartree energy:
```{math}
E_{h} = m_e\left(\frac{e^2}{4\pi\,\epsilon_0\,\hbar} \right)^2 \simeq 27.211\dots \text{eV}
```
If we consider the Hydrogen atom Hamiltonian:
```{math}
\hat{H} = - \frac{\hbar^2}{2m}\nabla^2 - \frac{1}{4\pi\,\epsilon_0}\frac{e^2}{r} \rightarrow -\frac{1}{2} \nabla^2  - \frac{1}{r}
```
