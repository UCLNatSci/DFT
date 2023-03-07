# The Kohn-Sham equation

Although we wish to solve the Schrodinger equation here:

```{math}
\left(-\frac{\hbar^2\,\nabla^2}{2m} + V_{ext}({\bf r}) \right)\,\Psi({\bf r}) = E\,\Psi({\bf r})
```
where in Hydrogen actoms we pick $V_{ext} = -\frac{e^2}{r}$.  However once we move beyond one electron systems we find the system looking a little more complicated, 
there is a nucleus-electron attraction, electron-electron repulsion as well as electron spin considerations to add in also.  

Although we know from solving the Hydrogen atom Schrodinger equation analytically, electrons get arranged into *orbitals*, e.g. Hydrogen has an electronic configuration 
$1s^1$ and Carbon has an electronic configuration $1s^2\,2s^2\,2p^2$.  The use of orbital wavefunctions $\phi({\bf r})$ allows us to simplify the picture to solve for 
l;larger and larger atoms with many electrons.  This is known as the **Hartree approximation** and serves us well up to a point - it does not reflect the anti-symmetric 
nature of electrons (this is down to their spin).  The adjustment we make to the system is to consider the *Slater determinant*:

```{math}
\Psi({\bf r_1},s_1;\, {\bf r_2}, s_2;\,\dots) = 
\begin{vmatrix} 
\phi_1({\bf r_1},s_1) & \phi_1({\bf r_2},s_2) & \dots & \phi_1({\bf r_N},s_N) \\
\phi_2({\bf r_1},s_1) & \phi_2({\bf r_2},s_2) & \dots & \phi_2({\bf r_N},s_N) \\
\vdots & \vdots & \ddots & \vdots \\
\phi_N({\bf r_1},s_1) & \phi_N({\bf r_2},s_2) & \dots & \phi_N({\bf r_N},s_N) \
\end{vmatrix}
```

The Hartree-Fock equation is essentially the Schrodinger equation for the ** orbital ** energy levels:
```{math}
\left(-\frac{\hbar^2\,\nabla^2}{2m} + V_N({\bf r}) + V_H({\bf r}) + V_{x}^i({\bf r})\right)\,\phi_i({\bf r}) = \mathcal{E}_i\,\phi_i({\bf r})
```
wher $V_x^i({\bf r})$ is an orbital dependent term.  

The issue is that solving these one by one and then finding the total energy contriution thanks to *every* electron is very time consuming and computational impractical, thus 
we need to simplify.

We attempt to use an electron density field $\rho({\bf r})$ here to help solve these equations, this has the advantage of not treating the electgrons as separate entities, 
but rather as one large field that we can attempt to find - hence we will drop the individualk orbital exchange terms for one large overall exchange term 
$V_x^i \rightarrow V_{XC}$.  Another assumption is to think about only the **outer** electrons as the relevent ones, those in filled electrong shells (so called tight binding) do not really 
play a big role in the energies of outer electrons.

The electron charge density can be found as:
```{math}
\rho({\bf r}) = e\,\sum_{i,\,occuped}^N |\phi({bf r})|^2
```

Our Hatree potential term can be found from electromagnetism (Gauss's law in potentials):

```{math}
\nabla^2 V_H({\bf r}) = -4\pi e\,\rho({\bf r})
```

This can be solve for $V_H$ using a greens function:

```{math}
V_H(\bf r}) = \int \frac{e^2\,\rho({\bf r})}{|{\bf r - r'}|}\,r'\,\sin(\theta)\,\mathrm{d}r'\,\mathrm{\theta}\,\mathrm{d}\phi
```

One prblem with solving this in one dimension is that it won't converge - there will be a some sort of divergence in each of the terms, so we modify this expression slightly:
```{math}
V_H = \int \frac{e^2\,\rho({\bf r})}{\sqrt{|r - r'}|^2 + \epsilon}}\,r'\,\sin(\theta)\,\mathrm{d}r'\,\mathrm{theta}\,\mathrm{d}\phi
```
and add in a term $\epsilon \ll 1$.  

The exchange term looks like:
```{math}
V_{XC}(\rho({\bf r}) = -\left(\frac{3\,\rho}{\pi}\right)^{1/3}
```
this is after we applied the local density approximation, before this step it looks VERY complicated.


This leads us to the Kohn-Sham (KS) equations:

```{math}
\left(-\frac{\hbar^2}{2m} \nabla^2 + V_N({\bf r) + V_H({\bf r}) + V_{xc}(\rho{\bf r})\right)\,\phi_i({\bf r}) = \mathcal{E}_i\,\phi_i({\bf r})
```

We aim to solve this system iteratively:

```{figure} ../figures/DFTFlowChart.png
---
name: DFTFlow
---
A schematic for following the **self-consistent** field approach to solve the KS equations.
```

We all this the self-consistent field approach to solve the KS equations, 
