# The Kohn-Sham equation

Although we wish to solve the Schrodinger equation here:

```{math}
\left(-\frac{\hbar^2\,\nabla^2}{2m} + V_{ext}({\bf r}) \right)\,\Psi({\bf r}) = E\,\Psi({\bf r})
```
where in Hydrogen actoms we pick $V_{ext} = -\frac{e^2}{r}$.  However once we move beyond one electron systems we find the system looking a little more complicated, there is a 
nucleus-electron attraction, electron-electron repulsion as well as electron spin considerations to add in also.  

Although we know from solving the Hydrogen atom Schrodinger equation analytically, electrons get arranged into *orbitals*, e.g. Hydrogen has an electronic configuration $1s^1$ and Carbon has an electronic configuration $1s^2\,2s^2\,2p^2$.  The use of orbital wavefunctions $\phi({\bf r})$ allows us to simplify the picture to solve for larger and larger atoms with many electrons.  This is known as the **Hartree approximation** and serves us well up to a point - it does not reflect the anti-symmetric nature of electrons (this is down to their spin).  The adjustment we make to the system is to consider the *Slater determinant*:

```{math}
\Psi({\bf r_1},s_1;\, {\bf r_2}, s_2;\,\dots) = 
\begin{vmatrix} 
\phi_1({\bf r_1},s_1) & \phi_1({\bf r_2},s_2) & \dots & \phi_1({\bf r_N},s_N) \\
\phi_2({\bf r_1},s_1) & \phi_2({\bf r_2},s_2) & \dots & \phi_2({\bf r_N},s_N) \\
\vdots & \vdots & \ddots & \vdots \\
\phi_N({\bf r_1},s_1) & \phi_N({\bf r_2},s_2) & \dots & \phi_N({\bf r_N},s_N) \
\end{vmatrix}
```

The Hartree-Fock equation is essentially the Schrodinger equation for orbital energy levels:
```{math}
\left(-\frac{\hbar^2\,\nabla^2}{2m} + V_N({\bf r}) + V_H({\bf r}) + V_{x}^i({\bf r})\right)\,\phi_i({\bf r}) = \mathcal{E}_i\,\phi_i({\bf r})
```
wher $V_x^i({\bf r})$ is an orbital dependent term.  


This leads us to the Kohn-Sham (KS) equations:
