# The Kohn-Sham equation

Although we wish to solve the Schrodinger equation here:

```{math}
\left(-\frac{\hbar^2\,\nabla^2}{2m} + V_{ext}({\bf r}) \right)\,\Psi({\bf r}) = E\,\Psi({\bf r})
```
where in Hydrogen actoms we pick $V_{ext} = -\frac{e^2}{r}$.  However once we move beyond one electron systems we find the system looking a little more complicated, there is a 
nucleus-electron attraction, electron-electron repulsion as well as electron spin considerations to add in also.  

Although we know from solvin the Hydrogen atom Schrodinger equation analytically, 

The Hartree-Fock equation is essentially the Schrodinger equation for orbital energy levels:
```{math}
\left(-\frac{\hbar^2\,\nabla^2}{2m} + V_N({\bf r}) + V_H({\bf r}) + V_{x}^i({\bf r})\right)\,\phi_i({\bf r}) = \mathcal{E}_i\,\phi_i({\bf r})
```
wher $V_x^i({\bf r})$ is an orbital dependent term 


This leads us to the Kohn-Sham (KS) equations:
