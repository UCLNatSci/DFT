# Computational Methods

## Matrices

We can write the ODEs here in matrix form:
```{math}
\hat{H}\phi = \mathcal{E}\,\phi \Rightarrow H_{ij}\,\phi_j = E\,\phi_i
```

where here the Hamiltonian operator has the form:

```{math}
\hat{H} = \hat{T} + V
```

$T$ is the kinetic energy, which in differential form looks like:
```{math}
\hat{T} = -\frac{\hbar^2}{2m}\nabla^2
```

In Hartree units, we can write this in matrix form as:
```{math}
-\frac{1}{2}D^2
```
where $D$ is a differentiation operator, which follows:
```{math}
D \begin{pmatrix} f_0 \\ f_1 \\ \vdots \\ f_{N-1} \end{pmatrix} = \begin{pmatrix} f_0' \\ f_1' \\ \vdots \\ f_{n-1}'\end{pmatrix} 
```

$V$ is the potential operator, which in the differential form just multiplies the wavefunction $V\,\phi$, but in matrix form it will be some diagonal matrix:
```{math}
V_{ij} = \begin{pmatrix} V_0 & 0 & \dots & 0\\ 0 & V_1 & \dots & 0 \\ \vdots & \vdots & \ddots &  0 \\ 0 & 0 & 0 & V_{N-1}\end{pmatrix}
```


## Integrals

We often need to find an integral, for example the Hartree potential:

```{math}
V_H = \int \frac{\rho(r')}{|r-r'|}\,\mathrm{d}r'
```

we can write this in matrix form also:
```{math}
V_i = \int \frac{\rho(r_j)}{|r_i-r_j|}\,\mathrm{d}r
```

and to calculate the integral we can just apply a trapezium rule expression:

```{math}
I_i = \int f(x_i)\,\mathrm{d}x = \sum \Big[f(x_i)\,\Delta x\Big]
```

so here:
```{math}
V_i = \sum \left[\frac{\rho(r_j)\,\Delta r}{\sqrt{(r_i-r_j)^2}}\right]
```