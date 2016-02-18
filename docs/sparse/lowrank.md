# Low rank approximations

We define a covariance structure $K_{\tau,\tau}$ between the cells based on
their pseudotimes. The expression of each gene then follows a normal
distribution with gene-specific temporal variation $\psi_g$ and noise
$\omega_g$ $$\mathcal{N}(0, \psi_g K_{\tau,\tau} + \omega_g I)$$

Now we consider a set of inducing pseudotimes, $u$, that we will use in a
sparse approximation of this Gaussian process. Following Quinonero and
Rasmussen's notation for Snelson and Ghahramani's Sparse Gaussian Processes
using pseudoinputs, we have the approximation $$\mathcal{N}(0,
\psi_g(Q_{\tau,\tau} + \textrm{diag}[K_{\tau,\tau} - Q_{\tau,\tau}]) +
\omega_g I)$$ where $Q_{\tau,\tau} = K_{\tau,u}K_{u,u}^{-1}K_{u,\tau}$. We
choose Snelson and Ghahramanani's sparse approximation as it does not
underestimate the variance at pseudotimes away from the inducing inputs. We
prefer the variance to be overestimated to encourage the cell pseudotimes to
stay within the range of the inducing pseudotimes.

## Low rank calculations

Note that we choose $M=|u|$ to be smaller than $C=|\tau|$ so that the
approximation gives us computational savings.  In particular we can use the
Matrix Inversion Lemma to efficiently calculate the precision of the sparse
approximation and hence efficiently evaluate the likelihood. Suppose that we
have a $M \times C$ matrix $A_{u,\tau}$ such that $Q_{\tau,\tau}=A_{u,\tau}^T
A_{u,\tau}$ then the precision $$\bigg[\psi_g(Q_{\tau,\tau} +
\textrm{diag}[K_{\tau,\tau} - Q_{\tau,\tau}]) + \omega_g I\bigg]^{-1}$$ is the
inverse of a low rank update, $\psi_g A_{u,\tau}^T A_{u, \tau}$, to the easily invertible
diagonal matrix $B_g = \omega_g I + \psi_g \textrm{diag}[K_{\tau,\tau} -
Q_{\tau,\tau}]$.  The Matrix Inversion Lemma allows us to calculate this at
the $\mathcal{O}(M^3)$ computational cost of an inversion of a $M \times M$
matrix and some $\mathcal{O}(C^2 M)$ matrix multiplications:
$$
B_g^{-1} -
B_g^{-1} A_{u,\tau}^T
    \bigg[
        \frac{1}{\psi_g} I_m + A_{u,\tau} B_g^{-1} A_{u,\tau}^T
    \bigg]^{-1}
A_{u,\tau} B_g^{-1} $$
The Matrix Determinant Lemma allows us to cheaply
calculate the determinant of the covariance as $$ \textrm{det}(I_m+\psi_g
A_{u,\tau} B_g^{-1} A^T_{u,\tau}) \textrm{det}(B_g) $$

This $A_{u,\tau}$ is easy to find through the inverse of a Cholesky
decomposition of $K_{u,u}$ multiplied by $K_{\tau,u}$. If $K_{u,u}=R R^T$ then
$$A_{u,\tau} = R^{-1} K_{u,\tau}$$ Note also that the Cholesky decomposition of
$K_{u,u}$ and its inverse can be done ahead of time as the inducing points are
fixed (unlike the $\tau$ which are parameters to be inferred).

## Predictive mean and variance

The predictive mean for gene $g$ is
$$\frac{\psi_g}{\omega_g} K_{\star,u} \Sigma_g K_{u,\tau} x$$
and the predictive variance is
$$\psi_g(K_{\star,\star} - Q_{\star,\star} + K_{\star,u} \Sigma_g K_{u,\star})$$
where
$$\Sigma_g = \bigg[\frac{\psi_g}{\omega_g} K_{u,\tau} K_{\tau,u} + K_{u,u}\bigg]^{-1}$$
