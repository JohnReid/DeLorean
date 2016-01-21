Low rank approximations
-----------------------

We define a covariance structure $K_{\tau,\tau}$ between the cells based on
their pseudotimes. The expression of each gene then follows a normal
distribution with gene-specific temporal variation $\psi_g$ and noise
$\omega_g$ $$\mathcal{N}(0, \psi_g K_{\tau,\tau} + \omega_g I)$$

Now we consider a set of inducing pseudotimes, $u$, that we will use in a
sparse approximation of this Gaussian process. Following Quinonero and
Rasmussen's notation for Snelson and Ghahramani's Sparse Gaussian Processes
using pseudoinputs, we have the approximation $$\mathcal{N}(0,
\psi_g(Q_{\tau,\tau} - \textrm{diag}[Q_{\tau,\tau} - K_{\tau,\tau})]) +
\omega_g I)$$ where $Q_{\tau,\tau} = K_{\tau,u}K_{u,u}^{-1}K_{u,\tau}$. We
choose Snelson and Ghahramanani's sparse approximation as it does not
underestimate the variance at pseudotimes away from the inducing inputs. We
prefer the variance to be overestimated to encourage the cell pseudotimes to
stay within the range of the inducing pseudotimes.

Note that we choose $|u|$ to be smaller than $|\tau|$ so that the approximation
gives us computational savings.  In particular we can use the Matrix Inversion
Lemma to efficiently calculate the precision of the sparse approximation and
hence efficiently evaluate the likelihood. Suppose that we have a $|\tau|
\times k$ matrix $A_{u,\tau}$ such that $Q_{\tau,\tau}=A_{u,\tau}^T A_{u,\tau}$
then the precision $$\bigg(\psi_g(Q_{\tau,\tau} - \textrm{diag}[Q_{\tau,\tau} -
K_{\tau,\tau})]) + \omega_g I\bigg)^{-1}$$ is a low rank update, $\psi_g
Q_{\tau,\tau}$, to the easily invertible diagonal matrix $\psi_g
\textrm{diag}[Q_{\tau,\tau} - K_{\tau,\tau})] + \omega_g I$. The Matrix
Inversion Lemma allows us to calculate this at the computational cost of an
inversion of a $k \times k$ matrix.

This $A_{u,\tau}$ is easy to find through the inverse of a Cholesky
decomposition of $K_{u,u}$ multiplied by $K_{\tau,u}$.  Note also that the
Cholesky decomposition of $K_{u,u}$ can be done ahead of time as the inducing
points are fixed (unlike the $\tau$ which are parameters to be estimated).
Note that in general $k$ may be smaller than $|u|$ when $K_{u,u}$ is not of
full rank providing further computational savings.
