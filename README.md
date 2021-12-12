# Convex Clustering with ADMM implementation via R package

## Description of Package

This R package implements ADMM algorithm to solve the convex clustering problem. Given n data points \(x_i\), \(i=1,\dots, n\), in p dimensions, \(x_i\in \mathbb{R}^p\), the key idea behind the convex clustering model is that if two observations \(x_i\) and \(x_j\) belong to the same cluster, then their corresponding centroids \(u^*_i\) and \(u^*_j\) should be the same. The optimal solution \(U^*\) = [\(u^*_i\),...,\(u^*_n\)].

 It aims to minimize
\[
\text{minimize}_{U \in \mathbb{R}^{n \times p}}\frac{1}{2}\sum_{i=1}^n \|x_i - u_i\|_2^2 + \gamma \sum_{i<j}w_{ij}\|u_i - u_j\|
\]

where 

\[
w_{ij} = exp(-mu \times\|x_i - x_j\|_2^2)
\]

are the non-negative gaussian kernel weights set between the ith and jth points used for implementation in the package. \(gamma\) is the positive tuning parameter, \(u_i\) is the cluster centroid which is \(u^{th}\) column of matrix U associated to the point \(x_i\).The first term is the fidelity term while the second term is the regularization term to penalize the differences
between different centroids in order to make sure that centroids for observations in
the same cluster should be identical.

## Algorithm's Implementation

To minimize the above objective, ADMM algorithm which employs variable splitting in order to account for the shrinkage penalties in the convex clustering problem.

The minimization problem is made into a equivalent constrained problem as follows:

\[
\text{minimize}\frac{1}{2}\sum_{i=1}^n \|x_i-u_i\|_2^2 + \gamma\sum_{l\in\epsilon}w_l\|v_l\|
\]

subject to 

\[
u_{l1} - u_{l2} - v_l = 0
\]

where l is a centroid pair (\(l1\),\(l2\)) and \[\epsilon = \{l = (l1,l2) : w_l > 0\}  \]

# ADMM updates

**Augmented Lagrangian:** given $\tau > 0$ 

$$
L_{\tau}(U, V, H) = \frac{1}{2}\sum_{i=1}^n \|x_i-u_i\|_2^2 + \gamma\sum_{l\in\epsilon}w_l\|v_l\| + \sum_{l\in\epsilon}(H_l,v_l - u_{l1} + u_{l2}) + \frac{\tau}{2}\sum_{l\in\epsilon}(\| v_l - u_{l1} + u_{l2}\|_2)
$$


**Update of U:**
$$
U_i = \frac{1}{1 + n\tau}y_i +\frac{1}{1 + n\tau}\bar{x}
$$

where 

$$
y_i = x_i + \sum_{l1 = i}[H_l^t + \tau V_l^t] - \sum_{l2 = i}[H_l^t + \tau V_l^t]
$$

$$
U^{t+1} = \frac{1}{1 + n\tau}Y +\frac{1}{1 + n\tau}\bar{x}
$$

**Update of V:**

$$
V_l^{t+1} = \arg\min_{V_l}\{\frac{\gamma  w_l}{\tau}\|v_l\| + \frac{1}{2\tau}(\|V_l - (u_{l1} - u_{l2} - \frac{H_l}{\tau})\|_2^2\}
$$

$$
V_l^{t+1} = prox_\frac{\gamma  w_l}{\tau}(u_{l1}^{t+1} - u_{l2}^{t+1} + \frac{H^{t}}{\tau})
$$

**Update of H:** 
$$
H_l^{t+1} = H_l^t + \tau(v_l^{t+1} - u_{l1}^{t+1} + u_{l2}^{t+1})
$$

The corresponding references are:

Defeng Sun, K. Toh, Yancheng Yuan, "Convex Clustering: Model, Theoretical Guarantee and Efficient Algorithm", Computer Science, Mathematics J. Mach. Learn. Res, 4 October 2018 

E. C. Chi and K. Lange. Splitting methods for convex clustering. J. Computational and
Graphical Statistics, 24(4):994–1013, 2015.

Huangyue Chen, Lingchen Kong, Yan Li, "A Novel Convex Clustering Method for High-Dimensional Data Using Semiproximal ADMM", Mathematical Problems in Engineering, vol. 2020, Article ID 9216351, 12 pages, 2020.

This function takes n × p data matrix X, (optional) gamma, (optional) initial n × p matrix of cluster centers U, (optional) initial n × p matrix of cluster adjacency matrix V, (optional) initial n × p matrix of starting value of eta and (optional) eps to monitor convergence. It returns a list of U, V and eta set of matrices at convergence.

## Installation

devtools::install_github("sowmyakolluru/ConvexCluster")



