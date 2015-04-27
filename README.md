# ForBES

**ForBES** (standing for **For**ward-**B**ackward **E**nvelope **S**olver) is a MATLAB solver for
nonsmooth convex optimization problems.

Here is a performance comparison between ForBES, the fast forward-backward splitting method (also
known as fast proximal gradient method) and ADMM (alternating direction method of multipliers),
applied to a Lasso problem with 3K observations and 500K features, for a total of 7.5M nonzero coefficients.
<p align="center">
<img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/lasso_random_3e3_5e5_lambda_3e-1.png"
 width="640" height="349">
</p>

## Installation

Simply clone the git repository, or click on this [link](https://github.com/lostella/ForBES/archive/master.zip)
to download it as a zip archive and decompress the archive. Then move with the MATLAB command line to
the directory of ForBES. Compile all the *mex*-files required by simply hitting

```
> make
```

## How to use it

ForBES consists mainly of two MATLAB routines, `minfbe` and `miname`.
In order to use them one must provide a description of the problem in a MATLAB
structure and (optionally) a set of options:

```
out = minfbe(prob, opt);
out = miname(prob, opt);
```

Structure `prob` contains attributes describing the details of the problem, such as objective
terms and constraints, while `opt` describes, e.g., details on the algorithm to use, termination
criteria, the level of verbosity, and so on. In the following we describe more in detail how to define
these structures. Output `out` will contain the results of the optimization process.

Examples on how to use `minfbe` and `miname` can be found in the [tests folder](https://github.com/lostella/ForBES/tree/master/tests). Furthermore, you can access the help file of the solvers directly from MATLAB with

```
> help minfbe
> help miname
```

## minfbe (convex composite problems)

We consider here problems in the form

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/minfbe_problem.png" alt="Convex composite problem"></p>

where *f1* is convex quadratic, *l* is a linear term and *f2* is any convex, twice continuously
differentiable function with Lipschitz continuous gradient. Any of these terms may be omitted in the
problem definition, in which case it is assumed to be identically zero.
Function *g* is a general proper, closed, convex function (possibly nonsmooth).
This form of includes many practical problems arising in several fields such as optimal
control, data analysis, machine learning, image and signal processing to name a few.

Since *f1* is quadratic, it is entirely specified by its Hessian and linear parts:

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/quad_fun.png" alt="Quadratic function"></p>

The generic nonlinear term *f2* is described by an appropriate function returning its value and gradient
(in this exact order) at any specified point.
The gradient may be computed only when the corresponding output argument is requested.
For example, the logistic function

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/log_loss.png" alt="Logistic function"></p>

can be defined in MATLAB as follows

```
function [fz, gradfz] = LogReg(z)
    pz = 1./(1+exp(-z));
    fz = -sum(log(pz));
    if nargout >= 2
        gradfz = (pz-1);
    end
end
```

The nonsmooth term *g(x)* is defined through its proximal mapping and the value of
*g* at the proximal point:

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/prox_definition.png"></p>

For example, if *g(x) = r||x||_1* (the L1-norm) then it is described in MATLAB as the following soft-thresholding
procedure:

```
function [z, v] = L1Norm(x, gam, r)
    uz = max(0, abs(x)-gam*r);
    if nargout >= 2
        v = r*sum(uz);
    end
    z = sign(x).*uz;
end
```

In the following table we summarize the attributes that may be used to define problem (1). Notice that **most of them are optional**.

Attribute | Type | Mandatory? | Default | What is it
----- | ---- | ---------- | ------- | ----------
`prob.x0` | vector | yes | - | The starting point for the algorithm.
`prob.Q` <br> `prob.q` | matrix (or function) and vector | no | Q = 0 <br> q = 0 | The Hessian and linear parts of the quadratic term *f1*.
`prob.A` <br> `prob.b` | matrix (or function) and vector | no | A = Id <br> b = 0 | The affine mapping with which *f1* is composed.
`prob.AT` | function | yes, if A is defined as a function | - | The function computing the adjoint of A.
`prob.lin` | vector | no | 0 | The linear term *l* in the objective.
`prob.f2` | function | no | the zero function | A procedure returning the value of *f2* (1st output) and its gradient (2nd output) at the specified point. It can optionally return also the Hessian as 3d ouput.
`prob.useHessian` | boolean, integer | no | 0 | A flag indicating whether f2 returns also the Hessian of *f2*.
`prob.C` <br> `prob.d` | matrix or function, vector | no | C = Id <br> d = 0 | The affine mapping with which *f2* is composed.
`prob.CT` | function | yes, if C is defined as a function | - | The function computing the adjoint of C.
`prob.Lf1` | real | no | computed numerically | The 2-norm of matrix A'QA.
`prob.Lf2` | real | no | computed numerically | The Lipschitz constant of the gradient of *f2*.
`prob.normC` | real | no | computed numerically | The 2-norm of matrix C.
`prob.g` | function | yes | - | A procedure that given *x* and *gamma* (in this order) computes the proximal point of *x* (1st output) with respect to *g* and stepsize *gamma*, and the value of *g* at the proximal point (2nd output)

**Example**: using LogReg and L1Norm defined above, we can test the `minfbe` onto the sparse logistic regression problem

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/log_reg.png" alt="Sparse logistic regression"></p>

as follows:

```
prob.C = diag(sparse(b))*A;
prob.f2 = @(x) LogReg(x);
prob.g = @(x, gam) L1Norm(x,gam,reg);
prob.x0 = zeros(n,1);
out = minfbe(prob);
```

The `out` structure will contain the results of the optimization process, including the computed solution
and some additional information like the progress of the algorithm during the iterations.

## miname (convex separable problems)

We consider now problems with linear equality constraints, of the following form:

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/miname_problem.png" alt="Equality constrained convex problem"></p>

with *f1* (if present) is strongly convex and quadratic on its domain, *f2* (if present) is strongly
convex and twice continuously differentiable in the interior of its domain, while *g* is proper, closed
and convex. Linear operators *A1, A2, B* need to be specified only if the corresponding term in the 
objective is present, and they are of appropriate dimension along with vector *c* in the constraints.

The problem is described by specifying the constraint and providing appropriate procedures for computing
the primal iterates (and the corresponding objective values) given a dual variable. Specifically:

<p align="center"><img src="https://raw.githubusercontent.com/lostella/ForBES/master/figures/miname_primal_updates.png"></p>

The following table summarizes the attributes defining problem (2). Again, notice that **many of them are optional**.

Attribute | Type | Mandatory? | Default | What is it
----- | ---- | ---------- | ------- | ----------
`prob.x1step` | function | no | - | Procedure minimizing *f1(w)- y'w* with respect to *w*, given *y* (since *f1* is quadratic, this procedure computes an affine mapping).
`prob.A1` | matrix or function | yes, if `prob.x1step` is defined | - | Matrix *A1* in the constraint.
`prob.A1T` | function | yes, if `prob.A1` is defined as a function | - | Procedure computing the adjoint of *A1*.
`prob.x2step` | function | no | - | Procedure minimizing *f2(w)- y'w* with respect to *w*, given *y*.
`prob.A2` | matrix or function | yes, if `prob.x2step` is defined | - | Matrix *A2* in the constraint.
`prob.A2T` | function | yes, if `prob.A2` is defined as a function | - | Procedure computing the adjoint of *A2*.
`prob.zsetp` | function | yes | - | Procedure minimizing the augmented Lagrangian with respect to *z* and computing *g(z)*.
`prob.B` | matrix | yes | - | Matrix *B* in the constraint.
`prob.c` | vector | yes | - | The right hand side of the constraint.

## Options

Optional settings may be enabled by specifying the correspondent fields in the `opt` structure passed
as second argument to `minfbe` and `miname`.

Attribute | Type | Default | What is it
----- | ---- | ------- | ----------
`opt.tolOpt` | scalar | 1e-5 | Tolerance on the optimality condition.
`opt.maxit` | integer | 100*n* | Maximum number of iterations.
`opt.method` | string | 'lbfgs' | Algorithm to use. Can select between: <br> 'sd' (steepest descent) <br> 'lbfgs' (limited memory BFGS) <br>  'cg-desc', 'cg-prp', 'cg-dyhs' (various CG algorithms) <br> 'bb' (Barzilai-Borwein).
`opt.variant` | string | 'global' | 'basic': Use the basic algorithm<br> 'global': Use the **global** variant<br> 'fast': Use the **fast** variant
`opt.linesearch` | string | method dependant | Line search strategy to use. Can select between: <br> 'armijo' (default for 'sd') <br> 'nonmonotone-armijo' (default for 'bb') <br> 'hager-zhang' (default for the rest) <br> 'lemarechal' <br> 'fletcher'

## Credits

ForBES is developed by Lorenzo Stella [`lorenzo.stella-at-imtlucca.it`] and Panos Patrinos [`panagiotis.patrinos-at-imtlucca.it`]. Any feedback, bug report or suggestion for future improvements is more than welcome. We recommend using the [issue tracker](https://github.com/lostella/ForBES/issues) to report bugs.
