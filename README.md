# Root-finding
simple code to solve algebraic equation and system of equations.

## Bisection method
The bisection algorithm is basically a binary search, and is based on the zeros theorem, i.e. if a function is good enough then there is a zero. Clearly then we more or less need to know where to look because it is necessary that the selected region contains zero. Once an interval has been chosen, the midpoint is sought and the function is evaluated at that point, depending on a condition, the midpoint becomes the new extreme of the interval, and so on the interval decreases (practically the function calculated at an extreme must have opposite sign compared to the same calculated in the other extreme)
The implementation of the method is shown in bisection.py

## Secant method
The method consists in constructing a sequence of points with the following criterion: given two initial points x_0,x_1, for each n>=1 the point x_{n+1} is the zero of the line passing through the points (x_{n- 1},f(x_{n-1})),(x_{n},f(x_{n}))

<img src="https://latex.codecogs.com/svg.image?x_{n&plus;1}=&space;x_{n}&space;-&space;\frac{x_{n}-x_{n-1}}{f(x_{n})-f(x_{n-1})}&space;f(x_{n})" title="https://latex.codecogs.com/svg.image?x_{n+1}= x_{n} - \frac{x_{n}-x_{n-1}}{f(x_{n})-f(x_{n-1})} f(x_{n})" />
The implementation of the method is shown in secanti.py


## Newton method
If we consider a x_0 very close to the solution we can expand in Taylor series and obtain:
<img src="https://latex.codecogs.com/svg.image?f(s)&space;=&space;0&space;=&space;f(x_0)&space;&plus;&space;(x_0-s)\frac{df}{dx}(x_0)&space;\hspace{5&space;mm}&space;\text{from&space;which}&space;\hspace{5&space;mm}&space;s&space;=&space;x_0&space;&plus;&space;\frac{f(x_0)}{\frac{df}{dx}(x_0)}" title="https://latex.codecogs.com/svg.image?f(s) = 0 = f(x_0) + (x_0-s)\frac{df}{dx}(x_0) \hspace{5 mm} \text{from which} \hspace{5 mm} s = x_0 + \frac{f(x_0)}{\frac{df}{dx}(x_0)}" />

which leads to the iterative method:

<img src="https://latex.codecogs.com/svg.image?x_{n&plus;1}&space;=&space;x_n&space;&plus;&space;\frac{f(x_n)}{\frac{df}{dx}(x_n)}" title="https://latex.codecogs.com/svg.image?x_{n+1} = x_n + \frac{f(x_n)}{\frac{df}{dx}(x_n)}" />

Obviously, in addition to the zeros of a single function, we can also solve a system; Basically newton raphson is like newton's rule seen above, only now x is a vector and instead of the derivative we need to compute the matrix of derivatives and invert it:

<img src="https://latex.codecogs.com/svg.image?\textbf{x}_{n&plus;1}&space;=&space;\textbf{x}_n&space;-&space;J(\textbf{x}_n)^{-1}F(\textbf{x}_n)" title="https://latex.codecogs.com/svg.image?\textbf{x}_{n+1} = \textbf{x}_n - J(\textbf{x}_n)^{-1}F(\textbf{x}_n)" />

Where J is defined as:

<img src="https://latex.codecogs.com/svg.image?J&space;=&space;\begin{bmatrix}&space;\dfrac{\partial&space;f_1}{\partial&space;x_1}&space;&&space;\cdots&space;&&space;\dfrac{\partial&space;f_1}{\partial&space;x_n}&space;\\&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\dfrac{\partial&space;f_m}{\partial&space;x_1}&space;&&space;\cdots&space;&&space;\dfrac{\partial&space;f_m}{\partial&space;x_n}&space;&space;\end{bmatrix}&space;\hspace{10&space;mm}&space;J_{ij}&space;=&space;\frac{\partial&space;f_i&space;(\mathbf&space;{x})}{\partial&space;x_j}" title="https://latex.codecogs.com/svg.image?J = \begin{bmatrix} \dfrac{\partial f_1}{\partial x_1} & \cdots & \dfrac{\partial f_1}{\partial x_n} \\ \vdots & \ddots & \vdots \\ \dfrac{\partial f_m}{\partial x_1} & \cdots & \dfrac{\partial f_m}{\partial x_n} \end{bmatrix} \hspace{10 mm} J_{ij} = \frac{\partial f_i (\mathbf {x})}{\partial x_j}" />

The implementation of the method in 1D is shown in tangenti.py, while in sistemi_non_linerai.py an implementation in the 2D case is shown since it is easy to invert a 2D matrix by hand. In the last one there is also a comparison with the functions of scipy

## The Kaczmarz method

Let Ax=b be a system of linear equations, let m be the number of rows of A, a_{i} be the i-th row of complex-valued matrix A, and let x^{0} be arbitrary complex-valued initial approximation to the solution of Ax=b. The method consists in iteratively calculating:

<img src="https://latex.codecogs.com/svg.image?x^{k&plus;1}&space;=&space;x^{k}&space;&plus;&space;\frac{b_{i}&space;-&space;\langle&space;a_{i},&space;x^{k}&space;\rangle}{\|&space;a_{i}&space;\|^2}&space;\overline{a_{i}}" title="https://latex.codecogs.com/svg.image?x^{k+1} = x^{k} + \frac{b_{i} - \langle a_{i}, x^{k} \rangle}{\| a_{i} \|^2} \overline{a_{i}}" />

where i=k mod m and \overline{a_i} denotes complex conjugation a_{i}

The implementation of the method is shown in kaczmarz.py, and there is also a comparison with the functions of numpy.

## Tridiagonal matrix algorithm
The tridiagonal matrix algorithm, is a simplified form of Gaussian elimination that can be used to solve tridiagonal systems of equations:

<img src="https://latex.codecogs.com/svg.image?\begin{bmatrix}&space;&space;&space;b_1&space;&&space;c_1&space;&&space;&space;&space;&space;&space;&space;&space;&space;&&space;&space;&space;&space;&space;&space;&space;&space;&&space;&space;0&space;&space;&space;&space;&space;&space;\\&space;&space;&space;a_2&space;&&space;b_2&space;&&space;c_2&space;&space;&space;&space;&&space;&space;&space;&space;&space;&space;&space;&space;&&space;&space;&space;&space;&space;&space;&space;&space;&space;\\&space;&space;&space;&space;&space;&space;&space;&&space;a_3&space;&&space;b_3&space;&space;&space;&space;&&space;\ddots&space;&&space;&space;&space;&space;&space;&space;&space;&space;&space;\\&space;&space;&space;&space;&space;&space;&space;&&space;&space;&space;&space;&space;&&space;\ddots&space;&&space;\ddots&space;&&space;c_{n-1}&space;\\&space;&space;&space;0&space;&space;&space;&&space;&space;&space;&space;&space;&&space;&space;&space;&space;&space;&space;&space;&space;&&space;a_n&space;&space;&space;&space;&&space;b_n\end{bmatrix}\begin{bmatrix}&space;&space;&space;x_1&space;&space;&space;&space;\\&space;&space;&space;x_2&space;&space;&space;&space;\\&space;&space;&space;x_3&space;&space;&space;&space;\\&space;&space;&space;\vdots&space;\\&space;&space;&space;x_n\end{bmatrix}=\begin{bmatrix}&space;&space;&space;d_1&space;&space;&space;&space;\\&space;&space;&space;d_2&space;&space;&space;&space;\\&space;&space;&space;d_3&space;&space;&space;&space;\\&space;&space;&space;\vdots&space;\\&space;&space;&space;d_n\end{bmatrix}" title="https://latex.codecogs.com/svg.image?\begin{bmatrix} b_1 & c_1 & & & 0 \\ a_2 & b_2 & c_2 & & \\ & a_3 & b_3 & \ddots & \\ & & \ddots & \ddots & c_{n-1} \\ 0 & & & a_n & b_n\end{bmatrix}\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ \vdots \\ x_n\end{bmatrix}=\begin{bmatrix} d_1 \\ d_2 \\ d_3 \\ \vdots \\ d_n\end{bmatrix}" />

For such systems, the solution can be obtained in O(n) operations instead of O(n^{3}) required by Gaussian elimination. 
The steps are:

<img src="https://latex.codecogs.com/svg.image?c'_i&space;=&space;\begin{cases}&space;&space;\cfrac{c_i}{b_i},&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;i&space;=&space;1,&space;\\&space;&space;\cfrac{c_i}{b_i&space;-&space;a_i&space;c'_{i&space;-&space;1}},&space;&&space;i&space;=&space;2,&space;3,&space;\dots,&space;n&space;-&space;1&space;\end{cases}\hspace{10&space;mm}d'_i&space;=&space;\begin{cases}&space;&space;\cfrac{d_i}{b_i},&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&space;&&space;i&space;=&space;1,&space;\\&space;&space;\cfrac{d_i&space;-&space;a_i&space;d'_{i&space;-&space;1}}{b_i&space;-&space;a_i&space;c'_{i&space;-&space;1}},&space;&&space;i&space;=&space;2,&space;3,&space;\dots,&space;n.&space;\end{cases}" title="https://latex.codecogs.com/svg.image?c'_i = \begin{cases} \cfrac{c_i}{b_i}, & i = 1, \\ \cfrac{c_i}{b_i - a_i c'_{i - 1}}, & i = 2, 3, \dots, n - 1 \end{cases}\hspace{10 mm}d'_i = \begin{cases} \cfrac{d_i}{b_i}, & i = 1, \\ \cfrac{d_i - a_i d'_{i - 1}}{b_i - a_i c'_{i - 1}}, & i = 2, 3, \dots, n. \end{cases}" />

The solution is then obtained by back substitution

<img src="https://latex.codecogs.com/svg.image?\\x_n&space;=&space;d'_n\\x_i&space;=&space;d'_i&space;-&space;c'_i&space;x_{i&space;&plus;&space;1},&space;\quad&space;i&space;=&space;n&space;-&space;1,&space;n&space;-&space;2,&space;\ldots,&space;1" title="https://latex.codecogs.com/svg.image?\\x_n = d'_n\\x_i = d'_i - c'_i x_{i + 1}, \quad i = n - 1, n - 2, \ldots, 1" />

The implementation of the method is shown in tri_sor.f

## Successive over-relaxation 
Successive over-relaxation (SOR) is a variant of the Gauss–Seidel method:

<img src="https://latex.codecogs.com/svg.image?x^{(k&plus;1)}_i&space;&space;=&space;(1-\omega)x^{(k)}_i&space;&plus;&space;\frac{\omega}{a_{ii}}&space;\left(b_i&space;-&space;\sum_{j<i}&space;a_{ij}x^{(k&plus;1)}_j&space;-&space;\sum_{j>i}&space;a_{ij}x^{(k)}_j&space;\right),\quad&space;i=1,2,\ldots,n&space;" title="https://latex.codecogs.com/svg.image?x^{(k+1)}_i = (1-\omega)x^{(k)}_i + \frac{\omega}{a_{ii}} \left(b_i - \sum_{j<i} a_{ij}x^{(k+1)}_j - \sum_{j>i} a_{ij}x^{(k)}_j \right),\quad i=1,2,\ldots,n " />

if \omega=1 we get Gauss–Seidel

The implementation of the method is shown in tri_sor.f

## Conjugate gradient

POI LO SCRIVO
