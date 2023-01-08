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

The implementation of the method is shown in sistemi_lineari.py, and there is also a comparison with the functions of numpy.
