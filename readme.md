# Numerical Methods Console Application

## Overview

This project is a comprehensive collection of numerical methods implemented in C++. It covers fundamental algorithms for solving linear and non-linear equations, interpolation, numerical differentiation and integration, differential equations, and curve fitting techniques.

---

## Project Structure

```text
rootX/
├── readme.md
│
├── A. Solution of Linear Equations/
│   ├── Gauss Elimination/
│   ├── Gauss Jordan/
│   ├── LU Decomposition/
│   ├── Matrix Inverse/
│   └── Iterative Methods/
│       ├── Jacobi/
│       └── GaussSeidel/
│
├── B. Solution of Non-Linear Equations/
│   ├── Bisection/
│   ├── False_Position/
│   ├── Secant/
│   └── Newton_Raphson/
│
├── C. Interpolation and Approximation/
│   ├── Newton Forward/
│   └── Newton Backward/
│
├── E. Solution of Differential Equations/
│   ├── Newton Forward Differentiation/
│   └── Runge Kutta/
│
├── F. Numerical Integration/
│   ├── Simpson 1 by 3/
│   └── Simpson 3 by 8/
│
└── G. Curve Fitting/
    ├── Linear/
    ├── Polynomial/
    └── Transcendental/
```

---

## Table of Contents

### A. Solution of Linear Equations

1.  [Gauss Elimination Method](#gauss-elimination)
2.  [Gauss Jordan Method](#gauss-jordan)
3.  [LU Decomposition Method](#lu-decomposition)
4.  [Jacobi Method](#jacobi)
5.  [GaussSeidel Method](#gaussseidel)

### B. Solution of Non-Linear Equations

1.  [Bisection Method](#bisection)
2.  [False Position Method](#false_position)
3.  [Secant Method](#secant)
4.  [Newton Raphson Method](#newton_raphson)

### C. Interpolation and Approximation

1.  [Newton Forward Interpolation](#newton-forward)
2.  [Newton Backward Interpolation](#newton-backward)

### E. Solution of Differential Equations

1.  [Newton Forward Differentiation](#newton-forward-differentiation)
2.  [Runge Kutta Method](#runge-kutta)

### F. Numerical Integration

1.  [Simpson 1 by 3](#simpson-1-by-3)
2.  [Simpson 3 by 8](#simpson-3-by-8)

### G. Curve Fitting Methods

1.  [Linear Equation](#linear)
2.  [Polynomial Equation](#polynomial)
3.  [Transcendental Equation](#transcendental)

---

## A. Solution of Linear Equations

<a id="gauss-elimination"></a>

### 1. Gauss Elimination Method

<a id="gauss-jordan"></a>

### 2. Gauss Jordan Method

<a id="lu-decomposition"></a>

### 3. LU Decomposition Method

**Theory**

1. **The main idea**: Given a square matrix (**A**), it will be rewritten as $A = LU$, where **L** is a lower triangular matrix and **U** is an upper triangular matrix.
2. **Using LU to solve AX = B**:
   Once we have $A = LU$, we can write $LUx = b$.
   We introduce a helper vector **Y**: $Ly = b$ and $Ux = y$.
   First, solve $Ly = b$ using **forward substitution**.
   Then, solve $Ux = y$ using **backward substitution**.
3. **Why pivoting is often needed**:
   Sometimes LU decomposition runs into trouble if a diagonal element (called a **pivot**) becomes zero. To fix this, we often swap rows to bring a better pivot into position. $PA = LU$. This is called **partial pivoting**.

**Algorithm**
For each column/step ($k = 1$) to ($n$):

1.  Compute the ($k$)-th row of (**U**)
2.  Compute the ($k$)-th column of (**L**)

Mathematically:
$$U_{k,j} = A_{k,j} - \sum_{s=1}^{k-1} L_{k,s} U_{s,j} \quad , \quad j = k, \dots, n$$
$$L_{i,k} = \frac{A_{i,k} - \sum_{s=1}^{k-1} L_{i,s} U_{s,k}}{U_{k,k}} \quad , \quad i = k+1, \dots, n$$

**Pseudocode**

```cpp
function solve_LU(Matrix A, Vector b):
    n = size(A)
    L = IdentityMatrix(n)
    U = ZeroMatrix(n)

    for i from 0 to n-1:
        pivot_row = find_max_in_column(i, from: i to n-1)
        swap_rows(A, i, pivot_row)
        swap_rows(b, i, pivot_row)

        for k from i to n-1:
            prod_sum = 0
            for j from 0 to i-1:
                prod_sum += L[i][j] * U[j][k]
            U[i][k] = A[i][k] - prod_sum

        for k from i to n-1:
            if i == k:
                L[i][i] = 1
            else:
                prod_sum = 0
                for j from 0 to i-1:
                    prod_sum += L[k][j] * U[j][i]
                L[k][i] = (A[k][i] - prod_sum) / U[i][i]

    for i from 0 to n-1:
        if abs(U[i][i]) < epsilon:
            return "No unique solution"

    Vector y(n)
    for i from 0 to n-1:
        sum = 0
        for j from 0 to i-1:
            sum += L[i][j] * y[j]
        y[i] = b[i] - sum

    Vector x(n)
    for i from n-1 down to 0:
        sum = 0
        for j from i+1 to n-1:
            sum += U[i][j] * x[j]
        x[i] = (y[i] - sum) / U[i][i]

    return x
```

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/A.%20Solution%20of%20Linear%20Equations/LU%20Decomposition)

Includes:

- C++ source code
- Input file
- Output file

**Further Study**

- [LU Decomposition - Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
- [LU Decomposition Method - GeeksforGeeks](https://www.geeksforgeeks.org/l-u-decomposition-system-linear-equations/)
- Numerical Methods for Engineers - Chapra & Canale

<a id="iterative-methods"></a>

### 4. Iterative Methods: Jacobi and Gauss-Seidel methods

<a id="jacobi"></a>

#### (i) Jacobi Iterative Method

<a id="gaussseidel"></a>

#### (ii) Gauss-Seidel Iterative Method

---

## B. Solution of Non-linear Equations

<a id="bisection"></a>

### 1. Bisection Method

<a id="false_position"></a>

### 2. False Position Method

<a id="secant"></a>

### 3. Secant Method

<a id="newton_raphson"></a>

### 4. Newton Raphson Method

## C. Interpolation and Approximation

<a id="newton-forward"></a>

### 1. Newton Forward Interpolation Method

**Theory**
Newton's Forward Interpolation is used to approximate the value of a function $f(x)$ at valid points, based on a set of known data points that are **equally spaced**.

The interpolation formula is given by:
$$y(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \frac{u(u-1)(u-2)}{3!} \Delta^3 y_0 + \dots$$
Where $u = \frac{x - x_0}{h}$.

**Algorithm**

1.  **Input**: Read ($n$) data points ($x, y$) and the value to interpolate ($value$).
2.  **Difference Table**: Construct the forward difference table.
3.  **Calculate u**: Compute $u = (value - x[0]) / (x[1] - x[0])$.
4.  **Compute Sum**: Apply formula.
5.  **Output**: Display the interpolated value.

**Pseudocode**

```cpp
function newton_forward(x[], y[][], n, value):
    for i from 1 to n-1:
        for j from 0 to n-i-1:
            y[j][i] = y[j+1][i-1] - y[j][i-1]

    sum = y[0][0]
    u = (value - x[0]) / (x[1] - x[0])

    for i from 1 to n-1:
        u_term = 1
        for k from 0 to i-1:
            u_term = u_term * (u - k)

        factorial = fact(i)
        sum = sum + (u_term * y[0][i]) / factorial

    return sum
```

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/C.%20Interpolation%20and%20Approximation/Newton%20Forward)

Includes:

- C++ source code
- Input file
- Output file

**Further Study**

- [Newton Forward Interpolation - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- Numerical Methods for Engineers - Chapra & Canale

<a id="newton-backward"></a>

### 2. Newton Backward Interpolation Method

**Theory**
Newton's Backward Interpolation is similar to the forward method but is more accurate for interpolating values near the **end** of the dataset.

The formula is given by:
$$y(x) = y_n + u \nabla y_n + \frac{u(u+1)}{2!} \nabla^2 y_n + \frac{u(u+1)(u+2)}{3!} \nabla^3 y_n + \dots$$
Where $u = \frac{x - x_n}{h}$.

**Algorithm**

1.  **Input**: Read ($n$) data points ($x, y$) and the values.
2.  **Difference Table**: Construct the backward difference table.
3.  **Calculate u**: Compute $u = (value - x[n-1]) / (x[1] - x[0])$.
4.  **Compute Sum**: Apply formula.
5.  **Output**: Display result.

**Pseudocode**

```cpp
function newton_backward(x[], y[][], n, value):
    for i from 1 to n-1:
        for j from n-1 down to i:
             y[j][i] = y[j][i-1] - y[j-1][i-1]

    sum = y[n-1][0]
    u = (value - x[n-1]) / (x[1] - x[0])

    for i from 1 to n-1:
        u_term = 1
        for k from 0 to i-1:
             u_term = u_term * (u + k)

        factorial = fact(i)
        sum = sum + (u_term * y[n-1][i]) / factorial

    return sum
```

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/C.%20Interpolation%20and%20Approximation/Newton%20Backward)

Includes:

- C++ source code
- Input file
- Output file

**Further Study**

- [Newton Backward Interpolation - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- Numerical Methods for Engineers - Chapra & Canale

---

## E. Solution of Differential Equations

<a id="newton-forward-differentiation"></a>

### 1. Newton Forward Differentiation

**Theory**

Newton's Forward Differentiation formula is used to compute the derivative of a function at a given point using a forward difference table. It is based on Newton's Forward Interpolation formula.

For equally spaced data points, the first derivative is approximated by:

$$f'(x) = \frac{1}{h} \left[ \Delta y_0 - \frac{2u-1}{2!} \Delta^2 y_0 + \frac{3u^2 - 6u + 2}{3!} \Delta^3 y_0 - \cdots \right]$$

Where:

- $h$ is the step size (interval between consecutive x values)
- $u = \frac{x - x_0}{h}$
- $\Delta y_0, \Delta^2 y_0, \ldots$ are forward differences

**Algorithm**

1. **Input:** Read $n$ data points $(x_i, y_i)$, the value at which derivative is needed, and the order of derivative.
2. **Difference Table:** Construct the forward difference table.
3. **Calculate $u$:** Compute $u = \frac{x - x_0}{h}$ where $h = x_1 - x_0$.
4. **Apply Formula:** Use the appropriate Newton's forward differentiation formula based on the order.
5. **Output:** Display the derivative value.

**Pseudocode**

```text
START

Read n
Read x[i], y[i][0]

FOR i = 1 to n-1
    FOR j = 0 to n-i-1
        y[j][i] = y[j+1][i-1] - y[j][i-1]
    END FOR
END FOR

Read value, order

h = x[1] - x[0]
u = (value - x[0]) / h

derivative = 0
P = 1

FOR k = 0 to n-1
    IF k > 0
        P = P × (u - (k-1))
    END IF

    IF k ≥ order
        Differentiate P (order times)
        derivative += (y[0][k] / k!) × P(u)
    END IF
END FOR

derivative = derivative / h^order

Print derivative

END
```

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/E.%20Solution%20of%20Differential%20Equations/Newton%20Forward%20Differentiation)

Includes:

- C++ source code
- Input file
- Output file

**Further Study**

- [Numerical Differentiation - Wikipedia](https://en.wikipedia.org/wiki/Numerical_differentiation)
- [Newton's Forward Difference Formula - MathWorld](https://mathworld.wolfram.com/ForwardDifference.html)
- Numerical Methods for Engineers - Chapra & Canale

<a id="runge-kutta"></a>

### 2. Runge Kutta Method

**Theory**
The Runge-Kutta method (specifically the fourth-order RK4) is a widely used technique for the approximate solution of ordinary differential equations (ODEs).
$$y_{n+1} = y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

**Algorithm**

1.  **Define function**: $f(x, y)$.
2.  **Input**: Initial $x_0, y_0$, target $x_n$, and step size $h$.
3.  **Iterate**: Until current $x$ reaches target $x_n$, calculate slopes $K$ and update $y$.
4.  **Output**: Final value of $y$.

**Pseudocode**

```cpp
function solve_rk4(x0, y0, xn, h):
    n_steps = (xn - x0) / h
    y = y0
    x = x0

    for i from 0 to n_steps:
        if x >= xn: break

        k1 = h * f(x, y)
        k2 = h * f(x + h/2, y + k1/2)
        k3 = h * f(x + h/2, y + k2/2)
        k4 = h * f(x + h, y + k3)

        k_avg = (k1 + 2*k2 + 2*k3 + k4) / 6.0

        y = y + k_avg
        x = x + h

    return y
```

**Implementation**

👉 [View Code & Files](https://github.com/suaib022/nm-demo/tree/main/E.%20Solution%20of%20Differential%20Equations/Runge%20Kutta)

Includes:

- C++ source code
- Input file
- Output file

**Further Study**

- [Runge-Kutta Methods - Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
- [Runge-Kutta Method - GeeksforGeeks](https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/)
- Numerical Methods for Engineers - Chapra & Canale

---

## F. Numerical Integration

<a id="simpson-1-by-3"></a>

### 1. Simpson’s 1/3 Rule

<a id="simpson-3-by-8"></a>

### 2. Simpson’s 3/8 Rule

## G. Curve Fitting Methods

<a id="linear"></a>

### 1. Linear Equation

<a id="polynomial"></a>

### 2. Polynomial Equation

<a id="transcendental"></a>

### 3. Transcendental Equation

---

## Environment Setup

### Requirements

- C++ Compiler (GCC / MinGW / Clang)
- C++11 or later

### Compile & Run

```bash
g++ code.cpp -o run
./run
```

---

## Sample Run

### Input

```text
Enter number of equations: 3
Enter coefficients...
```

### Output

```text
Solution:
x = 1.0
y = 2.0
z = 3.0
```

---

## Contribution

### Our Team Contribution

| Full Name              | GitHub Username                                                     | Roll Number |
| ---------------------- | ------------------------------------------------------------------- | ----------- |
| MD SUAIB AHMED SAFI    | [suaib022](https://github.com/suaib022)                             | 2207115     |
| ASHRAFUR RAHMAN NIHAD  | [ARN101](https://github.com/ARN101)                                 | 2207116     |
| DADHICHI SAREKR SHAYON | [Dadhichi-Sarker-Shayon](https://github.com/Dadhichi-Sarker-Shayon) | 2207118     |

- **Member 1:** Method implementation & testing
- **Member 2:** Mathematical theory & PDFs
- **Member 3:** Repository structure & documentation

---

## Future Contributions

We welcome contributions! 🙌

**How to contribute:**

1. Fork the repository
2. Create a new branch
3. Add method / optimize code / improve documentation
4. Submit a pull request

---

## References

1. Chapra, S. C., & Canale, R. P. _Numerical Methods for Engineers_
2. Burden & Faires, _Numerical Analysis_
3. NPTEL Lecture Notes on Numerical Methods
4. _Numerical Methods for Engineers_ — Steven Chapra
5. MIT OpenCourseWare (Numerical Analysis)
6. NPTEL Numerical Methods
