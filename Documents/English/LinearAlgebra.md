# LinearAlgebra

## Introduction

Integrates content related to linear algebra. Currently includes: a matrix class template.

## Dependencies

Please ensure that your C++ version is at least C++14.

## Matrix Class Template

### 1. Location

This type is located in the `QMath` namespace within the `LinearAlgebra.hpp` file. You can access it via `QMath::Matrix<...>`.

### 2. Template Parameters

Its template parameter definition is as follows:

```cpp
template <typename Type> class Matrix;
```

### 3. Constructors

The main constructors are as follows:

```cpp
// Func1
constexpr Matrix(): n(0), m(0), data(nullptr) {}

// Func2
constexpr Matrix(size_t n, size_t m, const Type& x = Type());

// Func3
Matrix(const std::vector<std::vector<Type>>& x);
```

- **Func1** constructs an empty matrix.
- **Func2** allows you to specify the number of rows, the number of columns, and an initial value.
- **Func3** constructs a matrix from a `std::vector<std::vector<...>>`. Please ensure that every inner `vector` in the outer `vector` has the same length.

There are also some other common constructors that are not worth listing here.

### 4. Member Functions and Overloaded Operators

| Member Function / Operation | Description | Time Complexity (For member functions or functions involving only one matrix, $N$ and $M$ denote the number of rows and columns. For operators involving two matrices, $N$ and $M$ denote the number of rows and columns of the left matrix, while $P$ and $Q$ denote the number of rows and columns of the right matrix) | Example |
|:-:|:-:|:-:|:-:|
| `Matrix& operator=(const Matrix& o)` | Copies from `o` and returns an lvalue reference. | $O(NM+PQ)$ | `a = b` |
| `Matrix& operator=(const std::vector<std::vector<Type>>& o)` | Copies from a `std::vector` and returns an lvalue reference, similar to **Func3** of the constructor. | $O(NM+PQ)$ | `a = b` |
| `Matrix operator+(const Matrix& o)` | Matrix addition. Please ensure $N=P$ and $M=Q$. | $O(NM)$ | `a + b` |
| `Matrix operator+(const Type &o)` | Adds the scalar `o` to every element of the matrix and returns the result. | $O(NM)$ | `a + o` |
| `Matrix operator-(const Matrix& o)` | Matrix subtraction. Please ensure $N=P$ and $M=Q$. | $O(NM)$ | `a - b` |
| `Matrix operator-(const Type &o)` | Subtracts the scalar `o` from every element of the matrix and returns the result. | $O(NM)$ | `a - o` |
| `Matrix& operator+=(const Matrix& o)` | $this \gets this + o$. Please ensure $N=P$ and $M=Q$. | $O(NM)$ | `a += b` |
| `Matrix& operator+=(const Type &o)` | Adds the scalar `o` to every element of the matrix and returns an lvalue reference. | $O(NM)$ | `a += o` |
| `Matrix& operator-=(const Matrix& o)` | $this \gets this - o$. Please ensure $N=P$ and $M=Q$. | $O(NM)$ | `a -= b` |
| `Matrix& operator-=(const Type &o)` | Subtracts the scalar `o` from every element of the matrix and returns an lvalue reference. | $O(NM)$ | `a -= o` |
| `Matrix operator*(const Matrix& o)` | Matrix multiplication. Please ensure $P=M$. | $O(NMQ)$ | `a * b` |
| `Matrix operator*(const Type& o)` | Multiplies every element of the matrix by the scalar `o` and returns the result. | $O(NM)$ | `a * o` |
| `Matrix& operator*=(const Matrix& o)` | $this \gets this \times o$. Please ensure $P=M$. | $O(NMQ)$ | `a *= b` |
| `Matrix& operator*=(const Type &o)` | Multiplies every element of the matrix by the scalar `o` and returns an lvalue reference. | $O(NM)$ | `a *= o` |
| `Matrix operator%(const Matrix& o)` | Multiplies corresponding elements of two matrices with the same shape and returns the result. Please ensure $N=P$ and $M=Q$. | $O(NM)$ | `a % b` |
| `Matrix& operator%=(const Matrix& o)` | Multiplies corresponding elements of two matrices with the same shape, stores the result in the left-hand operand, and returns an lvalue reference. Please ensure $N=P$ and $M=Q$. | $O(NM)$ | `a %= b` |
| `bool operator==(const Matrix& o)` | Returns `true` if the two matrices have identical shape, size, and content; otherwise returns `false`. | $O(NM+PQ)$ | `a == b` |
| `bool operator!=(const Matrix& o)` | Returns `false` if the two matrices have identical shape, size, and content; otherwise returns `true`. | $O(NM+PQ)$ | `a != b` |
| `(const) Type& operator()(size_t row, size_t col)` | Returns the value at row `row` and column `col` of the matrix (by reference). | $O(1)$ | `a(0, 1)` |
| `(const) Type* operator[](size_t row)` | Returns a pointer to the beginning of row `row` of the matrix. Each row of the matrix is stored contiguously in memory. | $O(1)$ | `a[0]` |
| `size_t N()` | Returns the number of rows of the matrix. | $O(1)$ | `a.N()` |
| `size_t M()` | Returns the number of columns of the matrix. | $O(1)$ | `a.M()` |
| `Matrix transpose()` | Returns the transpose of the matrix. | $O(NM)$ | `a.transpose()` |
| `Matrix& transposeSelf()` | Transposes the matrix in place and returns an lvalue reference to itself. | $O(NM)$ | `a.transposeSelf()` |
| `Matrix& resize(size_t n, size_t m, Type x = Type())` | Changes the shape of the matrix. Excess elements are cleared, and missing elements are filled with `x`. | $O(NM)$ | `a.resize(128, 256, 2)` |
| `Matrix& assgin(size_t n, size_t m, Type x = Type())` | Clears the current matrix and constructs an $n \times m$ matrix with all elements initialized to `x`. | $O(NM)$ | `a.assgin(128, 256, 2)` |
| `Matrix applyFunction(Type (*func)(Type))` / `Matrix applyFunction(Type (*func)(const Type&))` | Applies `func` to every element $x$ of the matrix and returns the resulting new matrix. | $O(NM)$ | `a.applyFunction(std::sqrt)` |
| `Matrix& applyFunctionSelf(Type (*func)(Type))` / `Matrix& applyFunctionSelf(Type (*func)(const Type&))` | Applies `func` to every element $x$ of the matrix in place and returns an lvalue reference to itself. | $O(NM)$ | `a.applyFunctionSelf(std::sqrt)` |