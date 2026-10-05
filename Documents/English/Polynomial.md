# Polynomial

## Introduction

Integrates polynomial-related functionalities. Currently includes: a polynomial class template.

### Overview of Advantages

1. Excellent Generality

    Through simple template parameter modifications, it can switch seamlessly between FFT and NTT modes, and can even be used together with a fully implemented arbitrary-precision class.

2. Excellent Speed

    While maintaining generic support, it can enter the [first page of the optimal solutions for the Luogu template problem](https://www.luogu.com.cn/record/list?pid=P3803&orderBy=1&status=&page=1), whether in FFT mode or NTT mode. The overall execution time displayed by Luogu can stably remain below 490ms. The visible submissions that are faster than it all use special optimizations targeting the modulus or the `int` data type, and therefore do not support generic types.

3. Refined Encapsulation


    It fully embraces object-oriented programming principles. This class inherits from `std::vector` and is encapsulated in the style of the STL, allowing seamless integration with the STL. This template does not expose any member functions that are useless to users.

4. Powerful Functionality

    It also implements polynomial `exp`, `ln`, `sqrt`, and other operations. It additionally uses FWT and FMT to implement OR, AND, and XOR convolutions.

6. Hybrid Use of Recursive and Iterative Implementations

    The technical implementation is relatively complex. See “Technical Implementation Sharing” at the end for details.

## Dependencies

Please ensure that the C++ version is greater than or equal to C++17.

## Polynomial Class Template

### I. Location

This type is located in the `QMath` namespace of the `Polynomial.hpp` file. You can invoke it through `QMath:Polynomial<...>`.

### II. Template Parameters

Here is its template parameter definition:

```cpp
template <
    typename T = double, // use double as the default data storage type
    typename TDFT = std::complex<double>, // use std::complex<double> as the default transform type
    auto UR = expn // unit root calculation function, accepting a size_t input value n and returning the n-th root of unity of type TDFT. Due to the special design of the Polynomial class, there is no need to implement a function for calculating the inverse of the unit root.
    auto T2TDFT = T2TFFT, // provide a conversion function from T (the storage type) to TDFT (the transform type); 2 is a homophone of "to"
    auto TDFT2T = TFFT2T, // same as above
> class Polynomial : public std::vector<T>; // inherits from std::vector
```

Given the above template parameter functions, it determines whether the transform is FFT or NTT, thereby supporting multiple types of transforms within a single template. ****Note: When the `T` type and `TDFT` type are the same, calls to `T2TDFT & TDFT2T` are automatically omitted. In this case, they can be passed arbitrarily or simply omitted altogether.****

The `QMath` namespace in the `Polynomial` file provides a set of functions required by the above template parameters and uses them as default parameters. They are respectively: `QMath::T2TFFT`, `QMath::TFFT2T`, and `QMath::expn`.

Their implementations are:

```cpp
namespace QMath {
    namespace detail {
        constexpr double PI2 = 6.283185307179586476925286766559005768394338798750211641949889;
    }
    std::complex<double> expn(size_t n) { return std::complex<double>(std::cos(detail::PI2 / n), std::sin(detail::PI2 / n)); }

    std::complex<double> T2TFFT(double x) { return std::complex<double>(x, 0); }

    double TFFT2T(const std::complex<double>& x) { return x.real(); }
}
```

Very concise; these are simply the implementations required by ordinary FFT.

Of course, this class template also supports NTT, but it requires an automatic modular arithmetic type and a corresponding set of implementations. It is recommended to use `QMath::Modular` from this repository. Click [here](./Modular.md) to read about `QMath::Modular`.

****Note: If you use another automatic modular arithmetic class template when testing the performance of the Polynomial class template, you may not achieve the expected performance.****

Below is an implementation of a `Polynomial` type definition based on the modulus `998244353`:

```cpp
#include "Modular"
using Mint = QMath::Modular<>; // define an automatic modular arithmetic type and use the default template parameters

inline Mint T2TNTT(Mint x) { return x; }
inline Mint TNTT2T(Mint x) { return x; }
inline Mint NTT(size_t x) { return Mint(3).fPow(998244352 / x); }

typedef QMath::Polynomial<Mint, Mint, NTT, T2NTT, NTT2T> Poly;
```

### III. Constructors

The main constructor is of the following form:

```cpp
Polynomial(size_t n = 0, const T& init = T()): std::vector<T>(n, init);
```

It accepts two parameters, representing the initial length of the polynomial and the initial coefficient of each term, respectively. The default initial length is 0, and the coefficient is `T()`.

There are also some other common constructors that do not need to be shown here.

### IV. Member Functions and Overloaded Operators

|Member Function / Operation|Purpose|Time Complexity ($n$ denotes the polynomial length)|Example|
|:-:|:-:|:-:|:-:|
|Everything inherited from `std::vector`|Omitted|Omitted|Omitted|
|`void clearCache()`|This template performs $O(n)$ preprocessing of the roots of unity required for its operations and stores them in the form of `thread_local static`, shared among objects of the same type. You can use this member function to clear the stored roots of unity.|$O(n)$|`a.clearCache()`|
|`Polynomial& operator=(const Polynomial& o)`|Copies and returns an lvalue reference|$O(n)$|`a = b`|
|`Polynomial derivative()`|Returns its own derivative polynomial|$O(n)$|`a.derivative()`|
|`Polynomial& derivativeSelf()`|Replaces itself with its own derivative polynomial and returns its own lvalue reference|$O(n)$|`a.derivativeSelf()`|
|`Polynomial integral()`|Returns its own integral polynomial|$O(n)$|`a.integral()`|
|`Polynomial& integralSelf()`|Replaces itself with its own integral polynomial and returns its own lvalue reference (with the constant term set to zero)|$O(n)$|`a.integralSelf()`|
|`T calc(const T &x)`|Evaluates itself at $x$ and returns the result|$O(n)$|`a.calc(1)`|
|`T calcDerivative(const T &x)`|Evaluates its derivative polynomial at $x$ and returns the result|$O(n)$|`a.calcDerivative(2)`|
|`T calcIntegral(const T &x)`|Evaluates its integral polynomial at $x$ and returns the result|$O(n)$|`a.calcIntegral(3)`|
|`Polynomial operator+(const Polynomial& o) const`|Polynomial addition; automatically extends to the longer of the two lengths and returns the result|$O(n)$|`a = a + b`|
|`Polynomial& operator+=(const Polynomial& o)`|Polynomial addition assignment; automatically extends to the longer of the two lengths and returns an lvalue reference|$O(n)$|`a += b`|
|`Polynomial operator-(const Polynomial& o) const`|Polynomial subtraction; automatically extends to the longer of the two lengths and returns the result|$O(n)$|`a = a - b`|
|`Polynomial& operator-=(const Polynomial& o)`|Polynomial subtraction assignment; automatically extends to the longer of the two lengths and returns an lvalue reference|$O(n)$|`a -= b`|
|`Polynomial operator*(const T& o) const`|Multiplies the polynomial by a scalar and returns the result|$O(n)$|`a = a * 4`|
|`Polynomial& operator*=(const T& o)`|Multiplies the polynomial by a scalar and returns an lvalue reference|$O(n)$|`a *= 4`|
|`Polynomial operator/(const T& o) const`|Divides the polynomial by a scalar and returns the result|$O(n)$|`a = a / 5`|
|`Polynomial& operator/=(const T& o)`|Divides the polynomial by a scalar and returns an lvalue reference|$O(n)$|`a /= 5`|
|`Polynomial operator*(const Polynomial& o) const`|Polynomial multiplication; returns the result|$O(nlog{n})$|`a = a * b`|
|`Polynomial& operator*=(const Polynomial& o)`|Polynomial multiplication assignment; returns an lvalue reference|$O(nlog{n})$|`a *= b`|
|`Polynomial& operator%=(size_t n)`|Similar to `.resize(n)`, but more intuitive by using the modular definition of a polynomial (regardless of whether the original length is greater than or less than $n$, it is uniformly extended to $n$), and returns its own lvalue reference|$O(n)$|`a %= 6`|
|`Polynomial operator%(size_t n) const`|Does not modify itself and returns the result of `.resize(n)`|$O(n)$|`a = a % 7`|
|`Polynomial inv() const`|Returns its multiplicative inverse|$O(nlog{n})$|`a.inv()`|
|`Polynomial& invSelf()`|Replaces itself with its multiplicative inverse and returns its own lvalue reference|$O(nlog{n})$|`a.invSelf()`|
|`Polynomial operator/(const Polynomial& o) const`|Uses the polynomial multiplicative inverse to implement division and returns the result|$O(nlog{n})$|`a = a / b`|
|`Polynomial& operator/=(const Polynomial& o)`|Uses the polynomial multiplicative inverse to implement division and returns its own lvalue reference|$O(nlog{n})$|`a /= b`|
|`template<auto LN = __return0> Polynomial ln() const`|Returns the result of taking $ln$ of itself (if it is not guaranteed that the constant term is $1$, a function for calculating $ln$ of type `T` must be provided; if none is provided, a function that always returns `T(0)` is supplied by default)|$O(nlog{n})$|`a.ln()`|
|`template<auto LN = __return0> Polynomial& lnSelf()`|Replaces itself with the result of taking $ln$ of itself and returns its own lvalue reference (if it is not guaranteed that the constant term is $1$, a function for calculating $ln$ of type `T` must be provided; if none is provided, a function that always returns `T(0)` is supplied by default)|$O(nlog{n})$|`a.lnSelf()`|
|`template<auto EXP = __return1> Polynomial exp() const`|Returns the result of taking $exp$ of itself (if it is not guaranteed that the constant term is $0$, a function for calculating $exp$ of type `T` must be provided; if none is provided, a function that always returns `T(1)` is supplied by default)|$O(nlog{n})$|`a.exp<exp>()`|
|`template<auto EXP = __return1> Polynomial& expSelf()`|Replaces itself with the result of taking $exp$ of itself and returns its own lvalue reference (if it is not guaranteed that the constant term is $0$, a function for calculating $exp$ of type `T` must be provided; if none is provided, a function that always returns `T(1)` is supplied by default)|$O(nlog{n})$|`a.expSelf<exp>()`|
|`template<auto SQRT = __return1> Polynomial sqrt()`|Returns its square root (if it is not guaranteed that the constant term is $1$, a function for calculating the square root of type `T` must be provided; if none is provided, a function that always returns `T(1)` is supplied by default)|$O(nlog{n})$|`a.sqrt()`|
|`template<auto SQRT = __return1> Polynomial& sqrtSelf()`|Replaces itself with its square root and returns its own lvalue reference (if it is not guaranteed that the constant term is $1$, a function for calculating the square root of type `T` must be provided; if none is provided, a function that always returns `T(1)` is supplied by default)|$O(nlog{n})$|`a.sqrtSelf()`|
|`template<typename U> Polynomial pow(U n, U m = U(-1))`|Computes its own $n$-th power (the constant term does not need to be guaranteed to be $1$, but must not be $0$. If the exponent needs to be reduced modulo a modulus to reduce its range, then $n$ should be passed as the exponent modulo $MOD$, and $m$ should be passed as the exponent modulo $varphi(MOD)$; otherwise, you may ignore $m$)|$O(nlog{n})$|`a.pow(10)`|
|`template<typename U> Polynomial& powSelf(U n, U m = U(-1))`|Replaces itself with its own $n$-th power and returns its own lvalue reference (the constant term does not need to be guaranteed to be $1$, but must not be $0$. If the exponent needs to be reduced modulo a modulus to reduce its range, then $n$ should be passed as the exponent modulo $MOD$, and $m$ should be passed as the exponent modulo $varphi(MOD)$; otherwise, you may ignore $m$)|$O(nlog{n})$|`a.powSelf(9982, 44353)`|
|`Polynomial operator\|(const Polynomial& o) const`|Computes OR convolution and returns the result, automatically extending to the smallest power of two greater than or equal to the lengths of both polynomials|$O(nlog{n})$|`a = a | b`|
|`Polynomial& operator\|=(const Polynomial& o)`|Computes OR convolution and returns an lvalue reference, automatically extending to the smallest power of two greater than or equal to the lengths of both polynomials|$O(nlog{n})$|`a |= b`|
|`Polynomial operator&(const Polynomial& o) const`|Computes AND convolution and returns the result, automatically extending to the smallest power of two greater than or equal to the lengths of both polynomials|$O(nlog{n})$|`a = a & b`|
|`Polynomial& operator&=(const Polynomial& o)`|Computes AND convolution and returns an lvalue reference, automatically extending to the smallest power of two greater than or equal to the lengths of both polynomials|$O(nlog{n})$|`a &= b`|
|`Polynomial operator^(const Polynomial& o)`|Computes XOR convolution and returns the result, automatically extending to the smallest power of two greater than or equal to the lengths of both polynomials|$O(nlog{n})$|`a = a ^ b`|
|`Polynomial& operator^=(const Polynomial& o)`|Computes XOR convolution and returns an lvalue reference, automatically extending to the smallest power of two greater than or equal to the lengths of both polynomials|$O(nlog{n})$|`a ^= b`|
|`friend std::ostream& operator<<(std::ostream& os, const Polynomial& p)`|Output stream overload. Output format: `f(x) = ax ^ 0 + bx ^ 1 + cx ^ 2......`|$O(n)$|`std::cout << a`|

## Technical Implementation Sharing

After reading "[Are Recursive Algorithms Really Slower Than Iterative Ones?](https://www.luogu.com.cn/article/tfoqhji5)", I found that recursive algorithms are actually more) Cache-friendly; the main issue is the overhead of recursion. The solution proposed in this article is to use template recursion so that the compiler expands the recursive functions at compile time. The effect is indeed significant. However, there are two problems:

* After adding this optimization, the measured performance on the Luogu template problem was: the overall FFT execution time decreased from around 580 ms to around 520 ms, while the overall NTT execution time deteriorated from around 480 ms to around 510 ms.

    Analysis of the cause: The reason iterative implementations are Cache-unfriendly is that each time the level changes, the polynomial array has to be traversed again, causing the tail of the polynomial to be evicted from the cache and the beginning to be loaded into it. The recursive implementation avoids this problem through its excellent recursive ordering. NTT uses `int` for storage, and `int` occupies little space. Under the data range of the Luogu template problem, almost the entire polynomial can fit into the L2 / L1 cache, so Cache optimization has little effect. FFT, however, uses `std::complex<double>` for storage, which occupies 4 times as much space as `int`, so the cache cannot hold it all, making Cache optimization much more effective.

    Solution: Use iteration for short sequences and recursion for long sequences. But is this the only way? Is that where we stop? Although the transform implementations of NTT / FFT can be divided into iterative and recursive forms, what they fundamentally do is the same; they only differ in their access order. Why can't iterative and recursive methods be used simultaneously on the same sequence? Therefore, I used a hybrid of iterative and recursive methods: first use the recursive implementation. Each time a recursion is entered, read the length of the current recursive interval. If the length is large, continue recursively. If the length is smaller than a threshold (4096 in this code), so small that even an iterative implementation can fit entirely into the cache, stop the recursion and implement everything remaining iteratively. The recursive structure may then look like this:

    ```
       ----------------
     --------    --------
    ----  ----  ----  ----
    ++++  ++++  ++++  ++++
    ++++  ++++  ++++  ++++
    ++++  ++++  ++++  ++++
    ```

    The `-` sections are intervals with large lengths and use recursion. The `+` sections are very small intervals and are completed directly using iteration. As can be seen, the top of the recursion tree is computed recursively, while the bottom is computed iteratively. Although the traversal methods are different, they perform the same operations, so correctness is unaffected. The [step-by-step tutorial on obtaining the optimal solution for the polynomial toolkit]([https://www.luogu.com.cn/article/k9j38kqv](https://www.luogu.com.cn/article/k9j38kqv), also](https://www.luogu.com.cn/article/k9j38kqv)貌似也提到了这个东西，但比较笼统，不知道意思是长的原序列全迭代，短的原序列全递归，还是把递归树拆开，同长度原序列递归迭代混用。) seems to mention this as well, but rather vaguely. I do not know whether it means that long original sequences are entirely iterative and short original sequences are entirely recursive, or that the recursion tree is split so that iterative and recursive methods are mixed for original sequences of the same length.

* The supported polynomial length range is limited because expanding recursive functions for long polynomials during compilation causes a dramatic increase in the compiled file size. Therefore, extremely long polynomials cannot be supported. Although this is acceptable under the modulus 998244353, it violates the principle of generality.


    Solution: If the polynomial is found to be too long, abandon the optimization and directly use iterative computation for the highest and longest few levels of the recursive structure. Once the problem size is reduced to a range that recursion can handle, recursion begins. The maximum recursive length supported by this code is $2^{24}$, and the Luogu template problem does not exceed this limit.

### Optimization Results

|Algorithm|NTT|FFT|
|:-:|:-:|:-:|
|Pure iterative version|https://www.luogu.com.cn/record/264397472 478ms|https://www.luogu.com.cn/record/264424326 572ms|
|Hybrid version|https://www.luogu.com.cn/record/264423088 477ms|https://www.luogu.com.cn/record/264422355 473ms|

As can be seen, while not affecting NTT speed, FFT has been optimized to a performance level comparable to NTT.

### Some Other Minor Optimizations Used

* If `TDFT` = `std::complex<T>`, the implementation automatically uses the three-for-two optimization of the Fast Fourier Transform. If `std::complex<>` is not used, this optimization is unavailable. It is strongly recommended to use `std::complex<>`!!!

* If `T` = `TDFT`, the calls to `T2TDFT()` and `TDFT2T()` are omitted. These two functions may be passed or omitted arbitrarily; it is recommended not to pass them.

* If the multiplication operation is of the form `a *= a` or `a = a * a`, only one transform is performed during the forward transform.
