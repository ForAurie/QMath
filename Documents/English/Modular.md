# Modular

## Introduction

Integrates automatic modular arithmetic types. Currently available: automatic modular arithmetic class template.

## Dependencies

Please ensure your C++ version is at least C++14.

## Modular Class Template

### 1. Location

This type resides in the `QMath` namespace within the `Modular.hpp` file. You can access it via `QMath::Modular<...>`.

### 2. Template Parameters

Here is its template parameter definition:

```cpp
template <
    typename T = int, // The type used to store data. Be sure to use a signed type and ensure its value range is no smaller than (-modulus * 2 + 1, modulus * 2 - 1).
    const T MOD = 998244353, // Default modulus
    typename MCT = unsigned long long // During multiplication, values will be temporarily cast to this type. Please ensure that a cast exists between T and MCT, and that its value range is no smaller than [0, (modulus - 1) ^ 2]. An unsigned type is recommended.
> class Modular;
```

### 3. Constructor

The primary constructor has the following form:
```cpp
constexpr Modular(const T& v = T());
```

You may pass one argument representing the initial value of this type, which defaults to `T()`.

Additionally, there are some commonly used constructors that are not necessary to display here.

### 4. Member Functions and Overloaded Operators

| Member Function / Operation | Description | Time Complexity ($M$ denotes the modulus) | Example |
|:-:|:-:|:-:|:-:|
| `Modular& operator=(const Modular& other)` | Copies and returns an lvalue reference | $O(1)$ | `a = b` |
| `Modular& operator=(const T& v)` | Constructs, assigns, and then returns an lvalue reference | $O(1)$ | `a = 8` |
| `const T& getVal()` | Returns the corresponding value of type `T` | $O(1)$ | `a.getVal()` |
| `Modular& setVal(const T& v)` | Constructs, assigns, and then returns an lvalue reference | $O(1)$ | `a.setVal(8)` |
| `Modular operator-()` | Returns the negation of itself | $O(1)$ | `-a` |
| `friend Modular operator+(const Modular& a, const Modular& b)` | Computes $a+b$, automatically applies modular reduction, and returns the result | $O(1)$ | `a + b` |
| `friend Modular operator-(const Modular& a, const Modular& b)` | Computes $a-b$, automatically applies modular reduction, and returns the result | $O(1)$ | `a - b` |
| `friend Modular operator*(const Modular& a, const Modular& b)` | Computes $a\times b$, automatically applies modular reduction, and returns the result | $O(1)$ | `a * b` |
| `friend Modular operator/(const Modular& a, const Modular& b)` | Computes $a \div b$, automatically applies modular reduction, and returns the result (when using this feature, please ensure the modulus is prime and the divisor is non-zero) | $O(\log{M})$ | `a / b` |
| `Modular& operator+=(const Modular& other)` | $a\gets a+b$, automatically applies modular reduction, and returns an lvalue reference | $O(1)$ | `a += b` |
| `Modular& operator-=(const Modular& other)` | $a\gets a-b$, automatically applies modular reduction, and returns an lvalue reference | $O(1)$ | `a -= b` |
| `Modular& operator*=(const Modular& other)` | $a\gets a\times b$, automatically applies modular reduction, and returns an lvalue reference | $O(1)$ | `a *= b` |
| `Modular& operator/=(const Modular& other)` | $a\gets a\div b$, automatically applies modular reduction, and returns an lvalue reference (when using this feature, please ensure the modulus is prime and the divisor is non-zero) | $O(\log{M})$ | `a /= b` |
| `template <typename U> Modular fPow(U exp)` | Uses fast exponentiation to compute $a^U$, automatically applies modular reduction, and returns the result | $O(\log{U})$ | `a.fPow(8)` |
| `template <typename U> Modular& fPowSelf(U exp)` | Uses fast exponentiation to compute $a\gets a^U$, automatically applies modular reduction, and returns a reference to itself | $O(\log{U})$ | `a.fPowSelf(8)` |
| `Modular inv()` | Computes the modular multiplicative inverse of itself and returns the result (when using this feature, please ensure the modulus is prime and the operand is non-zero) | $O(\log{M})$ | `a.inv()` |
| `Modular& invSelf()` | Sets itself to its modular multiplicative inverse and returns an lvalue reference to itself (when using this feature, please ensure the modulus is prime and the operand is non-zero) | $O(\log{M})$ | `a.invSelf()` |
| `friend bool operator==(const Modular& a, const Modular& b)` | Returns true if both values are equal, otherwise returns false | $O(1)$ | `a==b` |
| `friend bool operator!=(const Modular& a, const Modular& b)` | Returns true if both values are not equal, otherwise returns false | $O(1)$ | `a!=b` |
| `friend std::ostream& operator<<(std::ostream& os, const Modular& m)` | Outputs it as a value of type `T` | $O(1)$ | `std::cout << a` |
| `std::istream& operator>>(std::istream& is, Modular& m)` | Reads input as a value of type `T` | $O(1)$ | `std::cin >> a` |