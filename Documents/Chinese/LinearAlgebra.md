# LinearAlgerbra

## 简介

集成线性代数相关内容，目前有：矩阵类模板。

## 依赖

请确保 C++ 版本大于等于 C++14。

## Matrix 类模板

### 一、位置

该类型位于 `LinearAlgebra` 文件的 `QMath` 命名空间中，你可以通过 `QMath:Matrix<...>` 调用。

### 二、模板参数

这是它的模板参数定义：

```cpp
template <typename Type> class Matrix;
```

### 三、构造函数

主要的构造函数形如：
```cpp
// Func1
constexpr Matrix(): n(0), m(0), data(nullptr) {}

// Func2
constexpr Matrix(size_t n, size_t m, const Type& x = Type());

// Func3
Matrix(const std::vector<std::vector<Type>>& x);
```

Func1 会构造一个空矩阵。

Func2 支持传入矩阵的高、宽和初始值。

Func3 支持从 `std::vector<std::vector<...>>` 构造一个矩阵，请确保第一层 `vector` 里的每一个 `vector` 长度均一样。

除此之外还有一些没必要展示的常用构造函数。

### 四、成员函数和重载运算符

|成员函数 / 运算|作用|时间复杂度（对于成员函数或只涉及一个矩阵的函数，$N,M$ 表示矩阵的行数和列数。对于涉及两个矩阵的运算符 $N,M$ 表示左矩阵的行数和列数，$P,Q$ 表示右矩阵的行数和列数）|示例|
|:-:|:-:|:-:|:-:|
|`Matrix& operator=(const Matrix& o)`|拷贝并返回左值引用|$O(NM+PQ)$|`a = b`|
|`Matrix& operator=(const std::vector<std::vector<Type>>& o)`|从 `std::vector` 拷贝并返回左值引用，和构造函数的 Func3相同|$O(NM+PQ)$|`a = b`|
|`Matrix operator+(const Matrix& o)`|矩阵加法，请保证 $N=P, M=Q$|$O(NM)$|`a + b`|
|`Matrix operator+(const Type &o)`|矩阵的每一个位置都加上数字 o，返回结果|$O(NM)$|`a + o`|
|`Matrix operator-(const Matrix& o)`|矩阵减法，请保证 $N=P, M=Q$|$O(NM)$|`a - b`|
|`Matrix operator-(const Type &o)`|矩阵的每一个位置都减去数字 o，返回结果|$O(NM)$|`a - o`|
|`Matrix& operator+=(const Matrix& o)`|$this\gets this+o$，请保证 $N=P, M=Q$|$O(NM)$|`a += b`|
|`Matrix& operator+=(const Type &o)`|给自己的每一个位置都加上数字 o，返回左值引用|$O(NM)$|`a += o`|
|`Matrix& operator-=(const Matrix& o)`|$this\gets this-o$，请保证 $N=P, M=Q$|$O(NM)$|`a -= b`|
|`Matrix& operator-=(const Type &o)`|给自己的每一个位置都减去数字 o，返回左值引用|$O(NM)$|`a -= o`|
|`Matrix operator*(const Matrix& o)`|矩阵乘法，请保证 $P=M$|$O(NMQ)$|`a * b`|
|`Matrix operator*(const Type& o)`|矩阵中每个位置乘数字 o，返回结果|$O(NM)$|`a * o`|
|`Matrix& operator*=(const Matrix& o)`|$this\gets this\times o$，请保证 $P=M$|$O(NMQ)$|`a * b`|
|`Matrix& operator*=(const Type &o)`|给自己的每一个位置都乘数字 o，返回左值引用|$O(NM)$|`a *= o`|
|`Matrix operator%(const Matrix& o)`|将两个大小形状相同的矩阵中的每个数对应相乘，返回结果，请保证 $N=P,M=Q$|$O(NM)$|`a % b`|
|`Matrix& operator%=(const Matrix& o)`|将两个大小形状相同的矩阵中的每个数对应相乘，并将结果存在运算符左侧的对象，返回左直引用，请保证 $N=P,M=Q$|$O(NM)$|`a %= b`|
|`bool operator==(const Matrix& o)`|两个矩阵形状、大小、内容完全相同则返回真否则返回假|$O(NM+PQ)$|`a == b`|
|`bool operator!=(const Matrix& o)`|两个矩阵形状、大小、内容完全相同则返回假否则返回真|$O(NM+PQ)$|`a != b`|
|`(const) Type& operator()(size_t row, size_t col)`|返回矩阵第 $row$ 行第 $col$ 列的值（带引用）。|$O(1)$|`a(0, 1)`|
|`(const) Type* operator[](size_t row)`|返回矩阵第 $row$ 行的起始指针，矩阵的每一行在内存中连续存储|$O(1)$|`a[0]`|
|`size_t N()`|返回矩阵的行数|`O(1)`|`a.N()`|
|`size_t M()`|返回矩阵的列数|`O(1)`|`a.M()`|
|`Matrix transpose()`|返回自己转置的结果|$O(NM)$|`a.transpose()`|
|`Matrix& transposeSelf()`|令自己转置并返回自己的左值引用|$O(NM)$|`a.transposeSelf()`|
|`Matrix& resize(size_t n, size_t m, Type x = Type())`|改变矩阵的形状，多余地方清空，不足地方补 $x$|$O(NM)$|`a.resize(128, 256, 2)`|
|`Matrix& assgin(size_t n, size_t m, Type x = Type())`|将自己清空，并构造一个 $n$ 行 $m$ 列初始值均为 $x$ 的矩阵|$O(NM)$|`a.assgin(128, 256, 2)`|
|`Matrix applyFunction(Type (*func)(Type)) / Matrix applyFunction(Type (*func)(const Type&))`|矩阵中的每一个数 $x$ 变为 $func(x)$，然后返回新的结果|$O(NM)$|`a.applyFunction(std::sqrt)`|
|`Matrix& applyFunctionSelf(Type (*func)(Type)) / Matrix& applyFunctionSelf(Type (*func)(const Type&))`|令矩阵中的每一个数 $x$ 变为 $func(x)$，然后返回自己的左值引用|$O(NM)$|`a.applyFunctionSelf(std::sqrt)`|