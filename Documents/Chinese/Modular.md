# Modular

## 简介

集成自动取模类型，目前有：自动取模类模板。

## 依赖

请确保 C++ 版本大于等于 C++14。

## Modualr 类模板

### 一、位置

该类型位于 `Modular.hpp` 文件的 `QMath` 命名空间中，你可以通过 `QMath:Modular<...>` 调用。

### 二、模板参数

这是它的模板参数定义：

```cpp
template <
    typename T = int, // 存储数据使用的类型，务必使用有符号类型且保证其值域不小于 (-模数 * 2 + 1, 模数 * 2 - 1)
    const T MOD = 998244353, // 默认模数
    typename MCT = unsigned long long // 在计算乘法时会临时强转到这个类型，请确保 T 和 MCT 之间存在强转定义，且其值域不小于 [0, (模数 - 1) ^ 2]，建议使用无符号类型。
> class Modular；
```

### 三、构造函数

主要的构造函数形如：
```cpp
constexpr Modular(const T& v = T());
```

可以传入一个参数，表示该类型的初始值，默认为 `T()`

除此之外还有一些没必要展示的常用构造函数。

### 四、成员函数和重载运算符

|成员函数 / 运算|作用|时间复杂度（$M$ 表示模数）|示例|
|:-:|:-:|:-:|:-:|
|`Modular& operator=(const Modular& other)`|拷贝并返回左值引用|$O(1)$|`a = b`|
|`Modular& operator=(const T& v)`|构造、赋值然后返回左值引用|$O(1)$|`a = 8`|
|`const T& getVal()`|返回对应的 `T` 类型值|$O(1)$|`a.getVal()`|
|`Modular& setVal(const T& v)`|构造、赋值然后返回左值引用|$O(1)$|`a.setVal(8)`|
|`Modular operator-()`|返回自己的相反数|$O(1)$|`-a`|
|`friend Modular operator+(const Modular& a, const Modular& b)`|计算 $a+b$、自动取模然后返回结果|$O(1)$|`a + b`|
|`friend Modular operator-(const Modular& a, const Modular& b)`|计算 $a-b$、自动取模然后返回结果|$O(1)$|`a - b`|
|`friend Modular operator*(const Modular& a, const Modular& b)`|计算 $a\times b$、自动取模然后返回结果|$O(1)$|`a * b`|
|`friend Modular operator/(const Modular& a, const Modular& b)`|计算 $a \div b$、自动取模然后返回结果（使用此功能时请确保模数是质数，且除数不为 0）|$O(\log{M})$|`a / b`|
|`Modular& operator+=(const Modular& other)`|$a\gets a+b$、自动取模然后返回左值引用|$O(1)$|`a += b`|
|`Modular& operator-=(const Modular& other)`|$a\gets a-b$、自动取模然后返回左值引用|$O(1)$|`a -= b`|
|`Modular& operator*=(const Modular& other)`|$a\gets a\times b$、自动取模然后返回左值引用|$O(1)$|`a *= b`|
|`Modular& operator/=(const Modular& other)`|$a\gets a\div b$、自动取模然后返回左值引用（使用此功能时请确保模数是质数，且除数不为 0）|$O(\log{M})$|`a /= b`|
|`template <typename U> Modular fPow(U exp)`|使用快速幂算法计算 $a^U$、自动取模然后返回结果|$O(\log{U})$|`a.fPow(8)`|
|`template <typename U> Modular& fPowSelf(U exp)`|使用快速幂算法计算 $a\gets a^U$、自动取模然后返回自己的引用|$O(\log{U})$|`a.fPowSelf(8)`|
|`Modular inv()`|计算自己的乘法逆元并返回结果（使用此功能时请确保模数是质数，且计算对象不为 0）|$O(\log{M})$|`a.inv()`|
|`Modular& invSelf()`|令自己为自己的乘法逆元并返回自己的左值引用（使用此功能时请确保模数是质数，且计算对象不为 0）|$O(\log{M})$|`a.invSelf()`|
|`friend bool operator==(const Modular& a, const Modular& b)`|两者数值相等则返回真否则返回假|$O(1)$|`a==b`|
|`friend bool operator!=(const Modular& a, const Modular& b)`|两者数值不相等则返回真否则返回假|$O(1)$|`a!=b`|
|`friend std::ostream& operator<<(std::ostream& os, const Modular& m)`|将其视作 `T` 类型输出|$O(1)$|`std::cout << a`|
|`std::istream& operator>>(std::istream& is, Modular& m)`|将其视作 `T` 类型输入|$O(1)$|`std::cin >> a`|