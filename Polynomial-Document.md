## 简介

封装了一个多项式模板类，已上传至 [GitHub 仓库的 Polynomial 文件](https://github.com/ForAurie/QMath/blob/main/Polynomial)，`include` 该文件或粘贴至代码开头即可使用，最低支持 C++17，如果你发现了 BUG 或你有一些建议，欢迎在 GitHub 提出，我看到后会及时修改。

本模板类具有以下几大优点：

1. 优良的通用性

    可以通过简单的模板参数修改，在 FFT 和 NTT 模式之间无缝切换，你甚至可以把定义完整的高精度类和这个模板搭配使用。

2. 优良的速度

    再保证支持泛型的情况下，它可以进入[模板题最优解第一页](https://www.luogu.com.cn/record/list?pid=P3803&orderBy=1&status=&page=1)，无论是 FFT 模式，还是 NTT 模式，都可以稳定在 490ms 以内，比它速度更优的提交都针对模数、`int` 数据类型做了特殊优化，不完全支持泛型。

3. 精致的封装

    充分发挥面向对象编程思想，该类继承自 `std::vector`，仿照 STL 进行封装，可以和 STL 无缝衔接，~~你甚至可以把它当 `std::vector` 用。。。用法和 `std::vector` 一模一样。~~ 本模板不会暴露出任何对用户没有用的成员函数。

4. 强大的功能

    它顺便实现了多项式 `exp`、`ln`、`sqrt` 等功能。还使用 FWT、FMT 实现了与、或、异或卷积，不久的将来会添加集合幂级数相关功能。

6. 递归实现与迭代实现的混合运用

    技术实现相关，介绍完用法后会详细讨论。

## 用法

以下内容涉及的代码过于细碎，难以测试，因此并没有经过充分的测试，只是为了展示基本的使用方法，你可以当作伪代码，仓库里的代码是经过了测试的，至少通过了每个模板题。

它位于 `Polynomial` 文件，`QMath` 命名空间。

### 一、模板参数和构造函数介绍

这是它的模板参数定义：

```cpp
template <
    typename T = double, // 默认使用 double 类型存储数据
    typename TDFT = std::complex<double>, // 默认使用 std::complex<double> 做变换
    auto UR = expn // 单位根计算函数，接受一个 size_t 输入值类型 n，返回 TDFT 类型的 n 次单位根，由于 Polynomial 类的特殊设计，所以无需实现求单位根的倒数这一功能。
    auto T2TDFT = T2TFFT, // 需要提供 T（存储时的） 到 TDFT（变换时的类型）的转换函数，2 为 to 的谐音
    auto TDFT2T = TFFT2T, // 同上
> class Polynomial : public std::vector<T>; // 继承 std::vector
```

如果给定上述模板参数的函数，就无需确定变换是 FFT 还是 NTT，因此实现了一个模板里同时支持多种变换的方式。**当 T 类型和 TDFT 类型相同时，会自动省去转换函数的调用，此时可以随意传入或干脆不传入。**

`Polynomial` 文件的 `QMath` 命名空间提供了一组上述模板参数需要的函数，并作为默认函数，它们的分别为：`QMath::T2TFFT`、`QMath::TFFT2T`、`QMath::expn`。

它们的实现为：

```cpp
constexpr double PI2 = 6.283185307179586476925286766559005768394338798750211641949889;
std::complex<double> expn(size_t n) { return std::complex<double>(std::cos(PI2 / n), std::sin(PI2 / n)); }
std::complex<double> T2TFFT(double x) { return std::complex<double>(x, 0); }
double TFFT2T(const std::complex<double>& x) { return x.real(); }
```

非常简短，就是普通 FFT 所需的实现。

当然，这个模板也支持 NTT，不过需要准备一个自动取模类型，[这里提供了一个](https://github.com/ForAurie/QMath/blob/main/Modular)和 Polynomial 在同一个仓库的类型，**如果你在测试 Polynomial 模板类性能时不使用这个自动取模模板类，而是使用别的，可能由于自动取模模板类过慢而达不到预期效果**。

该自动取模类型位于 `Modular` 文件的 `QMath` 命名空间。模板参数定义为：

```cpp
template <
    typename T = int, // 存储数据使用的类型，务必使用有符号类型且保证其值域大于模数的一倍，不然会出错
    const T MOD = 998244353, // 默认模数
    typename MCT = unsigned long long // 计算乘法时转换的临时类型，在计算乘法时会临时强转到这个类型，请确保 T 和 MCT 之间存在强转定义且其值域大于模数的平方，建议使用无符号类型。
> class Modular；

// 如果要构造该类型可以这样
QMath::Modular<int, 998244353, unsigned long long> a(1); // 初始化为 1
// 如果要修改值建议使用显式函数来降低用错的风险
a.setVal(2);
// 访问数值请使用显式函数
a.getVal();
```

该类型使用快速幂实现单 $\log$ 除法，要想使除法功能保持正确，请确保模数为质数。该类型重载了输入输出流。不久的将来会提供编译期自动检测模数是否是质数的功能和 `exgcd` 实现求逆元，之前其实实现过 `exgcd` 求逆元，但是手贱给删了。

下面提供一组基于模数 `998244353` 的 Polynomial 类型定义的实现：

```cpp
using Mint = QMath::Modular<>; // 定义自动取模类型并使用默认模板参数
inline Mint T2TNTT(Mint x) { return x; }
inline Mint TNTT2T(Mint x) { return x; }
inline Mint NTT(size_t x) { return Mint(3).fPow(998244352 / x); }

typedef QMath::Polynomial<Mint, Mint, NTT, T2NTT, NTT2T> Poly;
```

给出模板题的代码：

```cpp
#include "Polynomial"
#include "Modular"
// 为保证此处代码的简洁，不再将模板粘过来而是使用 include 代替。
#include <bits/stdc++.h>

using namespace std;
using Mint = QMath::Modular<>;
inline Mint T2TNTT(Mint x) { return x; }
inline Mint TNTT2T(Mint x) { return x; }
inline Mint NTT(size_t x) { return Mint(3).fPow(998244352 / x); }

int main() {
    ios::sync_with_stdio(false), cin.tie(nullptr);
    int n, m;
    cin >> n >> m;
    QMath::Polynomial<Mint, Mint, NTT, T2TNTT, TNTT2T> a(n + 1), b(n + 1, Mint(0)); // 构造函数和 vector 是一样的，这两种都可以。
    for (auto& i : a) cin >> i; // 精致封装的体现，可以像 vector 一样使用，同时 Modular 重载了输入输出流，如果使用 cin、cout 可以直接输入输出。
    for (auto& i : b) cin >> i;
    // 你甚至可以对 a 离散化：
    // sort(a.begin(), a.end());
    // a.erase(unique(a.begin(), a.end()), a.end());
    a *= b;
    for (auto i : a) cout << i << ' ';
    return 0; // 返回零就通过了模板题！超短代码（不算模板） AC 多项式乘法模板并获得极快的速度！
}
```

需要注意的是，本模板会 $O(n)$ 预处理要用到的单位根（由于使用特殊优化，需要快速随机访问任意单位根，所以不能 $O(\log{n})$ 预处理，好像也能？只是有点麻烦？不过这样干精度误差可能比较大？对 FFT 不友好？改天仔细试试，如果优化明显就更新到代码仓库），并以 `thread_local static` 的形式存储下来和相同类型共享，你可以通过成员函数 `clearCache()` 清空存储下来的单位根。

### 二、其成员函数和重载的运算符介绍


|成员函数 / 运算|作用|复杂度（$n$ 表示多项式长度）|示例|
|:-:|:-:|:-:|:-:|
|所有从 `std::vector` 继承过来的东西|略|略|略|
|`void clearCache()`|清空预处理的单位根|$O(n)$|`a.clearCache()`|
|`Polynomial& operator=(const Polynomial& o)`|拷贝并返回拷贝到的对象的引用|$O(n)$|`a = b`|
|`Polynomial derivative()`|返回自己的导数多项式|$O(n)$|`a.derivative()`|
|`Polynomial& derivativeSelf()`|令自己为自己的导数多项式并返回左值引用（常数项补零）|$O(n)$|`a.derivativeSelf()`|
|`Polynomial integral()`|返回自己的积分多项式|$O(n)$|`a.integral()`|
|`Polynomial& integralSelf()`|令自己为自己的积分多项式并返回自己的引用（常数项补零）|$O(n)$|`a.integralSelf()`|
|`T calc(const T &x)`|计算 $x$ 带入自己的结果|$O(n)$|`a.calc(1)`|
|`T calcDerivative(const T &x)`|计算 $x$ 带入自己的导数多项式的结果|$O(n)$|`a.calcDerivative(2)`|
|`T calcIntegral(const T &x)`|计算 $x$ 带入自己的积分多项式的结果|$O(n)$|`a.calcIntegral(3)`|
|`Polynomial operator+(const Polynomial& o) const`|多项式加法，自动拉伸成两者较长长度，返回结果|$O(n)$|`a = a + b`|
|`Polynomial& operator+=(const Polynomial& o)`|多项式加等于，自动拉伸成两者较长长度，返回左值引用|$O(n)$|`a += b`|
|`Polynomial operator-(const Polynomial& o) const`|多项式减法，自动拉伸成两者较长长度，返回结果|$O(n)$|`a = a - b`|
|`Polynomial& operator-=(const Polynomial& o)`|多项式减等于，自动拉伸成两者较长长度，返回左值引用|$O(n)$|`a -= b`|
|`Polynomial operator*(const T& o) const`|多项式乘数字，返回结果|$O(n)$|`a = a * 4`|
|`Polynomial& operator*=(const T& o)`|多项式乘等于数字，返回左值引用|$O(n)$|`a *= 4`|
|`Polynomial operator/(const T& o) const`|多项式除数字，返回结果|$O(n)$|`a = a / 5`|
|`Polynomial& operator/=(const T& o)`|多项式除等于数字，返回左值引用|$O(n)$|`a /= 5`|
|`Polynomial operator*(const Polynomial& o) const`|多项式乘多项式，返回结果|$O(n\log{n})$|`a = a * b`|
|`Polynomial& operator*=(const Polynomial& o)`|多项式乘等于多项式，返回左值引用|$O(n\log{n})$|`a *= b`|
|`Polynomial& operator%=(size_t n)`|类似于 `.resize(n)`，借助多项式的模定义，更加直观（无论原长度比 $n$ 大还是比 $n$ 小，都统一拉伸到 $n$）|$O(n)$|`a %= 6`|
|`Polynomial operator%(size_t n) const`|不修改自己，返回一个 `.resize(n)` 的结果|$O(n)$|`a = a % 7`|
|`Polynomial inv() const`|返回自己的乘法逆|$O(n\log{n})$|`a.inv()`|
|`Polynomial& invSelf()`|令自己为自己的乘法逆，返回自己的引用|$O(n\log{n})$|`a.invSelf()`|
|`Polynomial operator/(const Polynomial& o) const`|使用多项式乘法逆实现|$O(n\log{n})$|`a = a / b`|
|`Polynomial& operator/=(const Polynomial& o)`|使用多项式乘法逆实现|$O(n\log{n})$|`a /= b`|
|`template<auto LN = __return0> Polynomial ln() const`|返回自己求 $\ln$ 的结果（若不保证零次项为 $1$ 则需传入一个用于 `T` 类型求 $\ln$ 的函数，若不传入，默认提供恒返回 `T(0)` 的函数）|$O(n\log{n})$|`a.ln()`|
|`template<auto LN = __return0> Polynomial& lnSelf()`|令自己为自己求 $\ln$ 的结果（若不保证零次项为 $1$ 则需传入一个用于 `T` 类型求 $\ln$ 的函数，若不传入，默认提供恒返回 `T(0)` 的函数），返回自己的引用|$O(n\log{n})$|`a.lnSelf()`|
|`template<auto EXP = __return1> Polynomial exp() const`|返回自己求 $\exp$ 的结果（若不保证零次项为 $0$ 则需传入一个用于 `T` 类型求 $\exp$ 的函数，若不传入，默认提供恒返回 `T(1)` 的函数）|$O(n\log{n})$|`a.exp<exp>()`|
|`template<auto EXP = __return1> Polynomial& expSelf()`|令自己为自己求 $\exp$ 的结果（若不保证零次项为 $0$ 则需传入一个用于 `T` 类型求 $\exp$ 的函数，若不传入，默认提供恒返回 `T(1)` 的函数），返回自己的引用|$O(n\log{n})$|`a.expSelf<exp>()`|
|`template<auto SQRT = __return1> Polynomial sqrt() const`|返回自己的平方根（若不保证零次项为 $1$ 则需传入一个用于 `T` 类型求平方根的函数，若不传入，默认提供恒返回 `T(1)` 的函数）|$O(n\log{n})$|`a.sqrt()`|
|`template<auto SQRT = __return1> Polynomial& sqrtSelf()`|令自己为自己的平方根（若不保证零次项为 $1$ 则需传入一个用于 `T` 类型求平方根的函数，若不传入，默认提供恒返回 `T(1)` 的函数），返回自己的引用|$O(n\log{n})$|`a.sqrtSelf()`|
|`template<typename U> Polynomial pow(U n, U m = U(-1))`|求自己的 $n$ 次幂，无需保证零次项为 $1$ 但必须保证零次项不为 $0$，若在求次幂时需要对指数取模来缩小范围，则 $n$ 应传入指数对模数 $MOD$ 取模的结果，$m$ 应传入指数对 $\varphi(MOD)$ 取模的结果，否则你可以忽略 $m$|$O(n\log{n})$|`a.pow(10)`|
|`template<typename U> Polynomial& powSelf(U n, U m = U(-1))`|令自己为自己的 $n$ 次幂，无需保证零次项为 $1$ 但必须保证零次项不为 $0$，若在求次幂时需要对指数取模来缩小范围，则 $n$ 应传入指数对模数 $MOD$ 取模的结果，$m$ 应传入指数对 $\varphi(MOD)$ 取模的结果，否则你可以忽略 $m$|$O(n\log{n})$|`a.powSelf(9982, 44353)`|
|`Polynomial operator\|(const Polynomial& o) const`|求或卷积并返回结果，自动拉伸至最小的大于等于两者的二的整次幂|$O(n\log{n})$|`a = a \| b`|
|`Polynomial& operator\|=(const Polynomial& o)`|求或卷积并返回左值引用，自动拉伸至最小的大于等于两者的二的整次幂|$O(n\log{n})$|`a \|= b`|
|`Polynomial operator&(const Polynomial& o) const`|求与卷积并返回结果，自动拉伸至最小的大于等于两者的二的整次幂|$O(n\log{n})$|`a = a & b`|
|`Polynomial& operator&=(const Polynomial& o)`|求与卷积并返回左值引用，自动拉伸至最小的大于等于两者的二的整次幂|$O(n\log{n})$|`a &= b`|
|`Polynomial operator^(const Polynomial& o) const`|求异或卷积并返回结果，自动拉伸至最小的大于等于两者的二的整次幂|$O(n\log{n})$|`a = a ^ b`|
|`Polynomial& operator^=(const Polynomial& o)`|求异或卷积并返回左值引用，自动拉伸至最小的大于等于两者的二的整次幂|$O(n\log{n})$|`a ^= b`|
|`friend std::ostream& operator<<(std::ostream& os, const Polynomial& p)`|输出流重载，输出形式：`f(x) = ax ^ 0 + bx ^ 1 + cx ^ 2......`|$O(n)$|`std::cout << a`|

## 技术实现分享

阅读了“[递归算法真的要比迭代慢吗？](https://www.luogu.com.cn/article/tfoqhji5)”，发现递归算法实际更 Cache-friendly，只是递归开销大，此文解决递归开销大的方案是使用模板递归让编译器在编译阶段展开递归函数，效果确实是明显的。但有两个问题：

* 实测加上这一优化后模板题的表现：FFT 从 580 ms 左右优化至 520 ms 左右，但 NTT 从 480 ms 左右劣化到 510 ms 左右。

    问题原因分析：迭代式 Cache-unfriendly 的原因是因为每次换层都要重新开始遍历多项式数组，需要把多项式末尾的从缓存移出，开头的移进去（能否分奇偶层使用不同方向遍历呢？改天试试）。递归式则由于优秀的递归顺序避免了这个问题，NTT 使用 `int` 存储，`int` 占空间小，在题目数据范围下几乎整个多项式都能被扔到 L2 / L1 缓存里面，因此优化 Cache 对其效果不明显，反而会，但 FFT 使用了 `std::complex<double>` 存储，空间占用是 `int` 的 $4$ 倍，缓存放不下，因此 Cache 对其优化效果明显。

    解决方案：可以短序列用迭代、长序列用递归。但只能这样吗？止步于此了吗？虽然 NTT / FFT 的变换实现分为迭代式和递归式，但它们做的事情本质相同，只是访问顺序不同，为什么不能在一个序列上同时使用迭代和递归呢？因此我的混合使用了迭代式和递归式：先使用递归实现，进入递归时判断目前递归区间的长度，若长度很大，则继续递归，如果长度小于一个阈值（本代码取 $4096$），短到即便是使用迭代式实现也能完全卡进缓存，则停止递归剩下的全部迭代实现，此时递归结构可能如下：

    ```
       ----------------
     --------    --------
    ----  ----  ----  ----
    ++++  ++++  ++++  ++++
    ++++  ++++  ++++  ++++
    ++++  ++++  ++++  ++++
    ```

    `-` 部分区间长度大，采用递归。`+` 部分区间很小，直接迭代完成。可以看到递归树顶部是递归算的，底部是迭代算的，虽然遍历方式不一样，但做的事情都一样，因此不影响正确性。[夺取多项式全家桶最优解手把手教程](https://www.luogu.com.cn/article/k9j38kqv) 貌似也提到了这个东西，但比较笼统，不知道意思是长的原序列全迭代，短的原序列全递归，还是把递归树拆开，同长度原序列递归迭代混用。

* 支持的多项式长度范围有限，因为编译期间展开长多项式的递归函数会导致编译文件体积暴增，因此不能支持特别长的多项式，虽在 $998244353$ 下可以接受，但是违背了泛型思想。

    解决方案：如果发现多项式过长，则在递归结构的最高最长的几层放弃优化，直接使用迭代计算，待问题规模缩小至递归可以接受时开始递归，本代码可递归长度的上限是 $2^{24}$，洛谷模板题其实没有超过这个上限。

**这是最终的优化效果**

|算法|NTT|FFT|
|:-:|:-:|:-:|
|纯迭代版|https://www.luogu.com.cn/record/264397472 478ms|https://www.luogu.com.cn/record/264424326 572ms|
|混用版|https://www.luogu.com.cn/record/264423088 477ms|https://www.luogu.com.cn/record/264422355 473ms|

可以看到在不影响 NTT 速度的情况下 FFT 被优化到了和 NTT 相同的速度，测试代码和最终仓库代码可以有所不同，但没有对速度影响明显的修改，大多是码风修改或格式修改。

**一些其它用到的小优化**

* 若 `TDFT` = `std::complex<T>` 则自动使用快速傅里叶变换的三变二优化，如果你不使用 `std::complex<>` 则享受不到此优化，建议使用 `std::complex<>`！！！

* 若 `T` = `TDFT`，则省略掉 `T2TDFT()` 和 `TDFT2T()` 函数的使用，你可以随意传入或不传入这两个函数，建议别传。

* 若乘法运算时形如 `a *= a` 或 `a = a * a`，则在正变换时只会做一次变换。