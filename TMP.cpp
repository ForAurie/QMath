        inline void DFTRecursion(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r) {
            auto &unitRoots = get_unit_roots();
            std::vector<size_t> stk(1, 0b10);
            const size_t n = r - l; stk.reserve(log2Floor(n));
            while (stk.size()) { // 这个玩意通过阴间的进出栈条件保证了长度为 1 的永远不会入栈。。。（<=64 会跑暴力），size_t 最低为作标记，剩余作下标
                auto t = stk.back(); stk.pop_back();
                if (t & 1) stk.push_back(t << 1);
                else {
                    int d = log2Floor(t); size_t id = (t ^ (1ull << d)) >> 1; --d;
                    const size_t len = n >> d;
                    const auto &w = unitRoots[id];
                    const auto bk = l + len * (++id);
                    for (auto i = bk - len, j = bk - (len >> 1); j != bk; ++i, ++j) {
                        auto &u = *i, &v = *j; auto temp = v * w;
                        std::tie(u, v) = std::pair{u + temp, u - temp};
                    }
                    if (len > 2) stk.push_back(t | 1), stk.push_back(t << 1);
                }
            }
        }
        inline void IDFTRecursion(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r) {
            auto &unitRoots = get_unit_roots();
            std::vector<std::pair<size_t, unsigned char>> stk(1, std::make_pair(1, 0));
            const size_t n = r - l; stk.reserve(log2Floor(n));
            while (stk.size()) {
                auto [t, flag] = stk.back(); stk.pop_back();
                if (flag == 0) {
                    const int d = log2Floor(t); const size_t len = n >> d;
                    if (len == 2) {
                        size_t id = (t ^ (1ull << d));
                        const auto& w = unitRoots[id];
                        const auto bk = l + len * (++id);
                        for (auto i = bk - len, j = bk - (len >> 1); j != bk; ++i, ++j) {
                            auto &u = *i, &v = *j;
                            std::tie(u, v) = std::pair{u + v, (u - v) * w};
                        }
                    } else stk.push_back(std::make_pair(t, 1)), stk.push_back(std::make_pair(t << 1, 0));
                } else if (flag == 1) stk.push_back(std::make_pair(t, 2)), stk.push_back(std::make_pair(t << 1 | 1, 0));
                else {
                    const int d = log2Floor(t); const size_t len = n >> d; size_t id = (t ^ (1ull << d));
                    const auto& w = get_unit_roots()[id];
                    const auto bk = l + len * (++id);
                    for (auto i = bk - len, j = bk - (len >> 1); j != bk; ++i, ++j) {
                        auto &u = *i, &v = *j;
                        std::tie(u, v) = std::pair{u + v, (u - v) * w};
                    }
                }
            }
            std::reverse(l + 1, r);
        }
































































        template<size_t n>
        void __DFT(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r, size_t id = 0) {
            if constexpr (n == 1) return;
            else {
                auto mid = l + ((r - l) >> 1);
                for (auto i = l, j = mid, w = get_unit_roots()[id]; j != r; ++i, ++j) {
                    auto &u = *i, &v = *j; auto temp = v * w;
                    std::tie(u, v) = std::pair{u + temp, u - temp};
                }
                id <<= 1;
                constexpr size_t m = n >> 1;
                __DFT<m>(l, mid, id);
                __DFT<m>(mid, r, id | 1);
            }
        }
        template<size_t n>
        void __IDFT(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r, size_t id = 0) {
            if constexpr (n == 1) return;
            else {
                constexpr size_t m = n >> 1;
                auto mid = l + ((r - l) >> 1);
                __IDFT<m>(l, mid, id << 1);
                __IDFT<m>(mid, r, id << 1 | 1);
                for (auto w = get_unit_roots()[id]; mid != r; ++l, ++mid) {
                    auto &u = *l, &v = *mid;
                    std::tie(u, v) = std::pair{u + v, (u - v) * w};
                }
            }
        }
        inline void DFTinline(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r, size_t id = 0) {
            size_t n = r - l;
            switch (n) {case 128:{__DFT<128>(l, r, id);break;}case 256:{__DFT<256>(l, r, id);break;}case 512:{__DFT<512>(l, r, id);break;}case 1024:{__DFT<1024>(l, r, id);break;}case 2048:{__DFT<2048>(l, r, id);break;}case 4096:{__DFT<4096>(l, r, id);break;}case 8192:{__DFT<8192>(l, r, id);break;}case 16384:{__DFT<16384>(l, r, id);break;}case 32768:{__DFT<32768>(l, r, id);break;}case 65536:{__DFT<65536>(l, r, id);break;}case 131072:{__DFT<131072>(l, r, id);break;}case 262144:{__DFT<262144>(l, r, id);break;}case 524288:{__DFT<524288>(l, r, id);break;}case 1048576:{__DFT<1048576>(l, r, id);break;}case 2097152:{__DFT<2097152>(l, r, id);break;}case 4194304:{__DFT<4194304>(l, r, id);break;}case 8388608:{__DFT<8388608>(l, r, id);break;}case 16777216:{__DFT<16777216>(l, r, id);break;}}
        }
        inline void IDFTinline(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r, size_t id = 0) {
            size_t n = r - l;
            switch (n) {case 128:{__IDFT<128>(l, r, id);break;}case 256:{__IDFT<256>(l, r, id);break;}case 512:{__IDFT<512>(l, r, id);break;}case 1024:{__IDFT<1024>(l, r, id);break;}case 2048:{__IDFT<2048>(l, r, id);break;}case 4096:{__IDFT<4096>(l, r, id);break;}case 8192:{__IDFT<8192>(l, r, id);break;}case 16384:{__IDFT<16384>(l, r, id);break;}case 32768:{__IDFT<32768>(l, r, id);break;}case 65536:{__IDFT<65536>(l, r, id);break;}case 131072:{__IDFT<131072>(l, r, id);break;}case 262144:{__IDFT<262144>(l, r, id);break;}case 524288:{__IDFT<524288>(l, r, id);break;}case 1048576:{__IDFT<1048576>(l, r, id);break;}case 2097152:{__IDFT<2097152>(l, r, id);break;}case 4194304:{__IDFT<4194304>(l, r, id);break;}case 8388608:{__IDFT<8388608>(l, r, id);break;}case 16777216:{__IDFT<16777216>(l, r, id);break;}}
            // std::reverse(l + 1, r);
        }
        // Iteration
        inline void DFT(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r) {
            auto& unitRoots = get_unit_roots();
            for(size_t l2 = (r - l) >> 1; l2 > 1024; l2 >>= 1)
                for(auto i = l, ww = unitRoots.begin(); i != r; i += l2 << 1, ++ww) {
                    auto &w = *ww;
                    for(auto uu = i, vv = i + l2; uu != i + l2; ++uu, ++vv) {
                        auto &u = *uu, &v = *vv; auto temp = v * w;
                        std::tie(u, v) = std::pair{u + temp, u - temp};
                    }
                }
            size_t id = 0, L = std::min((size_t) 2048, (size_t) (r - l));
            for (auto i = l; i != r; i += L, ++id) DFTinline(i, i + L, id);
        }
        inline void IDFT(std::vector<TFFT>::iterator l, std::vector<TFFT>::iterator r) {
            auto& unitRoots = get_unit_roots(); const size_t n = r - l;
            size_t id = 0, L = std::min((size_t) 2048, n);
            for (auto i = l; i != r; i += L, ++id) IDFTinline(i, i + L, id);
            for(size_t l2 = 2048; l2 < n; l2 <<= 1)
                for(auto i = l, ww = unitRoots.begin(); i != r; i += l2 << 1, ++ww) {
                    auto &w = *ww;
                    for(auto uu = i, vv = i + l2; uu != i + l2; ++uu, ++vv) {
                        auto &u = *uu, &v = *vv;
                        std::tie(u, v) = std::pair{u + v, (u - v) * w};
                    }
                }
            std::reverse(l + 1, r);
        }


        inline void DFTIteration(typename std::vector<TDFT>::iterator l, typename std::vector<TDFT>::iterator r, size_t id = 0, size_t low = 2) {
            auto& unitRoots = get_unit_roots();
            low = std::max((size_t) 1, low >> 1);
            for(size_t l2 = (r - l) >> 1; l2 >= low; l2 >>= 1, id <<= 1)
                for(auto i = l, ww = unitRoots.begin() + id; i != r; i += l2 << 1, ++ww) {
                    auto &w = *ww;
                    for(auto uu = i, vv = i + l2; uu != i + l2; ++uu, ++vv) {
                        auto &u = *uu, &v = *vv; auto temp = v * w;
                        std::tie(u, v) = std::pair{u + temp, u - temp};
                    }
                }
        }
        inline void IDFTIteration(typename std::vector<TDFT>::iterator l, typename std::vector<TDFT>::iterator r, size_t id = 0, size_t low = 2) {
            auto& unitRoots = get_unit_roots(); const size_t n = r - l;
            id <<= log2Floor(n) - 1;
            for(size_t l2 = std::max((size_t) 1, low >> 1); l2 < n; l2 <<= 1, id >>= 1)
                for(auto i = l, ww = unitRoots.begin() + id; i != r; i += l2 << 1, ++ww) {
                    auto &w = *ww;
                    for(auto uu = i, vv = i + l2; uu != i + l2; ++uu, ++vv) {
                        auto &u = *uu, &v = *vv;
                        std::tie(u, v) = std::pair{u + v, (u - v) * w};
                    }
                }
        }
        #define LEN 2048
        inline void DFT(typename std::vector<TDFT>::iterator l, typename std::vector<TDFT>::iterator r) {
            DFTIteration(l, r, 0, LEN << 1);
            size_t idx = 0;
            for (auto i = l; i != r; i += LEN, ++idx) DFTIteration(i, i + LEN, idx, 2);
        }
        inline void IDFT(typename std::vector<TDFT>::iterator l, typename std::vector<TDFT>::iterator r) {
            size_t idx = 0;
            for (auto i = l; i != r; i += LEN, ++idx) IDFTIteration(i, i + LEN, idx, 2);
            IDFTIteration(l, r, 0, LEN << 1);
            std::reverse(l + 1, r);
        }
        #undef LEN





        static void ensure_precomputed(size_t n) {
            auto rt = get_unit_roots();
            get_root_size() = n;
            rt.resize(n);
            const auto wn = calcUR(n);
            for (size_t i = 1; i < n; i++) rt[i] = rt[i - 1] * wn;
        }
        inline 
        // IterationRadix4
        inline void DFT(typename std::vector<TDFT>::iterator l, typename std::vector<TDFT>::iterator r) {
            size_t len = r - l;
            if (log2Floor(len) & 1) {
                len >>= 1;
                for (auto i = l, j = l + len; j != r; ++i, ++j) {
                    auto &u = *i, &v = *j;
                    std::tie(u, v) = std::pair{u + v, u - v};
                }
            }
            const auto& unitRoots = get_unit_roots();
            const auto I = calcUR(4);
            for (size_t l4 = len >> 2; l4; len >>= 2, l4 >>= 2) {
                for (auto i = l, ww = unitRoots.begin(); i != r; i += len, ++ww) {
                    const auto &w1 = unitRoots[bit_rev];
                    const auto w2 = w1 * w1, w3 = w2 * w1;
                    for (
                        auto p0 = i, p1 = p0 + l4, p2 = p1 + l4, p3 = p2 + l4;
                        p0 != i + l4;
                        ++p0, ++p1, ++p2, ++p3
                    ) {
                        auto &x0 = *p0, &x1 = *p1, &x2 = *p2, &x3 = *p3;
                        x1 *= w1, x2 *= w2, x3 *= w3;
                        std::tie(x0, x2) = std::pair{x0 + x2, x0 - x2};
                        std::tie(x1, x3) = std::pair{x1 + x3, (x1 - x3) * I};
                        std::tie(x0, x1) = std::pair{x0 + x1, x0 - x1};
                        std::tie(x2, x3) = std::pair{x2 + x3, x2 - x3};
                    }
                }
            }
        }
        // IterationRadix4
        inline void IDFT(typename std::vector<TDFT>::iterator l, typename std::vector<TDFT>::iterator r) {
            auto& unitRoots = get_unit_roots(); const size_t n = r - l;
            size_t len = 4;
            const auto I_INV = TDFT(1) / calcUR(4);
            for(size_t l4 = 1; len <= n; l4 <<= 2, len <<= 2)
                for(auto i = l, ww = unitRoots.begin(); i != r; i += len, ++ww) {
                    const auto &w1 = TDFT(1) / *ww;
                    const auto w2 = w1 * w1, w3 = w2 * w1;
                    for (
                        auto p0 = i, p1 = p0 + l4, p2 = p1 + l4, p3 = p2 + l4;
                        p0 != i + l4;
                        ++p0, ++p1, ++p2, ++p3
                    ) { 
                        auto &x0 = *p0, &x1 = *p1, &x2 = *p2, &x3 = *p3;
                        std::tie(x0, x1) = std::pair{x0 + x1, x0 - x1};
                        std::tie(x2, x3) = std::pair{x2 + x3, (x2 - x3) * I_INV};
                        std::tie(x0, x2) = std::pair{x0 + x2, (x0 - x2) * w2};
                        std::tie(x1, x3) = std::pair{(x1 + x3) * w1, (x1 - x3) * w3};
                    }
                }
            if (len != n) {
                len >>= 2;
                for (auto i = l, j = l + len; j != r; ++i, ++j) {
                    auto &u = *i, v = *j;
                    std::tie(u, v) = std::pair{u + v, u - v};
                }
            }
            // reverse(l + 1, r);
        }