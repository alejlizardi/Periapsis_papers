// rho_tool.cpp — fast exact rho(sigma) + experiment drivers.
// Conventions match core.py / the paper: rho counts VERTICES on a properly
// colored (alternating) simple path in H_sigma. Colors: 0 = P (positions),
// 1 = V (values).
//
// Modes:
//   single <comma-perm>                 : print n, b, irr, rho, rho(quotient)
//   exh <n> [minonly]                   : exhaustive S_n joint (b,rho) stats
//   stack <comma-module> <k>            : stacked family value
//   modsearch <maxlen> <maxsum> <budget>: module slope search
//   liftstress <n> <samples> <seed>     : test rho(sigma) >= rho(tau_sigma)
//   advsearch <n> <samples> <seed>      : adversarial random block families
//
// Build: g++ -O2 -fopenmp -march=native -o rho_tool rho_tool.cpp

#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef uint64_t u64;

struct H {
    int n;
    u64 adj[2][64]; // adj[c][v] = bitmask of c-colored neighbors
};

static H build(const vector<int>& sigma) {
    H h;
    int n = (int)sigma.size();
    h.n = n;
    for (int c = 0; c < 2; c++)
        for (int v = 0; v < n; v++) h.adj[c][v] = 0;
    vector<int> pos(n);
    for (int i = 0; i < n; i++) pos[sigma[i]] = i;
    for (int i = 0; i + 1 < n; i++) {
        h.adj[0][i] |= 1ULL << (i + 1);
        h.adj[0][i + 1] |= 1ULL << i;
    }
    for (int v = 0; v + 1 < n; v++) {
        int a = pos[v], b = pos[v + 1];
        h.adj[1][a] |= 1ULL << b;
        h.adj[1][b] |= 1ULL << a;
    }
    return h;
}

// Alternation-aware reachability upper bound: from state (v, need color c),
// how many unvisited vertices could still be appended (ignoring revisit
// constraints beyond 'vis'). Returns popcount of reachable vertex set.
static int reach_bound(const H& h, int v, int c, u64 vis) {
    u64 seen[2] = {0, 0}; // seen[col] = vertices reached needing color col
    u64 frontier[2] = {0, 0};
    u64 acc = 0;
    u64 first = h.adj[c][v] & ~vis;
    frontier[1 - c] = first;
    seen[1 - c] = first;
    acc |= first;
    bool change = true;
    while (change) {
        change = false;
        for (int col = 0; col < 2; col++) {
            u64 f = frontier[col];
            frontier[col] = 0;
            while (f) {
                int w = __builtin_ctzll(f);
                f &= f - 1;
                u64 nx = h.adj[col][w] & ~vis & ~seen[1 - col];
                if (nx) {
                    seen[1 - col] |= nx;
                    frontier[1 - col] |= nx;
                    acc |= nx;
                    change = true;
                }
            }
        }
    }
    return __builtin_popcountll(acc);
}

struct DFS {
    const H* h;
    int best;
    long long nodes, budget; // budget<0 = unlimited
    bool aborted;
    void run(int v, int c, u64 vis, int len, bool prune) {
        if (aborted) return;
        if (++nodes > budget && budget > 0) { aborted = true; return; }
        if (len > best) best = len;
        u64 cand = h->adj[c][v] & ~vis;
        if (!cand) return;
        if (prune && len + reach_bound(*h, v, c, vis) <= best) return;
        while (cand) {
            int w = __builtin_ctzll(cand);
            cand &= cand - 1;
            run(w, 1 - c, vis | (1ULL << w), len + 1, prune);
        }
    }
};

// exact rho; prune=true uses reachability bound (worth it for n >~ 14)
static int rho(const H& h, bool prune = false, long long budget = -1,
               bool* ok = nullptr) {
    DFS d{&h, 1, 0, budget, false};
    for (int v = 0; v < h.n; v++)
        for (int c = 0; c < 2; c++)
            d.run(v, c, 1ULL << v, 1, prune);
    if (ok) *ok = !d.aborted;
    return d.best;
}

static int rho_perm(const vector<int>& s, bool prune = false,
                    long long budget = -1, bool* ok = nullptr) {
    return rho(build(s), prune, budget, ok);
}

// ---- block machinery ----
static vector<pair<int,int>> blocks_of(const vector<int>& s) {
    int n = (int)s.size();
    vector<pair<int,int>> out;
    int i = 0;
    while (i < n) {
        int j = i;
        if (j + 1 < n && abs(s[j + 1] - s[j]) == 1) {
            int eps = s[j + 1] - s[j];
            while (j + 1 < n && s[j + 1] - s[j] == eps) j++;
        }
        out.push_back({i, j});
        i = j + 1;
    }
    return out;
}

static vector<int> quotient_of(const vector<int>& s) {
    auto bl = blocks_of(s);
    int b = (int)bl.size();
    vector<int> lows(b), idx(b), rank_(b);
    for (int k = 0; k < b; k++) {
        int lo = INT_MAX;
        for (int i = bl[k].first; i <= bl[k].second; i++) lo = min(lo, s[i]);
        lows[k] = lo;
    }
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int c) { return lows[a] < lows[c]; });
    for (int r = 0; r < b; r++) rank_[idx[r]] = r;
    return rank_;
}

static vector<int> stack_perm(const vector<int>& module, int k) {
    vector<int> sizes;
    sizes.push_back(2);
    for (int r = 0; r < k; r++)
        for (int x : module) sizes.push_back(x);
    sizes.push_back(3);
    sizes.push_back(2);
    vector<int> s;
    int base = 0;
    for (int sz : sizes) {
        for (int v = base + sz - 1; v >= base; v--) s.push_back(v);
        base += sz;
    }
    return s;
}

static vector<int> parse_perm(const string& str) {
    vector<int> out;
    string cur;
    for (char ch : str) {
        if (ch == ',') { out.push_back(stoi(cur)); cur.clear(); }
        else cur += ch;
    }
    if (!cur.empty()) out.push_back(stoi(cur));
    return out;
}

static string show(const vector<int>& s) {
    string o = "[";
    for (size_t i = 0; i < s.size(); i++) {
        if (i) o += ",";
        o += to_string(s[i]);
    }
    return o + "]";
}

// ---- exhaustive mode ----
static void exhaustive(int n) {
    // joint histogram cnt[b][rho], plus exemplars of per-b minima
    vector<vector<u64>> cnt(n + 1, vector<u64>(n + 1, 0));
    vector<int> minrho_b(n + 1, INT_MAX);
    vector<vector<int>> exemplar(n + 1);
    long long chunks = (long long)n * (n - 1);
#pragma omp parallel
    {
        vector<vector<u64>> lc(n + 1, vector<u64>(n + 1, 0));
#pragma omp for schedule(dynamic)
        for (long long ch = 0; ch < chunks; ch++) {
            int a = (int)(ch / (n - 1)), bvel = (int)(ch % (n - 1));
            vector<int> rest;
            for (int x = 0; x < n; x++) if (x != a) rest.push_back(x);
            int bval = rest[bvel];
            rest.erase(rest.begin() + bvel);
            vector<int> s(n);
            s[0] = a; s[1] = bval;
            do {
                for (int i = 2; i < n; i++) s[i] = rest[i - 2];
                H h = build(s);
                int r = rho(h);
                int bb = (int)blocks_of(s).size();
                lc[bb][r]++;
                if (r < minrho_b[bb]) {
#pragma omp critical
                    if (r < minrho_b[bb]) { minrho_b[bb] = r; exemplar[bb] = s; }
                }
            } while (next_permutation(rest.begin(), rest.end()));
        }
#pragma omp critical
        for (int i = 0; i <= n; i++)
            for (int j = 0; j <= n; j++) cnt[i][j] += lc[i][j];
    }
    printf("n=%d joint (b,rho) table; rows=b, min rho per b, count at min, exemplar\n", n);
    int globalmin = INT_MAX;
    for (int bb = 1; bb <= n; bb++) {
        u64 tot = 0; int mn = -1; u64 atmin = 0;
        for (int r = 1; r <= n; r++) {
            tot += cnt[bb][r];
            if (cnt[bb][r] && mn < 0) { mn = r; atmin = cnt[bb][r]; }
        }
        if (!tot) continue;
        globalmin = min(globalmin, mn);
        printf("b=%2d total=%12llu minrho=%2d count@min=%10llu ex=%s\n",
               bb, (unsigned long long)tot, mn, (unsigned long long)atmin,
               show(exemplar[bb]).c_str());
    }
    printf("rho_min(%d) = %d\n", n, globalmin);
    // full histogram dump for downstream analysis
    printf("HIST b rho count\n");
    for (int bb = 1; bb <= n; bb++)
        for (int r = 1; r <= n; r++)
            if (cnt[bb][r])
                printf("H %d %d %llu\n", bb, r, (unsigned long long)cnt[bb][r]);
}

// ---- module search ----
static void modsearch(int maxlen, int maxsum, long long budget) {
    // enumerate modules over sizes {1,2,3,4} (canonical: skip rotations? no —
    // stacking is order-sensitive, keep all), compute slope estimate
    vector<vector<int>> mods;
    vector<int> cur;
    function<void(int)> gen = [&](int sum) {
        if (!cur.empty()) mods.push_back(cur);
        if ((int)cur.size() >= maxlen) return;
        for (int s = 1; s <= 4; s++) {
            if (sum + s > maxsum) continue;
            cur.push_back(s);
            gen(sum + s);
            cur.pop_back();
        }
    };
    gen(0);
    printf("modules to test: %zu\n", mods.size());
    struct Res { double slope; vector<int> mod; int r1, r2, n1, n2; };
    vector<Res> results;
#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < (long long)mods.size(); i++) {
        const auto& M = mods[i];
        vector<int> s1 = stack_perm(M, 1), s2 = stack_perm(M, 2);
        if ((int)s2.size() > 62) continue;
        bool ok1, ok2;
        int r1 = rho_perm(s1, true, budget, &ok1);
        if (!ok1) continue;
        int r2 = rho_perm(s2, true, budget, &ok2);
        if (!ok2) continue;
        int msize = 0;
        for (int x : M) msize += x;
        double slope = double(r2 - r1) / msize;
#pragma omp critical
        results.push_back({slope, M, r1, r2, (int)s1.size(), (int)s2.size()});
    }
    sort(results.begin(), results.end(),
         [](const Res& a, const Res& b) { return a.slope < b.slope; });
    printf("top 40 lowest slope estimates (rho(k=2)-rho(k=1))/|M|:\n");
    for (int i = 0; i < (int)results.size() && i < 40; i++) {
        auto& R = results[i];
        printf("slope=%.4f M=%s rho@k1=%d(n=%d) rho@k2=%d(n=%d)\n", R.slope,
               show(R.mod).c_str(), R.r1, R.n1, R.r2, R.n2);
    }
}

// ---- random structured sigma generator ----
static vector<int> random_block_perm(int n, mt19937_64& rng) {
    // random composition of n into blocks of size 1..4, random value
    // arrangement of blocks, each block asc/desc at random
    vector<int> sizes;
    int rem = n;
    while (rem > 0) {
        int s = 1 + (int)(rng() % 4);
        if (s > rem) s = rem;
        sizes.push_back(s);
        rem -= s;
    }
    int b = (int)sizes.size();
    vector<int> order(b);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    // value base offset for each block: blocks get contiguous value ranges
    // in the shuffled order
    vector<int> base(b);
    int acc = 0;
    vector<int> rank_(b);
    for (int r = 0; r < b; r++) rank_[order[r]] = r;
    // careful: assign value ranges by rank in shuffled order
    vector<int> base_by_rank(b);
    {
        int a2 = 0;
        for (int r = 0; r < b; r++) {
            base_by_rank[r] = a2;
            a2 += sizes[order[r]];
        }
    }
    vector<int> s;
    for (int k = 0; k < b; k++) {
        int sz = sizes[k];
        int lo = base_by_rank[rank_[k]];
        bool desc = rng() & 1;
        if (sz == 1) s.push_back(lo);
        else if (desc)
            for (int v = lo + sz - 1; v >= lo; v--) s.push_back(v);
        else
            for (int v = lo; v < lo + sz; v++) s.push_back(v);
    }
    return s;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "mode required\n"); return 1; }
    string mode = argv[1];
    if (mode == "single") {
        vector<int> s = parse_perm(argv[2]);
        auto bl = blocks_of(s);
        vector<int> tau = quotient_of(s);
        bool irr = (int)bl.size() == (int)s.size();
        printf("n=%zu b=%zu irr=%d rho=%d rho_quotient=%d tau=%s\n", s.size(),
               bl.size(), irr ? 1 : 0, rho_perm(s, s.size() > 13),
               tau.size() > 1 ? rho_perm(tau, tau.size() > 13) : 1,
               show(tau).c_str());
    } else if (mode == "exh") {
        exhaustive(atoi(argv[2]));
    } else if (mode == "stack") {
        vector<int> M = parse_perm(argv[2]);
        int k = atoi(argv[3]);
        vector<int> s = stack_perm(M, k);
        printf("n=%zu sigma=%s\n", s.size(), show(s).c_str());
        printf("b=%zu rho=%d\n", blocks_of(s).size(), rho_perm(s, true));
    } else if (mode == "modsearch") {
        modsearch(atoi(argv[2]), atoi(argv[3]), argc > 4 ? atoll(argv[4]) : 200000000LL);
    } else if (mode == "liftstress") {
        int n = atoi(argv[2]);
        long long samples = atoll(argv[3]);
        u64 seed = argc > 4 ? atoll(argv[4]) : 12345;
        long long viol = 0, tested = 0;
#pragma omp parallel reduction(+ : viol, tested)
        {
#ifdef _OPENMP
            mt19937_64 rng(seed + 1000 * omp_get_thread_num());
#else
            mt19937_64 rng(seed);
#endif
#pragma omp for schedule(dynamic, 64)
            for (long long i = 0; i < samples; i++) {
                vector<int> s = random_block_perm(n, rng);
                vector<int> tau = quotient_of(s);
                if ((int)tau.size() == (int)s.size() || tau.size() < 2) continue;
                int rs = rho_perm(s, n > 13);
                int rt = rho_perm(tau, tau.size() > 13);
                tested++;
                if (rs < rt) {
                    viol++;
#pragma omp critical
                    printf("VIOLATION rho=%d rho_tau=%d sigma=%s\n", rs, rt,
                           show(s).c_str());
                }
            }
        }
        printf("liftstress n=%d tested=%lld violations=%lld\n", n, tested, viol);
    } else if (mode == "advsearch") {
        int n = atoi(argv[2]);
        long long samples = atoll(argv[3]);
        u64 seed = argc > 4 ? atoll(argv[4]) : 999;
        int globalmin = INT_MAX;
        vector<int> bestperm;
#pragma omp parallel
        {
#ifdef _OPENMP
            mt19937_64 rng(seed + 1000 * omp_get_thread_num());
#else
            mt19937_64 rng(seed);
#endif
#pragma omp for schedule(dynamic, 16)
            for (long long i = 0; i < samples; i++) {
                vector<int> s = random_block_perm(n, rng);
                bool ok;
                int r = rho_perm(s, true, 50000000, &ok);
                if (!ok) continue;
                if (r < globalmin) {
#pragma omp critical
                    if (r < globalmin) {
                        globalmin = r;
                        bestperm = s;
                        printf("new min rho=%d (n=%d, ratio %.4f) b=%zu sigma=%s\n",
                               r, n, double(r) / n, blocks_of(s).size(),
                               show(s).c_str());
                        fflush(stdout);
                    }
                }
            }
        }
        printf("advsearch n=%d done: min rho=%d ratio=%.4f\n", n, globalmin,
               double(globalmin) / n);
    } else {
        fprintf(stderr, "unknown mode\n");
        return 1;
    }
    return 0;
}
