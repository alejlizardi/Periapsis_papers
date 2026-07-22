// blockfree_odd_sweep.cpp — exhaustive gamma-BF check at odd n via the odd
// dichotomy (Lemma W-odd): recursively enumerate ALL block-free sigma (pruned
// generation) and test PC-Hamiltonicity by the (k,j) near-perfect-matching
// acyclicity dichotomy (validated against the DP oracle, 600/600).
// Prints every failing sigma (candidate counterexample — DP-verify those).
// usage: blockfree_odd_sweep <n>
#include <bits/stdc++.h>
using namespace std;

static int n, m;
static int sigma[32], invp[32];
static bool used[32];
static long long total = 0, fails = 0;

static int par[32];
static int find(int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; }

static bool pc_ham_odd() {
    for (int i = 0; i < n; i++) invp[sigma[i]] = i;
    for (int k = 0; k <= m; k++) for (int j = 0; j <= m; j++) {
        for (int i = 0; i < n; i++) par[i] = i;
        bool acyc = true;
        for (int t = 0; t < m && acyc; t++) {
            int a = 2 * t + (t >= k ? 1 : 0);
            int ra = find(a), rb = find(a + 1);
            if (ra == rb) { acyc = false; break; }
            par[ra] = rb;
        }
        for (int t = 0; t < m && acyc; t++) {
            int v = 2 * t + (t >= j ? 1 : 0);
            int ra = find(invp[v]), rb = find(invp[v + 1]);
            if (ra == rb) { acyc = false; break; }
            par[ra] = rb;
        }
        if (acyc) return true;
    }
    return false;
}

static void rec(int pos) {
    if (pos == n) {
        total++;
        if (!pc_ham_odd()) {
            fails++;
            for (int i = 0; i < n; i++) printf("%d%c", sigma[i], i + 1 < n ? ',' : '\n');
            fflush(stdout);
        }
        return;
    }
    for (int v = 0; v < n; v++) {
        if (used[v]) continue;
        if (pos && abs(v - sigma[pos - 1]) == 1) continue;   // block-free pruning
        used[v] = true; sigma[pos] = v;
        rec(pos + 1);
        used[v] = false;
    }
}

int main(int argc, char** argv) {
    n = atoi(argv[1]); m = (n - 1) / 2;
    if (n % 2 == 0) { fprintf(stderr, "n must be odd\n"); return 1; }
    rec(0);
    fprintf(stderr, "n=%d blockfree_total=%lld ham_fails=%lld\n", n, total, fails);
    return 0;
}
