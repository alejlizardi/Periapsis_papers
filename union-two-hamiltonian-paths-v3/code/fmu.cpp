// fmu.cpp — WP-J T2b: exact max-mu engine for G(n,sigma) and the first f(n) data.
//
// G(n,sigma) per paper Def def:Gns: vertices v_0..v_{n-1} (ids 0..n-1, top path)
// and u_0..u_{n-1} (ids n..2n-1, bottom path); matching M = {v_i u_{sigma(i)}}.
// mu(path) = #matching edges on the simple path; maxmu(sigma) = max over simple
// paths; f(n) = min_sigma maxmu(sigma).
//
// Exact per-sigma computation: bitmask DP over (vertex-set, endpoint) states —
// dp[mask][e] = max matching edges of a simple path with vertex set mask ending
// at e; feasible for 2n <= ~20.  Validated by mode checkdp (DP vs plain DFS
// exhaustive over S_n) and by the Conversion Lemma guard maxmu >= rho (mode
// slack, vs the validated 2-color rho engine copied verbatim from census.cpp).
//
// Build: g++ -O3 -fopenmp -march=native -o fmu fmu.cpp
#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef uint64_t u64;
typedef uint32_t u32;

// ---------- 2-color rho core (VERBATIM census.cpp port) for the slack guard ----------
struct H { int n; u64 adj[2][64]; };
static H build(const vector<int>& sigma) {
    H h; int n = (int)sigma.size(); h.n = n;
    for (int c = 0; c < 2; c++) for (int v = 0; v < n; v++) h.adj[c][v] = 0;
    vector<int> pos(n);
    for (int i = 0; i < n; i++) pos[sigma[i]] = i;
    for (int i = 0; i + 1 < n; i++) { h.adj[0][i] |= 1ULL<<(i+1); h.adj[0][i+1] |= 1ULL<<i; }
    for (int v = 0; v + 1 < n; v++) { int a=pos[v],b=pos[v+1]; h.adj[1][a]|=1ULL<<b; h.adj[1][b]|=1ULL<<a; }
    return h;
}
static int reach_bound(const H& h, int v, int c, u64 vis) {
    u64 seen[2]={0,0}, frontier[2]={0,0}, acc=0;
    u64 first = h.adj[c][v] & ~vis;
    frontier[1-c]=first; seen[1-c]=first; acc|=first;
    bool change=true;
    while (change) { change=false;
        for (int col=0; col<2; col++) {
            u64 f=frontier[col]; frontier[col]=0;
            while (f) { int w=__builtin_ctzll(f); f&=f-1;
                u64 nx=h.adj[col][w]&~vis&~seen[1-col];
                if (nx){ seen[1-col]|=nx; frontier[1-col]|=nx; acc|=nx; change=true; } } } }
    return __builtin_popcountll(acc);
}
struct DFS2 {
    const H* h; int best;
    void run(int v,int c,u64 vis,int len,bool prune){
        if(len>best)best=len;
        u64 cand=h->adj[c][v]&~vis; if(!cand)return;
        if(prune&&len+reach_bound(*h,v,c,vis)<=best)return;
        while(cand){int w=__builtin_ctzll(cand);cand&=cand-1;run(w,1-c,vis|(1ULL<<w),len+1,prune);} }
};
static int rho_perm(const vector<int>& s, bool prune=true){
    H h=build(s); DFS2 d{&h,1};
    for(int v=0;v<h.n;v++)for(int c=0;c<2;c++)d.run(v,c,1ULL<<v,1,prune);
    return d.best;
}

// ---------- G(n,sigma) adjacency ----------
// vertex ids: 0..n-1 top (v_i), n..2n-1 bottom (u_{j} = n+j).
struct G {
    int n, N;                 // N = 2n
    u32 adj[40];              // adjacency bitmask over N<=32 vertices
    u32 mate[40];             // mate[v] = bit of the matched partner (single bit) or 0
};
static G buildG(const vector<int>& sigma){
    G g; int n=(int)sigma.size(); g.n=n; g.N=2*n;
    for(int v=0;v<g.N;v++){ g.adj[v]=0; g.mate[v]=0; }
    for(int i=0;i+1<n;i++){ g.adj[i]|=1u<<(i+1); g.adj[i+1]|=1u<<i; }
    for(int j=0;j+1<n;j++){ g.adj[n+j]|=1u<<(n+j+1); g.adj[n+j+1]|=1u<<(n+j); }
    for(int i=0;i<n;i++){ int u=n+sigma[i];
        g.adj[i]|=1u<<u; g.adj[u]|=1u<<i;
        g.mate[i]=1u<<u; g.mate[u]=1u<<i; }
    return g;
}

// ---------- exact maxmu: bitmask DP ----------
// dp indexed by [mask][endpoint]: max mu over simple paths with vertex set=mask
// ending at endpoint. int8 values (mu <= n <= 127).
static int maxmu_dp(const G& g, vector<int8_t>& dp){
    int N=g.N; size_t SZ=((size_t)1<<N)*N;
    if(dp.size()<SZ) dp.assign(SZ,-1);
    else fill(dp.begin(),dp.begin()+SZ,-1);
    int best=0;
    // seed single vertices
    for(int v=0;v<N;v++) dp[((size_t)(1u<<v))*N+v]=0;
    for(u32 mask=1; mask< (1u<<N); mask++){
        u32 mm=mask;
        size_t base=(size_t)mask*N;
        while(mm){
            int e=__builtin_ctz(mm); mm&=mm-1;
            int8_t cur=dp[base+e];
            if(cur<0) continue;
            if(cur>best) best=cur;
            u32 nx=g.adj[e]&~mask;
            while(nx){
                int w=__builtin_ctz(nx); nx&=nx-1;
                int8_t add=(g.mate[e]>>w)&1;   // edge e-w is the matching edge iff w is e's mate
                u32 m2=mask|(1u<<w);
                int8_t& slot=dp[(size_t)m2*N+w];
                if(cur+add>slot) slot=cur+add;
            }
        }
    }
    return best;
}
// plain DFS brute (independent of the DP) for validation
struct BruteMu {
    const G* g; int best;
    void run(int v,u32 vis,int mu){
        if(mu>best)best=mu;
        u32 cand=g->adj[v]&~vis;
        while(cand){ int w=__builtin_ctz(cand); cand&=cand-1;
            run(w,vis|(1u<<w),mu+((g->mate[v]>>w)&1)); }
    }
};
static int maxmu_brute(const G& g){
    BruteMu b{&g,0};
    for(int v=0;v<g.N;v++) b.run(v,1u<<v,0);
    return b.best;
}

// ---------- helpers ----------
static string show(const vector<int>& s){ string o="["; for(size_t i=0;i<s.size();i++){if(i)o+=",";o+=to_string(s[i]);} return o+"]"; }
static vector<int> parse_perm(const string& str){ vector<int> o; string c; for(char ch:str){ if(ch==','){o.push_back(stoi(c));c.clear();}else c+=ch;} if(!c.empty())o.push_back(stoi(c)); return o; }

// ---------- modes ----------

// checkdp <n> : DP vs brute over ALL sigma in S_n (validation; n<=7 practical)
static void checkdp(int n){
    long long viol=0,tot=0;
    vector<int> s(n); iota(s.begin(),s.end(),0);
    vector<int8_t> dp;
    do{
        G g=buildG(s);
        int a=maxmu_dp(g,dp), b=maxmu_brute(g); tot++;
        if(a!=b){ viol++; printf("DPVIOL s=%s dp=%d brute=%d\n",show(s).c_str(),a,b); if(viol>10)return; }
    }while(next_permutation(s.begin(),s.end()));
    printf("# checkdp n=%d tested=%lld violations=%lld\n",n,tot,viol);
}

// exh <n> : exhaustive sweep of S_n; report f(n), the maxmu histogram, the
// minimizers, and the CONVERSION GUARD (any sigma with maxmu < rho -> BUG).
// Parallel over first-symbol shards.
static void exh(int n){
    int fmin=INT_MAX; long long guard_viol=0;
    map<int,long long> hist;
    vector<pair<vector<int>,int>> minimizers;
    #pragma omp parallel
    {
        map<int,long long> lh;
        int lmin=INT_MAX; vector<pair<vector<int>,int>> lmins; long long lviol=0;
        vector<int8_t> dp;
        #pragma omp for schedule(dynamic)
        for(int a=0;a<n;a++){
            vector<int> rest; for(int i=0;i<n;i++) if(i!=a) rest.push_back(i);
            sort(rest.begin(),rest.end());
            do{
                vector<int> s; s.push_back(a); s.insert(s.end(),rest.begin(),rest.end());
                G g=buildG(s);
                int m=maxmu_dp(g,dp);
                lh[m]++;
                if(m<lmin){ lmin=m; lmins.clear(); }
                if(m==lmin && lmins.size()<50) lmins.push_back({s,m});
                int r=rho_perm(s);
                if(m<r) { lviol++; printf("GUARDVIOL s=%s maxmu=%d rho=%d\n",show(s).c_str(),m,r); }
            }while(next_permutation(rest.begin(),rest.end()));
        }
        #pragma omp critical
        { for(auto&kv:lh)hist[kv.first]+=kv.second;
          if(lmin<fmin){ fmin=lmin; minimizers=lmins; }
          else if(lmin==fmin) minimizers.insert(minimizers.end(),lmins.begin(),lmins.end());
          guard_viol+=lviol; }
    }
    long long tot=0; for(auto&kv:hist)tot+=kv.second;
    printf("# exh n=%d total=%lld guard_violations=%lld\n",n,tot,guard_viol);
    printf("# f(%d) = %d\n",n,fmin);
    printf("# maxmu HIST:"); for(auto&kv:hist) printf(" %d:%lld",kv.first,kv.second); printf("\n");
    long long mc=0; for(auto&kv:hist){ if(kv.first==fmin) mc=kv.second; }
    printf("# minimizer count = %lld (showing <=50)\n",mc);
    for(auto& p:minimizers) printf("MIN %d %s b=? rho=%d\n",p.second,show(p.first).c_str(),rho_perm(p.first));
}

// slack <n> : exhaustive joint (rho, maxmu) histogram over S_n (n<=8 practical).
static void slack(int n){
    map<pair<int,int>,long long> joint;
    long long guard_viol=0;
    #pragma omp parallel
    {
        map<pair<int,int>,long long> lj; long long lv=0;
        vector<int8_t> dp;
        #pragma omp for schedule(dynamic)
        for(int a=0;a<n;a++){
            vector<int> rest; for(int i=0;i<n;i++) if(i!=a) rest.push_back(i);
            sort(rest.begin(),rest.end());
            do{
                vector<int> s; s.push_back(a); s.insert(s.end(),rest.begin(),rest.end());
                G g=buildG(s);
                int m=maxmu_dp(g,dp);
                int r=rho_perm(s);
                lj[{r,m}]++;
                if(m<r) lv++;
            }while(next_permutation(rest.begin(),rest.end()));
        }
        #pragma omp critical
        { for(auto&kv:lj)joint[kv.first]+=kv.second; guard_viol+=lv; }
    }
    printf("# slack n=%d guard_violations=%lld\n",n,guard_viol);
    printf("# (rho, maxmu) : count ; slack = maxmu - rho\n");
    for(auto&kv:joint) printf("J %d %d %lld\n",kv.first.first,kv.first.second,kv.second);
}

// single <perm> : maxmu + rho + slack for one sigma (DP if 2n<=24, else brute DFS with progress)
static void single(const vector<int>& s){
    G g=buildG(s);
    int m;
    if(g.N<=24){ vector<int8_t> dp; m=maxmu_dp(g,dp); }
    else m=maxmu_brute(g);
    printf("n=%zu maxmu=%d rho=%d slack=%d\n",s.size(),m,rho_perm(s),m-rho_perm(s));
}

int main(int argc,char**argv){
    if(argc<2){fprintf(stderr,"mode: checkdp|exh|slack|single\n");return 1;}
    string m=argv[1];
    if(m=="checkdp") checkdp(atoi(argv[2]));
    else if(m=="exh") exh(atoi(argv[2]));
    else if(m=="slack") slack(atoi(argv[2]));
    else if(m=="single") single(parse_perm(argv[2]));
    else {fprintf(stderr,"unknown mode\n");return 1;}
    return 0;
}
