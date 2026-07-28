// k3.cpp — WP-J k=3 census machinery.
// 2-color core (build/reach_bound/rho/decide) copied VERBATIM from the WP-I
// census.cpp port of the paper's rho_tool.cpp (validated zero-mismatch there;
// re-validated here against rho_tool exh histograms — see dsweep).
//
// Conventions:
//  * H_sigma convention (2-color): vertices = positions, blue = {i,i+1},
//    red = {sigma^-1(v), sigma^-1(v+1)}.  rho counts VERTICES.
//  * word convention (k=3, per NGT-21/C8): for word w, P(w) has edges
//    {w[i], w[i+1]}.  Pair (id, w): f(w) = rho(id, w) = rho_pc(w^-1).
//    Triple (id, s, t): three colors E(id), E(s), E(t).
//  * D(n) = {w : f(w) <= n-1}.  Since rho_pc is inv-invariant (D4, validated),
//    D(n) = {sigma : rho_pc(sigma) <= n-1} as a set.
//
// Build: g++ -O2 -fopenmp -march=native -o k3 k3.cpp
#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef uint64_t u64;

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

struct DFS {
    const H* h; int best; long long nodes, budget; bool aborted;
    void run(int v,int c,u64 vis,int len,bool prune){
        if(aborted)return; if(++nodes>budget&&budget>0){aborted=true;return;}
        if(len>best)best=len;
        u64 cand=h->adj[c][v]&~vis; if(!cand)return;
        if(prune&&len+reach_bound(*h,v,c,vis)<=best)return;
        while(cand){int w=__builtin_ctzll(cand);cand&=cand-1;run(w,1-c,vis|(1ULL<<w),len+1,prune);} }
};
static int rho(const H& h, bool prune=false, long long budget=-1, bool* ok=nullptr){
    DFS d{&h,1,0,budget,false};
    for(int v=0;v<h.n;v++)for(int c=0;c<2;c++)d.run(v,c,1ULL<<v,1,prune);
    if(ok)*ok=!d.aborted; return d.best;
}
static int rho_perm(const vector<int>& s, bool prune=false){ return rho(build(s),prune); }

struct Dec {
    const H* h; int t; bool found;
    void run(int v,int c,u64 vis,int len){
        if(found)return;
        if(len>=t){found=true;return;}
        u64 cand=h->adj[c][v]&~vis; if(!cand)return;
        if(len+reach_bound(*h,v,c,vis)<t)return;
        while(cand){int w=__builtin_ctzll(cand);cand&=cand-1;
            run(w,1-c,vis|(1ULL<<w),len+1); if(found)return;} }
};
static bool decide(const H& h, int t){
    if(t<=1)return true;
    Dec d{&h,t,false};
    for(int v=0;v<h.n&&!d.found;v++)for(int c=0;c<2&&!d.found;c++)d.run(v,c,1ULL<<v,1);
    return d.found;
}
static bool decide_perm(const vector<int>& s,int t){ return decide(build(s),t); }

// ---------- k=3 : triple (id, s, t) in WORD convention ----------
struct H3 { int n; u64 adj[3][64]; };
static H3 build3(const vector<int>& s, const vector<int>& t){
    H3 h; int n=(int)s.size(); h.n=n;
    for(int c=0;c<3;c++)for(int v=0;v<n;v++)h.adj[c][v]=0;
    for(int i=0;i+1<n;i++){ h.adj[0][i]|=1ULL<<(i+1); h.adj[0][i+1]|=1ULL<<i; }
    for(int i=0;i+1<n;i++){ int a=s[i],b=s[i+1]; h.adj[1][a]|=1ULL<<b; h.adj[1][b]|=1ULL<<a; }
    for(int i=0;i+1<n;i++){ int a=t[i],b=t[i+1]; h.adj[2][a]|=1ULL<<b; h.adj[2][b]|=1ULL<<a; }
    return h;
}
// reachability bound for 3 colors: BFS over (vertex, arriving-color) states.
static int reach3(const H3& h, int v, int lastc, u64 vis){
    u64 seen[3]={0,0,0}, acc=0; // seen[c] = vertices reached arriving via color c
    u64 frontier[3]={0,0,0};
    for(int c=0;c<3;c++) if(c!=lastc){ u64 f=h.adj[c][v]&~vis; frontier[c]|=f; seen[c]|=f; acc|=f; }
    bool change=true;
    while(change){ change=false;
        for(int c=0;c<3;c++){ u64 f=frontier[c]; frontier[c]=0;
            while(f){ int w=__builtin_ctzll(f); f&=f-1;
                for(int c2=0;c2<3;c2++) if(c2!=c){
                    u64 nx=h.adj[c2][w]&~vis&~seen[c2];
                    if(nx){ seen[c2]|=nx; frontier[c2]|=nx; acc|=nx; change=true; } } } } }
    return __builtin_popcountll(acc);
}
struct Dec3 {
    const H3* h; int t; bool found;
    void run(int v,int lastc,u64 vis,int len){
        if(found)return;
        if(len>=t){found=true;return;}
        if(len+reach3(*h,v,lastc,vis)<t)return;
        for(int c=0;c<3;c++){ if(c==lastc)continue;
            u64 cand=h->adj[c][v]&~vis;
            while(cand){ int w=__builtin_ctzll(cand); cand&=cand-1;
                run(w,c,vis|(1ULL<<w),len+1); if(found)return; } } }
};
static bool decide3(const H3& h, int t){
    if(t<=1)return true;
    Dec3 d{&h,t,false};
    for(int v=0;v<h.n&&!d.found;v++) d.run(v,3,1ULL<<v,1); // lastc=3: any first color
    return d.found;
}
// exact rho3 (for spot checks / deficit measurement)
struct DFS3 {
    const H3* h; int best;
    void run(int v,int lastc,u64 vis,int len){
        if(len>best)best=len;
        if(len+reach3(*h,v,lastc,vis)<=best)return;
        for(int c=0;c<3;c++){ if(c==lastc)continue;
            u64 cand=h->adj[c][v]&~vis;
            while(cand){ int w=__builtin_ctzll(cand); cand&=cand-1;
                run(w,c,vis|(1ULL<<w),len+1); } } }
};
static int rho3(const H3& h){
    DFS3 d{&h,1};
    for(int v=0;v<h.n;v++) d.run(v,3,1ULL<<v,1);
    return d.best;
}

// ---------- perm packing (n<=16: 4 bits per entry) ----------
static u64 pack(const vector<int>& s){ u64 x=0; for(size_t i=0;i<s.size();i++) x |= (u64)s[i]<<(4*i); return x; }
static vector<int> unpack(u64 x,int n){ vector<int> s(n); for(int i=0;i<n;i++) s[i]=(int)((x>>(4*i))&0xF); return s; }
static string show(const vector<int>& s){ string o="["; for(size_t i=0;i<s.size();i++){if(i)o+=",";o+=to_string(s[i]);} return o+"]"; }
static vector<int> parse_perm(const string& str){ vector<int> o; string c; for(char ch:str){ if(ch==','){o.push_back(stoi(c));c.clear();}else c+=ch;} if(!c.empty())o.push_back(stoi(c)); return o; }

// ---------- modes ----------

// dsweep <n> : sweep S_n; dump all sigma with rho_pc(sigma) <= n-1 (the deficient
// set D(n)), one per line "rho perm", plus the FULL rho histogram (cross-check vs
// rho_tool exh HIST).  Parallel over first-symbol shards.
static void dsweep(int n){
    vector<long long> hist(n+1,0);
    vector<pair<int,u64>> defic;
    #pragma omp parallel
    {
        vector<long long> lh(n+1,0);
        vector<pair<int,u64>> ld;
        #pragma omp for schedule(dynamic)
        for(int a=0;a<n;a++){
            vector<int> rest; for(int i=0;i<n;i++) if(i!=a) rest.push_back(i);
            sort(rest.begin(),rest.end());
            do {
                vector<int> s; s.push_back(a); s.insert(s.end(),rest.begin(),rest.end());
                if(decide_perm(s,n)) { lh[n]++; continue; }   // PC-Ham
                int r=rho_perm(s); lh[r]++;
                ld.push_back({r,pack(s)});
            } while(next_permutation(rest.begin(),rest.end()));
        }
        #pragma omp critical
        { for(int i=0;i<=n;i++)hist[i]+=lh[i]; defic.insert(defic.end(),ld.begin(),ld.end()); }
    }
    long long tot=0; for(auto v:hist)tot+=v;
    fprintf(stderr,"dsweep n=%d total=%lld |D|=%zu\n",n,tot,defic.size());
    printf("# dsweep n=%d total=%lld\n",n,tot);
    printf("# HIST:"); for(int r=1;r<=n;r++) if(hist[r]) printf(" %d:%lld",r,hist[r]); printf("\n");
    printf("# |D(%d)| = %zu\n",n,defic.size());
    sort(defic.begin(),defic.end());
    for(auto& p:defic) printf("%d %s\n",p.first,show(unpack(p.second,n)).c_str());
}

// tri <n> <dfile> : read D(n) ("rho perm" lines); check all ordered pairs (s,t)
// for s^-1 t in D; print every coset triangle {s,t,s^-1t}.
static void tri(int n, const char* dfile){
    ifstream in(dfile); string line; vector<vector<int>> D;
    while(getline(in,line)){
        if(line.empty()||line[0]=='#')continue;
        size_t sp=line.find(' ');
        string ps=line.substr(sp+1);
        ps=ps.substr(1,ps.size()-2); // strip [ ]
        D.push_back(parse_perm(ps));
    }
    unordered_set<u64> dset; for(auto& s:D) dset.insert(pack(s));
    printf("# tri n=%d |D|=%zu ordered-pairs=%lld\n",n,D.size(),(long long)D.size()*(long long)D.size());
    long long tricount=0;
    #pragma omp parallel for schedule(dynamic) reduction(+:tricount)
    for(int i=0;i<(int)D.size();i++){
        const vector<int>& s=D[i];
        vector<int> sinv(n); for(int k2=0;k2<n;k2++) sinv[s[k2]]=k2;
        for(int j=0;j<(int)D.size();j++){
            const vector<int>& t=D[j];
            vector<int> u(n); for(int k2=0;k2<n;k2++) u[k2]=sinv[t[k2]]; // s^-1 t
            if(dset.count(pack(u))){
                tricount++;
                #pragma omp critical
                printf("TRIANGLE s=%s t=%s sinvt=%s\n",show(s).c_str(),show(t).c_str(),show(u).c_str());
            }
        }
    }
    printf("# triangles(ordered pairs with s^-1 t in D) = %lld\n",tricount);
}

// trifull <n> <dfile> : for every ordered pair (s,t) in D x D with i<=j and
// s^-1 t in D (a coset triangle), run decide3(id,s,t,n) INLINE; print only
// deficits.  Counts: triangle pairs, checked triples, deficits.
static void trifull(int n, const char* dfile){
    ifstream in(dfile); string line; vector<vector<int>> D;
    while(getline(in,line)){
        if(line.empty()||line[0]=='#')continue;
        size_t sp=line.find(' ');
        string ps=line.substr(sp+1); ps=ps.substr(1,ps.size()-2);
        D.push_back(parse_perm(ps));
    }
    unordered_set<u64> dset; dset.reserve(D.size()*2);
    for(auto& s:D) dset.insert(pack(s));
    long long ND=D.size();
    fprintf(stderr,"trifull n=%d |D|=%lld pairs(i<=j)=%lld\n",n,ND,ND*(ND+1)/2);
    long long tricnt=0, defcnt=0;
    #pragma omp parallel for schedule(dynamic,16) reduction(+:tricnt,defcnt)
    for(long long i=0;i<ND;i++){
        const vector<int>& s=D[i];
        vector<int> sinv(n); for(int k2=0;k2<n;k2++) sinv[s[k2]]=k2;
        vector<int> u(n);
        for(long long j=i;j<ND;j++){
            const vector<int>& t=D[j];
            for(int k2=0;k2<n;k2++) u[k2]=sinv[t[k2]];
            if(dset.count(pack(u))){
                tricnt++;
                H3 h=build3(s,t);
                if(!decide3(h,n)){ defcnt++;
                    #pragma omp critical
                    { printf("DEF3 s=%s t=%s rho3=%d\n",show(s).c_str(),show(t).c_str(),rho3(h)); fflush(stdout); }
                }
            }
        }
    }
    printf("# trifull n=%d |D|=%lld triangle-pairs(i<=j)=%lld checked=%lld deficits=%lld\n",
           n,ND,tricnt,tricnt,defcnt);
}

// trisweep <n> : EXHAUSTIVE all ordered pairs (s,t) in S_n x S_n, 3-color decide
// at threshold n; count and print deficient triples.  (feasible n<=7)
static void trisweep(int n){
    long long defc=0, tot=0;
    vector<u64> perms;
    { vector<int> s(n); iota(s.begin(),s.end(),0);
      do{ perms.push_back(pack(s)); }while(next_permutation(s.begin(),s.end())); }
    long long P=perms.size();
    #pragma omp parallel for schedule(dynamic,64) reduction(+:defc,tot)
    for(long long i=0;i<P;i++){
        vector<int> s=unpack(perms[i],n);
        for(long long j=0;j<P;j++){
            vector<int> t=unpack(perms[j],n);
            tot++;
            H3 h=build3(s,t);
            if(!decide3(h,n)){ defc++;
                #pragma omp critical
                printf("DEF3 s=%s t=%s rho3=%d\n",show(s).c_str(),show(t).c_str(),rho3(h));
            }
        }
    }
    printf("# trisweep n=%d pairs=%lld deficient=%lld\n",n,tot,defc);
}

// tricheck <n> <dfile> : for every coset triangle found in D (or if none, exit),
// run the direct 3-color decision on (id,s,t).  Reads TRIANGLE lines from stdin.
static void tricheck(int n){
    string line;
    while(getline(cin,line)){
        if(line.rfind("TRIANGLE",0)!=0)continue;
        // parse s=[..] t=[..]
        size_t i1=line.find("s=["), i2=line.find("] t=["), i3=line.find("] sinvt=");
        vector<int> s=parse_perm(line.substr(i1+3,i2-i1-3));
        vector<int> t=parse_perm(line.substr(i2+5,i3-i2-5));
        H3 h=build3(s,t);
        bool ham=decide3(h,n);
        printf("TRIPLE s=%s t=%s ham3=%d%s\n",show(s).c_str(),show(t).c_str(),(int)ham,
               ham?"":" *** k=3 DEFICIT ***");
        if(!ham) printf("   rho3=%d\n",rho3(h));
    }
}

// randtri <n> <cnt> <seed> : random triples; verify B1: rho3(id,s,t) >= max(f(s),f(t),f(s^-1 t))
// where f(w)=rho(id,w)=rho_pc(w^-1). Engine cross-validation (2c vs 3c).
static void randtri(int n, long long cnt, u64 seed){
    mt19937_64 rng(seed);
    long long viol=0;
    for(long long it=0; it<cnt; it++){
        vector<int> s(n),t(n); iota(s.begin(),s.end(),0); iota(t.begin(),t.end(),0);
        shuffle(s.begin(),s.end(),rng); shuffle(t.begin(),t.end(),rng);
        vector<int> sinv(n),u(n); for(int i2=0;i2<n;i2++) sinv[s[i2]]=i2;
        for(int i2=0;i2<n;i2++) u[i2]=sinv[t[i2]];
        // f(w) = rho_pc(w^-1); rho_pc computed by 2-color engine on w^-1 (H convention)
        auto f=[&](const vector<int>& w){ vector<int> wi(n); for(int i3=0;i3<n;i3++) wi[w[i3]]=i3; return rho_perm(wi); };
        int fs=f(s), ft=f(t), fu=f(u);
        int r3=rho3(build3(s,t));
        if(r3 < max({fs,ft,fu})){ viol++;
            printf("B1VIOL s=%s t=%s rho3=%d fs=%d ft=%d fu=%d\n",show(s).c_str(),show(t).c_str(),r3,fs,ft,fu);
        }
    }
    printf("# randtri n=%d cnt=%lld B1-violations=%lld\n",n,cnt,viol);
}

int main(int argc,char**argv){
    if(argc<2){fprintf(stderr,"mode required: dsweep|tri|trisweep|tricheck|randtri|f|rho3\n");return 1;}
    string m=argv[1];
    if(m=="dsweep") dsweep(atoi(argv[2]));
    else if(m=="tri") tri(atoi(argv[2]),argv[3]);
    else if(m=="trifull") trifull(atoi(argv[2]),argv[3]);
    else if(m=="trisweep") trisweep(atoi(argv[2]));
    else if(m=="tricheck") tricheck(atoi(argv[2]));
    else if(m=="randtri") randtri(atoi(argv[2]),atoll(argv[3]),argc>4?atoll(argv[4]):42);
    else if(m=="f"){ vector<int> w=parse_perm(argv[2]); int n=(int)w.size();
        vector<int> wi(n); for(int i=0;i<n;i++) wi[w[i]]=i;
        printf("f=%d\n",rho_perm(wi)); }
    else if(m=="rho3"){ vector<int> s=parse_perm(argv[2]), t=parse_perm(argv[3]);
        printf("rho3=%d\n",rho3(build3(s,t))); }
    else {fprintf(stderr,"unknown mode\n");return 1;}
    return 0;
}
