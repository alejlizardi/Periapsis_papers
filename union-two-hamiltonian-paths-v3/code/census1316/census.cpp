// census.cpp — WP-I rho_min census engine (n=13,14).
// Ports rho() + reach_bound() VERBATIM from the paper artifact rho_tool.cpp
// (same H_sigma conventions: vertices = positions, blue = position edges,
// red = value edges; rho counts VERTICES on an alternating PC path).
// Adds: (1) threshold decision decide(sigma,t) with early exit + reach prune,
//       (2) D4 symmetry canonicalization (rev / comp / inv generators),
//       (3) validation / seed / census modes with checkpoint-resume.
//
// Build: g++ -O2 -fopenmp -march=native -o census census.cpp
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

// --- reach_bound: alternation-aware reachable-vertex upper bound (VERBATIM). ---
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

// --- exact rho (VERBATIM from rho_tool, prune optional). ---
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

// --- threshold decision: is there a PC path on >= t vertices? early exit. ---
struct Dec {
    const H* h; int t; bool found;
    void run(int v,int c,u64 vis,int len){
        if(found)return;
        if(len>=t){found=true;return;}
        u64 cand=h->adj[c][v]&~vis; if(!cand)return;
        if(len+reach_bound(*h,v,c,vis)<t)return;   // cannot possibly reach t
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

// ---------- blocks (for b(sigma) reporting) ----------
static int nblocks(const vector<int>& s){
    int n=(int)s.size(),i=0,b=0;
    while(i<n){int j=i; if(j+1<n&&abs(s[j+1]-s[j])==1){int eps=s[j+1]-s[j];while(j+1<n&&s[j+1]-s[j]==eps)j++;} b++; i=j+1;}
    return b;
}

// ---------- D4 symmetry ----------
static void rev_(const int* s,int n,int* o){ for(int i=0;i<n;i++)o[i]=s[n-1-i]; }
static void comp_(const int* s,int n,int* o){ for(int i=0;i<n;i++)o[i]=n-1-s[i]; }
static void inv_(const int* s,int n,int* o){ for(int i=0;i<n;i++)o[s[i]]=i; }

// fill up to 8 D4 images of s into img[8][n]; returns them (with possible dupes).
static void d4_images(const int* s,int n,int img[8][16]){
    int a[4][16];
    for(int i=0;i<n;i++)a[0][i]=s[i];
    rev_(a[0],n,a[1]); comp_(a[0],n,a[2]); rev_(a[2],n,a[3]); // e, rev, comp, rev*comp
    for(int k=0;k<4;k++){ for(int i=0;i<n;i++)img[k][i]=a[k][i]; inv_(a[k],n,img[4+k]); }
}
static int cmp_(const int* x,const int* y,int n){ for(int i=0;i<n;i++){if(x[i]<y[i])return -1;if(x[i]>y[i])return 1;} return 0; }

// is s the lexicographically smallest in its D4 orbit?
static bool is_canon(const int* s,int n){
    int img[8][16]; d4_images(s,n,img);
    for(int k=1;k<8;k++) if(cmp_(img[k],s,n)<0) return false;
    return true;
}
// orbit size (distinct images) of s under D4.
static int orbit_size(const int* s,int n){
    int img[8][16]; d4_images(s,n,img); int cnt=0;
    for(int k=0;k<8;k++){ bool dup=false; for(int j=0;j<k;j++) if(cmp_(img[k],img[j],n)==0){dup=true;break;} if(!dup)cnt++; }
    return cnt;
}

static string show(const vector<int>& s){ string o="["; for(size_t i=0;i<s.size();i++){if(i)o+=",";o+=to_string(s[i]);} return o+"]"; }
static vector<int> parse_perm(const string& str){ vector<int> o; string c; for(char ch:str){ if(ch==','){o.push_back(stoi(c));c.clear();}else c+=ch;} if(!c.empty())o.push_back(stoi(c)); return o; }

// ---------- modes ----------

// dumprho <n> : print exact rho for every sigma in S_n  (cross-check vs pc.py)
static void dumprho(int n){
    vector<int> s(n); iota(s.begin(),s.end(),0);
    do { printf("%s %d\n", ([&]{string o;for(int i=0;i<n;i++){if(i)o+=",";o+=to_string(s[i]);}return o;}()).c_str(), rho_perm(s)); }
    while(next_permutation(s.begin(),s.end()));
}

// checksym <n> : verify rho invariance under rev,comp,inv over all S_n (self, fast)
static void checksym(int n){
    long long viol=0, tot=0;
    vector<int> s(n); iota(s.begin(),s.end(),0);
    do {
        int r=rho_perm(s); tot++;
        int a[4][16], img[8][16]; d4_images(s.data(),n,img);
        for(int k=1;k<8;k++){ vector<int> t(img[k],img[k]+n); if(rho_perm(t)!=r){ viol++;
            printf("SYMVIOL k=%d s=%s r=%d rk=%d\n",k,show(s).c_str(),r,rho_perm(t)); if(viol>20)return; } }
        (void)a;
    } while(next_permutation(s.begin(),s.end()));
    printf("checksym n=%d tested=%lld violations=%lld\n",n,tot,viol);
}

// singlefile <n> : read comma-perms from stdin, print "rho b orbit canon" each (batch resolve)
static void singlefile(int n){
    string line;
    while(getline(cin,line)){
        if(line.empty())continue;
        vector<int> s=parse_perm(line);
        int nn=(int)s.size();
        printf("%d %d %d %d\n", rho_perm(s, nn>=12), nblocks(s), orbit_size(s.data(),nn), (int)is_canon(s.data(),nn));
    }
}

// rhofile <n> : read comma-perms from stdin (one/line), print exact rho each (batch ref)
static void rhofile(int n){
    string line;
    while(getline(cin,line)){
        if(line.empty())continue;
        vector<int> s=parse_perm(line);
        printf("%d\n", rho_perm(s, n>=12));
    }
}

// checkprune <n> : verify rho(prune=false)==rho(prune=true) for all s in S_n
static void checkprune(int n){
    long long viol=0,tot=0;
    vector<int> s(n); iota(s.begin(),s.end(),0);
    do { H h=build(s); int r0=rho(h,false), r1=rho(h,true); tot++;
        if(r0!=r1){viol++; printf("PRUNEVIOL s=%s unpruned=%d pruned=%d\n",show(s).c_str(),r0,r1); if(viol>20)return;}
    } while(next_permutation(s.begin(),s.end()));
    printf("checkprune n=%d tested=%lld violations=%lld\n",n,tot,viol);
}

// checkdecide <n> : verify decide(s,t) == (rho(s)>=t) for all s in S_n, all t in [1..n]
static void checkdecide(int n){
    long long viol=0, tot=0;
    vector<int> s(n); iota(s.begin(),s.end(),0);
    do { int r=rho_perm(s);
        for(int t=1;t<=n;t++){ bool d=decide_perm(s,t); bool exp=(r>=t); tot++;
            if(d!=exp){ viol++; printf("DECVIOL s=%s t=%d rho=%d decide=%d\n",show(s).c_str(),t,r,(int)d); if(viol>20)return; } }
    } while(next_permutation(s.begin(),s.end()));
    printf("checkdecide n=%d tested=%lld violations=%lld\n",n,tot,viol);
}

// randcheck <n> <samples> <seed> : decide(s,t) vs rho(s) on random s, all t; zero tol
static void randcheck(int n,long long samples,u64 seed){
    long long viol=0,tot=0;
#pragma omp parallel reduction(+:viol,tot)
    {
#ifdef _OPENMP
        mt19937_64 rng(seed+1000*omp_get_thread_num());
#else
        mt19937_64 rng(seed);
#endif
#pragma omp for schedule(dynamic,256)
        for(long long i=0;i<samples;i++){
            vector<int> s(n); iota(s.begin(),s.end(),0);
            for(int k=n-1;k>0;k--){int j=(int)(rng()%(k+1));swap(s[k],s[j]);}
            H h=build(s); int r=rho(h, n>=12);   // pruned (exact, admissible) — feasible at n=14
            for(int t=1;t<=n;t++){ bool d=decide(h,t); bool e=(r>=t); tot++;
                if(d!=e){ viol++;
#pragma omp critical
                    printf("RANDDECVIOL s=%s t=%d rho=%d decide=%d\n",show(s).c_str(),t,r,(int)d);} }
        }
    }
    printf("randcheck n=%d tested_perms=%lld checks=%lld violations=%lld\n",n,samples,tot,viol);
}

// seed <n> <samples> <seed> : find low rho via random + random-block perms
static vector<int> random_block_perm(int n, mt19937_64& rng){
    vector<int> sizes; int rem=n;
    while(rem>0){int s=1+(int)(rng()%4); if(s>rem)s=rem; sizes.push_back(s); rem-=s;}
    int b=(int)sizes.size(); vector<int> order(b); iota(order.begin(),order.end(),0);
    shuffle(order.begin(),order.end(),rng);
    vector<int> rank_(b); for(int r=0;r<b;r++)rank_[order[r]]=r;
    vector<int> base_by_rank(b); {int a2=0;for(int r=0;r<b;r++){base_by_rank[r]=a2;a2+=sizes[order[r]];}}
    vector<int> s;
    for(int k=0;k<b;k++){int sz=sizes[k],lo=base_by_rank[rank_[k]];bool desc=rng()&1;
        if(sz==1)s.push_back(lo); else if(desc)for(int v=lo+sz-1;v>=lo;v--)s.push_back(v); else for(int v=lo;v<lo+sz;v++)s.push_back(v);}
    return s;
}
static void seedscan(int n,long long samples,u64 seed){
    int gmin=INT_MAX; vector<int> best;
#pragma omp parallel
    {
#ifdef _OPENMP
        mt19937_64 rng(seed+7919*omp_get_thread_num());
#else
        mt19937_64 rng(seed);
#endif
        int lmin=INT_MAX; vector<int> lbest;
#pragma omp for schedule(dynamic,64)
        for(long long i=0;i<samples;i++){
            vector<int> s = (i&1)? random_block_perm(n,rng) : [&]{vector<int> t(n);iota(t.begin(),t.end(),0);for(int k=n-1;k>0;k--){int j=(int)(rng()%(k+1));swap(t[k],t[j]);}return t;}();
            int r=rho_perm(s,true);
            if(r<lmin){lmin=r;lbest=s;}
        }
#pragma omp critical
        if(lmin<gmin){gmin=lmin;best=lbest;}
    }
    printf("seed n=%d samples=%lld min_rho=%d b=%d sigma=%s\n",n,samples,gmin,nblocks(best),show(best).c_str());
}

// ---------- WP-K counting modes ----------
// countfull <n> : |D(n)| with NO symmetry — every sigma in S_n, deficient iff decide(sigma,n)==false.
// Also splits deficit: deep = #{rho <= n-2} = decide(sigma,n-1)==false (only evaluated on deficient).
static void countfull(int n){
    long long def=0, deep=0, tot=0;
    double t0=0;
#ifdef _OPENMP
    t0=omp_get_wtime();
#endif
#pragma omp parallel for schedule(dynamic,1) reduction(+:def,deep,tot)
    for(int a=0;a<n;a++){
        vector<int> rest; for(int x=0;x<n;x++) if(x!=a) rest.push_back(x);
        sort(rest.begin(),rest.end());
        vector<int> s(n); s[0]=a;
        do{
            for(int i=1;i<n;i++)s[i]=rest[i-1];
            tot++;
            H h=build(s);
            if(!decide(h,n)){ def++; if(!decide(h,n-1)) deep++; }
        }while(next_permutation(rest.begin(),rest.end()));
    }
#ifdef _OPENMP
    fprintf(stderr,"[countfull] wall=%.2fs\n",omp_get_wtime()-t0);
#endif
    printf("countfull n=%d total=%lld deficient=%lld deep=%lld\n",n,tot,def,deep);
}

// countsym <n> <c0lo> <c0hi> : orbit-weighted |D(n)| over canonical reps only.
// For each canonical sigma with decide(sigma,n)==false, add its TRUE orbit size (8/|stab|,
// computed by applying all 8 group elements and counting distinct images) — trap 1: counts
// are weighted by orbit size, never reps*8. deep split (t=n-1) fused in.
// Emits "PROG chunks_done chunks_total elapsed_s" on stderr every few chunks (heartbeat feed).
static void countsym(int n,int c0lo,int c0hi){
    vector<array<int,3>> chunks;
    for(int a=c0lo;a<c0hi;a++)for(int b=0;b<n;b++){if(b==a)continue;for(int c=0;c<n;c++){if(c==a||c==b)continue;chunks.push_back({a,b,c});}}
    long long scanned=0, canon=0, defreps=0, defw=0, deepreps=0, deepw=0;
    long long chunks_done=0;
    double t0=0;
#ifdef _OPENMP
    t0=omp_get_wtime();
#endif
#pragma omp parallel reduction(+:scanned,canon,defreps,defw,deepreps,deepw)
    {
#pragma omp for schedule(dynamic,1)
        for(long long ci=0;ci<(long long)chunks.size();ci++){
            int a=chunks[ci][0],b=chunks[ci][1],c=chunks[ci][2];
            vector<int> rest; for(int x=0;x<n;x++)if(x!=a&&x!=b&&x!=c)rest.push_back(x);
            sort(rest.begin(),rest.end());
            vector<int> s(n); s[0]=a;s[1]=b;s[2]=c;
            do {
                for(int i=3;i<n;i++)s[i]=rest[i-3];
                scanned++;
                if(!is_canon(s.data(),n)) continue;
                canon++;
                H h=build(s);
                if(decide(h,n)) continue;               // has a PC-Ham path: not deficient
                int ob=orbit_size(s.data(),n);
                defreps++; defw+=ob;
                if(!decide(h,n-1)){ deepreps++; deepw+=ob; }
            } while(next_permutation(rest.begin(),rest.end()));
#pragma omp critical
            {
                chunks_done++;
                if(chunks_done%16==0||chunks_done==(long long)chunks.size()){
#ifdef _OPENMP
                    fprintf(stderr,"PROG %lld %zu %.1f\n",chunks_done,chunks.size(),omp_get_wtime()-t0);
#else
                    fprintf(stderr,"PROG %lld %zu 0\n",chunks_done,chunks.size());
#endif
                    fflush(stderr);
                }
            }
        }
    }
#ifdef _OPENMP
    fprintf(stderr,"[countsym] wall=%.2fs\n",omp_get_wtime()-t0);
#endif
    printf("countsym n=%d c0=[%d,%d) scanned=%lld canon=%lld def_reps=%lld def_weighted=%lld deep_reps=%lld deep_weighted=%lld\n",
           n,c0lo,c0hi,scanned,canon,defreps,defw,deepreps,deepw);
}

// census <n> <K> <c0lo> <c0hi> [outfile]
//   collect ALL canonical sigma with rho(sigma) <= K, i.e. decide(sigma,K+1)==false.
//   parallel over (s[0],s[1],s[2]) chunks; only chunks with c0 in [c0lo,c0hi) run
//   (enables coarse checkpoint/resume by first-symbol ranges).
//   Prints survivors as: SURV rho b orbit sigma
static void census(int n,int K,int c0lo,int c0hi){
    // enumerate all perms with s[0]=a in [c0lo,c0hi); parallel over (a,b,c) triples
    vector<array<int,3>> chunks;
    for(int a=c0lo;a<c0hi;a++)for(int b=0;b<n;b++){if(b==a)continue;for(int c=0;c<n;c++){if(c==a||c==b)continue;chunks.push_back({a,b,c});}}
    long long survtot=0;
    int decT=K+1;
    vector<string> survlines;
#pragma omp parallel
    {
        vector<string> local;
        long long lsurv=0;
#pragma omp for schedule(dynamic,1)
        for(long long ci=0;ci<(long long)chunks.size();ci++){
            int a=chunks[ci][0],b=chunks[ci][1],c=chunks[ci][2];
            vector<int> rest; for(int x=0;x<n;x++)if(x!=a&&x!=b&&x!=c)rest.push_back(x);
            sort(rest.begin(),rest.end());
            vector<int> s(n); s[0]=a;s[1]=b;s[2]=c;
            do {
                for(int i=3;i<n;i++)s[i]=rest[i-3];
                if(!is_canon(s.data(),n)) continue;
                H h=build(s);
                if(decide(h,decT)) continue;            // rho>=K+1, not a survivor
                int r=rho(h,true);                       // exact (rho<=K)
                int ob=orbit_size(s.data(),n);
                char buf[256];
                snprintf(buf,sizeof buf,"SURV %d %d %d %s",r,nblocks(s),ob,show(s).c_str());
                local.push_back(buf); lsurv++;
            } while(next_permutation(rest.begin(),rest.end()));
        }
#pragma omp critical
        { for(auto& l:local)survlines.push_back(l); survtot+=lsurv; }
    }
    for(auto& l:survlines) printf("%s\n",l.c_str());
    // summarize
    int mn=INT_MAX; for(auto& l:survlines){int r=atoi(l.c_str()+5); if(r<mn)mn=r;}
    long long witn=0,witn_full=0;
    for(auto& l:survlines){int r,bb,ob;char t[8];sscanf(l.c_str(),"%s %d %d %d",t,&r,&bb,&ob); if(r==mn){witn++;witn_full+=ob;}}
    printf("census n=%d K=%d c0=[%d,%d) canon_survivors=%lld  rho_min_here=%d  canon_witnesses=%lld full_witnesses=%lld\n",
           n,K,c0lo,c0hi,survtot,mn,witn,witn_full);
}

int main(int argc,char**argv){
    if(argc<2){fprintf(stderr,"mode required\n");return 1;}
    string m=argv[1];
    if(m=="single"){vector<int> s=parse_perm(argv[2]);printf("n=%zu b=%d rho=%d orbit=%d canon=%d\n",s.size(),nblocks(s),rho_perm(s,s.size()>13),orbit_size(s.data(),(int)s.size()),(int)is_canon(s.data(),(int)s.size()));}
    else if(m=="decide"){vector<int> s=parse_perm(argv[2]);int t=atoi(argv[3]);printf("decide=%d rho=%d\n",(int)decide_perm(s,t),rho_perm(s,s.size()>13));}
    else if(m=="dumprho")dumprho(atoi(argv[2]));
    else if(m=="checksym")checksym(atoi(argv[2]));
    else if(m=="rhofile")rhofile(atoi(argv[2]));
    else if(m=="singlefile")singlefile(atoi(argv[2]));
    else if(m=="checkprune")checkprune(atoi(argv[2]));
    else if(m=="checkdecide")checkdecide(atoi(argv[2]));
    else if(m=="randcheck")randcheck(atoi(argv[2]),atoll(argv[3]),argc>4?atoll(argv[4]):12345);
    else if(m=="seed")seedscan(atoi(argv[2]),atoll(argv[3]),argc>4?atoll(argv[4]):999);
    else if(m=="census"){int n=atoi(argv[2]),K=atoi(argv[3]);int lo=argc>4?atoi(argv[4]):0,hi=argc>5?atoi(argv[5]):n;census(n,K,lo,hi);}
    else if(m=="countfull")countfull(atoi(argv[2]));
    else if(m=="countsym"){int n=atoi(argv[2]);int lo=argc>3?atoi(argv[3]):0,hi=argc>4?atoi(argv[4]):n;countsym(n,lo,hi);}
    else {fprintf(stderr,"unknown mode\n");return 1;}
    return 0;
}
