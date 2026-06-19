// search.cpp -- prescribed-automorphism search for a recursively differentiable
// quasigroup of order n (the construction method of the paper, Section 4).
//
// It looks for a Latin square L such that the derived table D[a][b]=L[b][L[a][b]] is
// also Latin, constrained to be invariant under a prescribed automorphism g. The search
// is a backtracker over g-orbit representatives with most-constrained-variable ordering,
// incremental pruning on the partial derivative table, and randomized restarts.
//
// This program is provided for reproducibility only: it is how the two tables were
// FOUND. Their correctness does not depend on it -- the printed tables are certified by
// the O(n^2) checks in rdiff_verify.py. To reproduce the tables in the paper:
//   order 14:  ./search 14 "cyc:0,1,2,3,4,5,6|7,8,9,10,11,12,13" 600000 1
//   order 18:  ./search 18 "cyc:0,1,2,3,4,5,6,7,8|9,10,11,12,13,14,15,16,17" 3600000 1
// (the exact table returned depends on the seed; any solution found is a valid witness).
//
// Build: g++ -O3 -march=native -std=c++17 search.cpp -o search
// Args:  search <n> <mode> <time_ms> <seed> [idem] [restart_base_nodes]
//   mode: raw | invol | ncyc | cyc:a,b,...|c,d,...
//
// Output (stderr): one status line. On FOUND, prints L's rows (stdout) and writes JSON.

#include <bits/stdc++.h>
using namespace std;

static int N, FULL;
static vector<vector<int>> L;
static vector<int> rowUsed, colUsed, dRowUsed, dColUsed;
static vector<vector<int>> Dval;
static vector<int> g; static bool haveAut; static int autOrder;
static long long deadline_ms; static chrono::steady_clock::time_point startTime;
static bool timedOut=false; static long long nodes=0; static bool solved=false;
static std::mt19937 rng;

static inline long long now_ms(){ return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-startTime).count(); }

static bool setD(int a,int b,int dv,vector<tuple<int,int,int>>&applied){
    if(Dval[a][b]!=-1) return Dval[a][b]==dv;
    int bit=1<<dv;
    if(dRowUsed[a]&bit) return false;
    if(dColUsed[b]&bit) return false;
    Dval[a][b]=dv; dRowUsed[a]|=bit; dColUsed[b]|=bit; applied.emplace_back(a,b,dv); return true;
}
static bool assign(int r,int c,int val,vector<tuple<int,int,int>>&applied){
    int bit=1<<val;
    if(rowUsed[r]&bit) return false;
    if(colUsed[c]&bit) return false;
    L[r][c]=val; rowUsed[r]|=bit; colUsed[c]|=bit;
    if(L[c][val]!=-1){ if(!setD(r,c,L[c][val],applied)) return false; }
    for(int a=0;a<N;++a) if(L[a][r]==c){ if(!setD(a,r,val,applied)) return false; }
    return true;
}
static void unassign(int r,int c,vector<tuple<int,int,int>>&applied){
    for(auto&[a,b,dv]:applied){ int bit=1<<dv; Dval[a][b]=-1; dRowUsed[a]&=~bit; dColUsed[b]&=~bit; }
    applied.clear();
    int v=L[r][c],bit=1<<v; rowUsed[r]&=~bit; colUsed[c]&=~bit; L[r][c]=-1;
}
struct CellSet{int r,c,val;vector<tuple<int,int,int>>applied;};
static bool assignOrbit(int r,int c,int val,vector<CellSet>&log){
    if(!haveAut){ CellSet cs;cs.r=r;cs.c=c;cs.val=val;
        if(L[r][c]!=-1) return L[r][c]==val;
        if(!assign(r,c,val,cs.applied)){unassign(r,c,cs.applied);return false;}
        log.push_back(std::move(cs)); return true; }
    int cr=r,cc=c,cv=val;
    for(int step=0;step<autOrder;++step){
        if(L[cr][cc]!=-1){ if(L[cr][cc]!=cv) return false; }
        else{ CellSet cs;cs.r=cr;cs.c=cc;cs.val=cv;
            if(!assign(cr,cc,cv,cs.applied)){unassign(cr,cc,cs.applied);return false;}
            log.push_back(std::move(cs)); }
        cr=g[cr];cc=g[cc];cv=g[cv];
        if(cr==r&&cc==c&&cv==val) break;
    }
    return true;
}
static void undoLog(vector<CellSet>&log,size_t from){ for(size_t i=log.size();i-->from;) unassign(log[i].r,log[i].c,log[i].applied); log.resize(from); }

static vector<pair<int,int>> repCells;
static void buildReps(){
    repCells.clear(); vector<vector<char>> seen(N,vector<char>(N,0));
    for(int r=0;r<N;++r)for(int c=0;c<N;++c){ if(seen[r][c])continue;
        if(haveAut){int cr=r,cc=c;for(int s=0;s<autOrder;++s){seen[cr][cc]=1;cr=g[cr];cc=g[cc];if(cr==r&&cc==c)break;}}
        else seen[r][c]=1;
        repCells.emplace_back(r,c); }
}

static long long restartLimit; static bool restartHit=false;

static bool backtrack(vector<CellSet>&log){
    if(solved) return true;
    if((++nodes & 0x3FFF)==0){ if(now_ms()>deadline_ms){timedOut=true;return false;} }
    if(timedOut) return false;
    if(nodes>restartLimit){ restartHit=true; return false; }   // trigger restart

    // MRV with random tie-break
    int bestCnt=1<<30; vector<pair<int,int>> bestCells;
    for(auto&pr:repCells){ int r=pr.first,c=pr.second; if(L[r][c]!=-1)continue;
        int avail=FULL&~rowUsed[r]&~colUsed[c]; int cnt=__builtin_popcount((unsigned)avail);
        if(cnt==0) return false;
        if(cnt<bestCnt){bestCnt=cnt;bestCells.clear();bestCells.push_back(pr);}
        else if(cnt==bestCnt) bestCells.push_back(pr);
    }
    if(bestCells.empty()){ for(int r=0;r<N;++r)for(int c=0;c<N;++c)if(L[r][c]==-1)return false; solved=true; return true; }
    auto pr = bestCells[ rng() % bestCells.size() ];
    int r=pr.first,c=pr.second;
    int avail=FULL&~rowUsed[r]&~colUsed[c];
    // randomized value order
    vector<int> vals; while(avail){int v=__builtin_ctz(avail);avail&=avail-1;vals.push_back(v);}
    shuffle(vals.begin(),vals.end(),rng);
    for(int val:vals){ size_t mark=log.size();
        if(assignOrbit(r,c,val,log)){ if(backtrack(log)) return true; }
        undoLog(log,mark);
        if(timedOut||restartHit) return false;
    }
    return false;
}

static void resetState(){ L.assign(N,vector<int>(N,-1)); rowUsed.assign(N,0);colUsed.assign(N,0);
    dRowUsed.assign(N,0);dColUsed.assign(N,0); Dval.assign(N,vector<int>(N,-1)); nodes=0; }

static void parseAut(const string&mode){
    g.assign(N,-1); haveAut=true;
    auto setCycle=[&](vector<int>cyc){int k=cyc.size();for(int i=0;i<k;++i)g[cyc[i]]=cyc[(i+1)%k];};
    if(mode=="raw"){haveAut=false;for(int i=0;i<N;++i)g[i]=i;autOrder=1;return;}
    else if(mode=="invol"){for(int i=0;i+1<N;i+=2)setCycle({i,i+1});if(N%2==1)g[N-1]=N-1;}
    else if(mode=="ncyc"){vector<int>c;for(int i=0;i<N;++i)c.push_back(i);setCycle(c);}
    else if(mode.rfind("cyc:",0)==0){ for(int i=0;i<N;++i)g[i]=i; string spec=mode.substr(4); size_t pos=0;
        while(pos<spec.size()){ size_t bar=spec.find('|',pos); string cyc=(bar==string::npos)?spec.substr(pos):spec.substr(pos,bar-pos);
            vector<int>elems;size_t p=0; while(p<cyc.size()){size_t cm=cyc.find(',',p);string num=(cm==string::npos)?cyc.substr(p):cyc.substr(p,cm-p);
                if(!num.empty())elems.push_back(stoi(num)); if(cm==string::npos)break;p=cm+1;}
            if(!elems.empty())setCycle(elems); if(bar==string::npos)break;pos=bar+1; } }
    else { fprintf(stderr,"unknown mode %s\n",mode.c_str()); exit(2); }
    vector<int>seen(N,0); for(int i=0;i<N;i++){if(g[i]<0||g[i]>=N||seen[g[i]]){fprintf(stderr,"bad perm\n");exit(2);}seen[g[i]]=1;}
    vector<int>vis(N,0);long long ord=1;auto lcm=[](long long a,long long b){return a/std::__gcd(a,b)*b;};
    for(int i=0;i<N;i++){if(vis[i])continue;int len=0,x=i;while(!vis[x]){vis[x]=1;x=g[x];len++;}ord=lcm(ord,len);}
    autOrder=(int)ord;
}

int main(int argc,char**argv){
    if(argc<5){fprintf(stderr,"usage: solve18 <n> <mode> <time_ms> <seed> [idem] [restart_base]\n");return 1;}
    N=atoi(argv[1]); string mode=argv[2]; long long tbox=atoll(argv[3]); unsigned seed=(unsigned)atoll(argv[4]);
    bool idem=false; long long restartBase=2000000;
    for(int i=5;i<argc;++i){string a=argv[i]; if(a=="idem")idem=true; else restartBase=atoll(a.c_str());}
    FULL=(N==32)?-1:((1<<N)-1);
    parseAut(mode);
    startTime=chrono::steady_clock::now(); deadline_ms=tbox;
    rng.seed(seed);
    resetState(); buildReps();

    int restart=0; restartLimit=restartBase;
    bool found=false;
    while(!found && !timedOut){
        restartHit=false;
        resetState();
        vector<CellSet> log;
        bool preOk=true;
        if(idem){ for(int i=0;i<N&&preOk;++i){ if(L[i][i]!=-1){if(L[i][i]!=i)preOk=false;continue;} if(!assignOrbit(i,i,i,log))preOk=false; } }
        if(preOk) found=backtrack(log);
        if(found&&solved) break;
        if(timedOut) break;
        // restart: grow the limit geometrically, advance seed
        restart++; restartLimit = restartBase * (1 + restart);  // linear-ish growth
        rng.seed(seed + 0x9E3779B9u*restart);
        if(restart%5==0) fprintf(stderr,"  [n=%d %s seed=%u] restart %d, limit=%lld, %lldms\n",N,mode.c_str(),seed,restart,restartLimit,now_ms());
    }

    long long elapsed=now_ms();
    string status = (found&&solved)?"FOUND":(timedOut?"TIMEOUT":"EXHAUSTED");
    fprintf(stderr,"n=%d mode=%s seed=%u idem=%d order=%d reps=%zu restarts=%d nodes=%lld elapsed_ms=%lld -> %s\n",
            N,mode.c_str(),seed,(int)idem,autOrder,repCells.size(),restart,nodes,elapsed,status.c_str());
    if(found&&solved){
        for(int r=0;r<N;++r){string s;for(int c=0;c<N;++c){s+=to_string(L[r][c]);if(c+1<N)s+=" ";}printf("%s\n",s.c_str());}
        for(string fn : {string("cpp_sol_n")+to_string(N)+"_s"+to_string(seed)+".json", string("cpp_sol_n")+to_string(N)+".json"}){
            FILE*f=fopen(fn.c_str(),"w"); fprintf(f,"{\"n\":%d,\"seed\":%u,\"mode\":\"%s\",\"L\":[",N,seed,mode.c_str());
            for(int r=0;r<N;++r){fprintf(f,"[");for(int c=0;c<N;++c)fprintf(f,"%d%s",L[r][c],c+1<N?",":"");fprintf(f,"]%s",r+1<N?",":"");}
            fprintf(f,"]}\n"); fclose(f);
        }
        fprintf(stderr,"wrote cpp_sol_n%d_s%u.json\n",N,seed);
    }
    return (found&&solved)?0:(timedOut?3:4);
}
