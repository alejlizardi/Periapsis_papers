"""GPU census for rho_min(n) via CuPy RawKernel (NVRTC; no Visual Studio).

Thread-per-sigma: device-side Lehmer unrank -> D4 canonical check -> adjacency
bitmasks -> iterative-DFS threshold decision (early exit at t vertices).

Kernels:
  decide_kernel(start,count,n,t)      : out[i]=decide(rank=start+i, t)   [RAW, for validation]
  decide_perms_kernel(perms,count,n,t): out[i]=decide(perms[i], t)       [RAW, arbitrary perms]
  census_kernel(start,count,n,t)      : out[i]=1 iff sigma canonical AND no PC path >= t
                                        (rho<=t-1)                       [production survivor flag]

Modes:
  selftest
  validate_s9                     : decide_kernel vs CPU rho over ALL S_9, all t (bit-identical)
  validate_perms <n> <samples>    : decide_perms_kernel vs CPU rho on random perms, all t
  equiv13                         : full GPU census n=13 K=9 -> survivor set must equal CPU's 62
  census <n> <K> [batch] [outdir] : full checkpointed census (canonical sigma with rho<=K)
  bench <n> <t> <nranks>          : throughput of census_kernel over first nranks ranks
"""
import sys, os, time, subprocess, json
import cupy as cp
import numpy as np

KERNEL_SRC = r'''
extern "C" {
__device__ void unrank(unsigned long long rank,int n,int* s){
    int avail[16]; for(int i=0;i<n;i++) avail[i]=i;
    unsigned long long f=1; for(int i=2;i<n;i++) f*=(unsigned long long)i;   // (n-1)!
    int m=n;
    for(int i=0;i<n;i++){
        int q=(int)(rank/f); rank%=f;
        s[i]=avail[q];
        for(int j=q;j<m-1;j++) avail[j]=avail[j+1];
        m--; if(i<n-1) f/=(unsigned long long)(n-1-i);
    }
}
__device__ bool lex_less(const int* x,const int* s,int n){
    for(int i=0;i<n;i++){ if(x[i]<s[i])return true; if(x[i]>s[i])return false; } return false;
}
__device__ bool is_canon(const int* s,int n){
    int a[4][16];
    for(int i=0;i<n;i++){ a[0][i]=s[i]; a[1][i]=s[n-1-i]; a[2][i]=n-1-s[i]; a[3][i]=n-1-s[n-1-i]; }
    for(int k=0;k<4;k++){
        if(k>0 && lex_less(a[k],s,n)) return false;
        int inv[16]; for(int i=0;i<n;i++) inv[a[k][i]]=i;
        if(lex_less(inv,s,n)) return false;
    }
    return true;
}
__device__ bool decide(const unsigned int* adj0,const unsigned int* adj1,int n,int t){
    if(t<=1) return true;
    int stC[16]; unsigned int stVis[16], stCand[16];
    for(int sv=0; sv<n; sv++) for(int sc=0; sc<2; sc++){
        int sp=0; stC[0]=sc; stVis[0]=(1u<<sv);
        stCand[0]=(sc==0?adj0[sv]:adj1[sv]) & ~stVis[0];
        while(sp>=0){
            unsigned int cand=stCand[sp];
            if(!cand){ sp--; continue; }
            int w=__ffs((int)cand)-1;
            stCand[sp]=cand & ~(1u<<w);
            int c=stC[sp];
            unsigned int nvis=stVis[sp] | (1u<<w);
            int nlen=sp+2;
            if(nlen>=t) return true;
            int nc=1-c; sp++;
            stC[sp]=nc; stVis[sp]=nvis;
            stCand[sp]=(nc==0?adj0[w]:adj1[w]) & ~nvis;
        }
    }
    return false;
}
__device__ void build_adj(const int* s,int n,unsigned int* adj0,unsigned int* adj1){
    int pos[16]; for(int i=0;i<n;i++) pos[s[i]]=i;
    for(int i=0;i<n;i++){adj0[i]=0;adj1[i]=0;}
    for(int i=0;i+1<n;i++){adj0[i]|=(1u<<(i+1));adj0[i+1]|=(1u<<i);}
    for(int v=0;v+1<n;v++){int a=pos[v],b=pos[v+1];adj1[a]|=(1u<<b);adj1[b]|=(1u<<a);}
}
__global__ void decide_kernel(unsigned long long start,unsigned long long count,int n,int t,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; unrank(start+i,n,s);
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    out[i]=decide(a0,a1,n,t)?1:0;
}
__global__ void decide_perms_kernel(const int* perms,unsigned long long count,int n,int t,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; for(int k=0;k<n;k++) s[k]=perms[i*n+k];
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    out[i]=decide(a0,a1,n,t)?1:0;
}
__global__ void census_kernel(unsigned long long start,unsigned long long count,int n,int t,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; unrank(start+i,n,s);
    if(!is_canon(s,n)){ out[i]=0; return; }
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    out[i]=decide(a0,a1,n,t)?0:1;
}
}
'''
_mod = cp.RawModule(code=KERNEL_SRC, options=('--std=c++11',))
K_decide = _mod.get_function("decide_kernel")
K_decide_perms = _mod.get_function("decide_perms_kernel")
K_census = _mod.get_function("census_kernel")

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS_EXE = os.path.join(HERE, "census.exe")
THREADS = 256

def factorial(n):
    f=1
    for i in range(2,n+1): f*=i
    return f

def unrank_host(rank,n):
    avail=list(range(n)); s=[0]*n; f=factorial(n-1)
    for i in range(n):
        q=rank//f; rank%=f; s[i]=avail.pop(q)
        if i<n-1: f//=(n-1-i)
    return s

def _launch_rank(kern,start,count,n,t):
    out=cp.empty(int(count),dtype=cp.uint8)
    blocks=(int(count)+THREADS-1)//THREADS
    kern((blocks,),(THREADS,),(np.uint64(start),np.uint64(count),np.int32(n),np.int32(t),out))
    return out

def cpu_rho_batch(perms,n):
    """Exact rho for a list of perms via census.exe rhofile (one process)."""
    inp="\n".join(",".join(map(str,p)) for p in perms)+"\n"
    out=subprocess.run([CENSUS_EXE,"rhofile",str(n)],input=inp,capture_output=True,text=True).stdout
    return np.array([int(x) for x in out.split()],dtype=np.int32)

def validate_s9():
    n=9; total=factorial(n)
    dump=subprocess.run([CENSUS_EXE,"dumprho",str(n)],capture_output=True,text=True).stdout.splitlines()
    assert len(dump)==total,(len(dump),total)
    cpu=np.array([int(l.rsplit(" ",1)[1]) for l in dump],dtype=np.int32)
    mism=0
    for t in range(1,n+1):
        gpu=cp.asnumpy(_launch_rank(K_decide,0,total,n,t)).astype(bool)
        exp=(cpu>=t); d=int((gpu!=exp).sum()); mism+=d
        if d:
            for idx in np.where(gpu!=exp)[0][:5]:
                print(f"  MISMATCH t={t} rank={idx} perm={unrank_host(int(idx),n)} gpu={int(gpu[idx])} cpu_rho={cpu[idx]}")
    print(f"[gpu validate_s9] all t, mismatches={mism}  {'PASS' if mism==0 else 'FAIL'}")
    return mism==0

def validate_perms(n,samples,seed=777):
    rng=np.random.default_rng(seed)
    perms=np.empty((samples,n),dtype=np.int32)
    base=np.arange(n)
    for i in range(samples):
        perms[i]=rng.permutation(base)
    cpu=cpu_rho_batch(perms.tolist(),n)
    d_perms=cp.asarray(perms.reshape(-1))
    mism=0
    for t in range(1,n+1):
        out=cp.empty(samples,dtype=cp.uint8)
        blocks=(samples+THREADS-1)//THREADS
        K_decide_perms((blocks,),(THREADS,),(d_perms,np.uint64(samples),np.int32(n),np.int32(t),out))
        gpu=cp.asnumpy(out).astype(bool); exp=(cpu>=t); dd=int((gpu!=exp).sum()); mism+=dd
        if dd:
            for idx in np.where(gpu!=exp)[0][:5]:
                print(f"  MISMATCH n={n} t={t} perm={perms[idx].tolist()} gpu={int(gpu[idx])} cpu_rho={cpu[idx]}")
    print(f"[gpu validate_perms] n={n} samples={samples} all t, mismatches={mism}  {'PASS' if mism==0 else 'FAIL'}")
    return mism==0

def gpu_census(n,K,batch=1<<26,outdir=None,report=True):
    """census_kernel over all n! ranks at t=K+1; return list of survivor perms (canonical, rho<=K).
    Checkpointed: writes cursor + survivor ranks to outdir if given."""
    total=factorial(n); t=K+1
    outdir=outdir or f"gpu_census_n{n}_K{K}"
    os.makedirs(outdir,exist_ok=True)
    cur_path=os.path.join(outdir,"cursor.txt"); surv_path=os.path.join(outdir,"survivor_ranks.txt")
    meta_path=os.path.join(outdir,"meta.json")
    def _write_atomic(path,text):
        tmp=path+".tmp"
        with open(tmp,"w") as f: f.write(text); f.flush(); os.fsync(f.fileno())
        for _ in range(20):            # Windows os.replace can transiently WinError 5 (AV/indexer lock)
            try:
                os.replace(tmp,path); return
            except PermissionError:
                time.sleep(0.25)
        # persistent lock: fall back to a direct write (non-atomic, but keeps the run advancing)
        with open(path,"w") as f: f.write(text); f.flush(); os.fsync(f.fileno())
        try: os.remove(tmp)
        except OSError: pass
    def _read_cursor():
        # tolerate a torn .tmp; the committed file (if present) is always intact
        for p in (cur_path,):
            if os.path.exists(p):
                try: return int(open(p).read().strip())
                except Exception: return 0
        return 0
    # param-guard: refuse to resume a checkpoint written for different (n,K,total)
    if os.path.exists(meta_path):
        try: meta=json.load(open(meta_path))
        except Exception: meta={}
        if meta and (meta.get("n")!=n or meta.get("K")!=K or meta.get("total")!=total):
            raise SystemExit(f"[abort] checkpoint {outdir} is for n={meta.get('n')} K={meta.get('K')} "
                             f"(total={meta.get('total')}), not n={n} K={K}. Use a fresh outdir.")
    _write_atomic(meta_path, json.dumps({"n":n,"K":K,"total":total,"batch":batch}))
    start=_read_cursor()
    if start and report: print(f"[resume] from rank {start}/{total} ({100*start/total:.2f}%)",flush=True)
    surv_ranks=[]
    if os.path.exists(surv_path):
        for tok in open(surv_path).read().split():
            try:
                r=int(tok)
            except ValueError:
                continue                 # skip a torn last token from an aborted append
            if r<start: surv_ranks.append(r)   # keep only survivors from fully-committed batches
    t0=time.time(); done_at_start=start
    while start<total:
        count=min(batch,total-start)
        out=_launch_rank(K_census,start,count,n,t)
        idx=cp.nonzero(out)[0]
        if idx.size:
            rs=(cp.asnumpy(idx).astype(np.uint64)+np.uint64(start)).tolist()
            surv_ranks.extend(int(r) for r in rs)
            with open(surv_path,"a") as f:   # append survivors, then fsync BEFORE advancing cursor
                for r in rs: f.write(f"{r}\n")
                f.flush(); os.fsync(f.fileno())
        start+=count
        _write_atomic(cur_path,str(start))   # atomic: cursor only advances once survivors are durable
        if report:
            el=time.time()-t0; rate=(start-done_at_start)/el if el>0 else 0
            eta=(total-start)/rate if rate>0 else 0
            print(f"  {start}/{total} ({100*start/total:.1f}%) surv={len(surv_ranks)} rate={rate/1e6:.1f}M/s eta={eta/60:.1f}min",flush=True)
    # dedup survivor ranks (in case of resume overlap)
    surv_ranks=sorted(set(surv_ranks))
    perms=[unrank_host(r,n) for r in surv_ranks]
    return perms

def resolve_and_report(n,K,perms,label=""):
    """CPU-resolve exact rho + b + orbit for each survivor perm; print rho_min + witnesses."""
    if not perms:
        print(f"=== GPU n={n} K={K}: NO survivors -> rho_min(n) > {K}"); return None,[]
    # batch resolve via census.exe singlefile: "rho b orbit canon" per line
    inp="\n".join(",".join(map(str,p)) for p in perms)+"\n"
    out=subprocess.run([CENSUS_EXE,"singlefile",str(n)],input=inp,capture_output=True,text=True).stdout.splitlines()
    assert len(out)==len(perms),(len(out),len(perms))
    raw=[]
    for p,o in zip(perms,out):
        rho,b,orb,can=(int(x) for x in o.split())
        raw.append((rho,b,orb,can,p))
    # SELF-HEAL: a valid census survivor is canonical AND has rho<=K. Drop anything else
    # (would only appear from checkpoint corruption); the CPU resolver is ground truth here.
    lines=[l for l in raw if l[3]==1 and l[0]<=K]
    dropped=[l for l in raw if not (l[3]==1 and l[0]<=K)]
    if dropped:
        print(f"  [self-heal] dropped {len(dropped)} survivor(s) failing canon/rho<=K "
              f"(checkpoint artifact); e.g. {dropped[0][:4]}")
    if not lines:
        print(f"=== GPU {label} n={n} K={K}: NO valid survivors -> rho_min(n) > {K}"); return None,[]
    mn=min(l[0] for l in lines)
    wit=[l for l in lines if l[0]==mn]
    full=sum(l[2] for l in wit)
    print(f"=== GPU {label} n={n} K={K} RESULT ===")
    print(f"rho_min({n}) = {mn}")
    print(f"survivors(rho<=K): {len(lines)}   argmin(rho={mn}) canonical: {len(wit)}  full-orbit: {full}")
    fromb={}
    for rho,b,orb,can,p in wit: fromb[b]=fromb.get(b,0)+orb
    print("witness full-count by b:",dict(sorted(fromb.items())))
    return mn,wit

if __name__=="__main__":
    m=sys.argv[1] if len(sys.argv)>1 else "selftest"
    if m=="selftest":
        print("device:",cp.cuda.runtime.getDeviceProperties(0)['name'].decode())
        print("unrank0 n=13:",unrank_host(0,13))
    elif m=="validate_s9":
        validate_s9()
    elif m=="validate_perms":
        validate_perms(int(sys.argv[2]),int(sys.argv[3]))
    elif m=="equiv13":
        t0=time.time()
        perms=gpu_census(13,9,outdir="gpu_census_n13_K9")
        print(f"gpu census n=13 wall={time.time()-t0:.1f}s survivors={len(perms)}")
        resolve_and_report(13,9,perms,label="equiv13")
    elif m=="census":
        n=int(sys.argv[2]); K=int(sys.argv[3])
        batch=int(sys.argv[4]) if len(sys.argv)>4 else (1<<26)
        outdir=sys.argv[5] if len(sys.argv)>5 else f"gpu_census_n{n}_K{K}"
        t0=time.time()
        perms=gpu_census(n,K,batch=batch,outdir=outdir)
        print(f"gpu census n={n} wall={time.time()-t0:.1f}s survivors={len(perms)}")
        resolve_and_report(n,K,perms,label="census")
    elif m=="bench":
        n=int(sys.argv[2]); t=int(sys.argv[3]); nr=int(sys.argv[4])
        t0=time.time(); out=_launch_rank(K_census,0,nr,n,t); cp.cuda.Stream.null.synchronize()
        dt=time.time()-t0
        print(f"[bench] n={n} t={t} nranks={nr} wall={dt:.3f}s rate={nr/dt/1e6:.1f}M sigma/s survivors={int(cp.count_nonzero(out))}")
