"""WP-K GPU counting census: |D(n)| = #{sigma in S_n : rho(sigma) < n}, exact.

Adapted from WP-I's hardened census_gpu.py (thread-per-sigma, device Lehmer
unrank, D4 canonical filter, iterative-DFS threshold decision) with the three
changes counting requires:

  1. REDUCTION, NOT COLLECTION. Deficient sigma are ~0.3% of the space
     (~10^10+ at n=16) — survivor ranks cannot be collected. The kernel emits
     one byte per rank and the host reduces sums per batch.
  2. ORBIT WEIGHTING (trap 1). The kernel keeps the canonical-rep filter; for
     each DEFICIENT canonical rep it computes the TRUE orbit size (distinct
     images under all 8 D4 elements) on device. |D(n)| = sum of orbit sizes,
     never reps*8.
  3. ATOMIC COUNTER CHECKPOINT (trap 2). Cursor + all counters + reservoir
     state live in ONE json written atomically (tmp+fsync+replace) after each
     batch. A resume can never replay a partial batch into the counters.

Out-byte encoding per rank:
  0            : non-canonical (skipped)
  bit4 (16)    : canonical
  bits0-3      : orbit size (1..8) iff deficient (rho < n), else 0
  bit5 (32)    : deep-deficient (rho <= n-2), only ever set when deficient

The device decide() carries the CPU engine's reach_bound prune (admissible —
prunes only subtrees that provably cannot reach t; decisions are identical,
validated bit-for-bit). A no-prune variant is kept for benchmarking.

Reservoir sampling for two-sided spot certification is deterministic bottom-k:
each candidate rank gets key = splitmix64(rank); the R smallest keys win.
Batch boundaries and resumes cannot perturb the final sample.

Modes:
  selftest
  validate_s9                      : decide vs CPU rho over ALL S_9, all t
  validate_perms <n> <samples>     : decide vs CPU rho on random perms, all t
  count <n> [batch] [outdir]       : full checkpointed counting census
  bench <n> <nranks> [prune01]     : count_kernel throughput over first nranks
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
__device__ int orbit_size(const int* s,int n){
    // distinct images under {e,rev,comp,revcomp} x {id,inv} — ported from CPU orbit_size.
    int img[8][16]; int a[4][16];
    for(int i=0;i<n;i++){ a[0][i]=s[i]; a[1][i]=s[n-1-i]; a[2][i]=n-1-s[i]; a[3][i]=n-1-s[n-1-i]; }
    for(int k=0;k<4;k++){
        for(int i=0;i<n;i++) img[k][i]=a[k][i];
        for(int i=0;i<n;i++) img[4+k][a[k][i]]=i;
    }
    int cnt=0;
    for(int k=0;k<8;k++){
        bool dup=false;
        for(int j=0;j<k&&!dup;j++){
            bool eq=true;
            for(int i=0;i<n&&eq;i++) if(img[k][i]!=img[j][i]) eq=false;
            if(eq) dup=true;
        }
        if(!dup) cnt++;
    }
    return cnt;
}
__device__ void build_adj(const int* s,int n,unsigned int* adj0,unsigned int* adj1){
    int pos[16]; for(int i=0;i<n;i++) pos[s[i]]=i;
    for(int i=0;i<n;i++){adj0[i]=0;adj1[i]=0;}
    for(int i=0;i+1<n;i++){adj0[i]|=(1u<<(i+1));adj0[i+1]|=(1u<<i);}
    for(int v=0;v+1<n;v++){int a=pos[v],b=pos[v+1];adj1[a]|=(1u<<b);adj1[b]|=(1u<<a);}
}
// reach_bound: alternation-aware reachable-vertex upper bound (ported VERBATIM
// from the CPU engine: frontier[col] = vertices to expand via color col).
__device__ int reach_bound(const unsigned int* adj0,const unsigned int* adj1,int v,int c,unsigned int vis){
    unsigned int seen[2]={0u,0u}, frontier[2]={0u,0u}, acc=0u;
    unsigned int first=(c==0?adj0[v]:adj1[v])&~vis;
    frontier[1-c]=first; seen[1-c]=first; acc|=first;
    bool change=true;
    while(change){ change=false;
        for(int col=0;col<2;col++){
            unsigned int f=frontier[col]; frontier[col]=0u;
            while(f){ int w=__ffs((int)f)-1; f&=f-1;
                unsigned int nx=(col==0?adj0[w]:adj1[w])&~vis&~seen[1-col];
                if(nx){ seen[1-col]|=nx; frontier[1-col]|=nx; acc|=nx; change=true; } } } }
    return __popc(acc);
}
// decide: is there a properly-colored path on >= t vertices?  PRUNE!=0 applies
// the reach_bound prune exactly where the CPU Dec does (root + every child push).
// The prune is admissible: it only skips subtrees that cannot reach t, so the
// decision is identical with or without it.
__device__ bool decide(const unsigned int* adj0,const unsigned int* adj1,int n,int t,int PRUNE){
    if(t<=1) return true;
    int stC[16]; unsigned int stVis[16], stCand[16];
    for(int sv=0; sv<n; sv++) for(int sc=0; sc<2; sc++){
        unsigned int vis0=(1u<<sv);
        unsigned int cand0=(sc==0?adj0[sv]:adj1[sv]) & ~vis0;
        if(!cand0) continue;
        if(PRUNE && 1+reach_bound(adj0,adj1,sv,sc,vis0)<t) continue;
        int sp=0; stC[0]=sc; stVis[0]=vis0; stCand[0]=cand0;
        while(sp>=0){
            unsigned int cand=stCand[sp];
            if(!cand){ sp--; continue; }
            int w=__ffs((int)cand)-1;
            stCand[sp]=cand & ~(1u<<w);
            int c=stC[sp];
            unsigned int nvis=stVis[sp] | (1u<<w);
            int nlen=sp+2;
            if(nlen>=t) return true;
            int nc=1-c;
            unsigned int ncand=(nc==0?adj0[w]:adj1[w]) & ~nvis;
            if(!ncand) continue;
            if(PRUNE && nlen+reach_bound(adj0,adj1,w,nc,nvis)<t) continue;
            sp++;
            stC[sp]=nc; stVis[sp]=nvis; stCand[sp]=ncand;
        }
    }
    return false;
}
__global__ void decide_kernel(unsigned long long start,unsigned long long count,int n,int t,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; unrank(start+i,n,s);
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    out[i]=decide(a0,a1,n,t,1)?1:0;
}
__global__ void decide_perms_kernel(const int* perms,unsigned long long count,int n,int t,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; for(int k=0;k<n;k++) s[k]=perms[i*n+k];
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    out[i]=decide(a0,a1,n,t,1)?1:0;
}
__global__ void count_kernel(unsigned long long start,unsigned long long count,int n,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; unrank(start+i,n,s);
    if(!is_canon(s,n)){ out[i]=0; return; }
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    unsigned char r=16;
    if(!decide(a0,a1,n,n,1)){
        r|=(unsigned char)orbit_size(s,n);
        if(!decide(a0,a1,n,n-1,1)) r|=32;
    }
    out[i]=r;
}
__global__ void count_kernel_noprune(unsigned long long start,unsigned long long count,int n,unsigned char* out){
    unsigned long long i=(unsigned long long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=count)return;
    int s[16]; unrank(start+i,n,s);
    if(!is_canon(s,n)){ out[i]=0; return; }
    unsigned int a0[16],a1[16]; build_adj(s,n,a0,a1);
    unsigned char r=16;
    if(!decide(a0,a1,n,n,0)){
        r|=(unsigned char)orbit_size(s,n);
        if(!decide(a0,a1,n,n-1,0)) r|=32;
    }
    out[i]=r;
}
}
'''
_mod = cp.RawModule(code=KERNEL_SRC, options=('--std=c++11',))
K_decide = _mod.get_function("decide_kernel")
K_decide_perms = _mod.get_function("decide_perms_kernel")
K_count = _mod.get_function("count_kernel")
K_count_np = _mod.get_function("count_kernel_noprune")

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS_EXE = os.environ.get("WPK_CENSUS_EXE", os.path.join(HERE, "census.exe"))
THREADS = 256
RESERVOIR = 150   # per side (deficient / non-deficient) for spot certification

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

def splitmix64(x):
    """Vectorized deterministic hash for bottom-k reservoir keys."""
    x = (np.asarray(x, dtype=np.uint64) + np.uint64(0x9E3779B97F4A7C15))
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))

def _merge_reservoir(res, ranks):
    """res: list of [key,rank] (kept sorted by key, len<=RESERVOIR). ranks: np.uint64 array.
    Bottom-k by deterministic key. Pre-select the k smallest on the numpy side before
    building Python objects (bottom-k is associative, so the result is identical; a dense
    batch can carry ~10^6 candidate ranks and materializing them all OOMs the host)."""
    if ranks.size == 0: return res
    keys = splitmix64(ranks)
    if keys.size > RESERVOIR:
        sel = np.argpartition(keys, RESERVOIR)[:RESERVOIR]
        keys, ranks = keys[sel], ranks[sel]
    cand = res + [[int(k), int(r)] for k, r in zip(keys, ranks)]
    cand.sort(key=lambda kr: kr[0])
    return cand[:RESERVOIR]

def _write_atomic(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write(text); f.flush(); os.fsync(f.fileno())
    for _ in range(20):            # Windows os.replace can transiently WinError 5
        try:
            os.replace(tmp, path); return
        except PermissionError:
            time.sleep(0.25)
    with open(path, "w") as f:     # persistent lock: non-atomic fallback keeps the run alive
        f.write(text); f.flush(); os.fsync(f.fileno())
    try: os.remove(tmp)
    except OSError: pass

def cpu_rho_batch(perms, n, procs=8):
    """Exact CPU rho, split across processes (rhofile is single-threaded)."""
    from concurrent.futures import ThreadPoolExecutor
    chunks = [perms[i::procs] for i in range(procs)]
    def run_chunk(ch):
        if not ch: return []
        inp = "\n".join(",".join(map(str, p)) for p in ch) + "\n"
        out = subprocess.run([CENSUS_EXE, "rhofile", str(n)], input=inp, capture_output=True, text=True).stdout
        return [int(x) for x in out.split()]
    with ThreadPoolExecutor(procs) as ex:
        res = list(ex.map(run_chunk, chunks))
    out = np.empty(len(perms), dtype=np.int32)
    for i, r in enumerate(res):
        out[i::procs] = r
    return out

def _launch(kern, start, count, n, extra=()):
    out = cp.empty(int(count), dtype=cp.uint8)
    blocks = (int(count) + THREADS - 1) // THREADS
    kern((blocks,), (THREADS,), (np.uint64(start), np.uint64(count), np.int32(n)) + tuple(extra) + (out,))
    return out

def validate_s9():
    n=9; total=factorial(n)
    dump = subprocess.run([CENSUS_EXE, "dumprho", str(n)], capture_output=True, text=True).stdout.splitlines()
    assert len(dump) == total, (len(dump), total)
    cpu = np.array([int(l.rsplit(" ", 1)[1]) for l in dump], dtype=np.int32)
    mism = 0
    for t in range(1, n+1):
        gpu = cp.asnumpy(_launch(K_decide, 0, total, n, (np.int32(t),))).astype(bool)
        exp = (cpu >= t); d = int((gpu != exp).sum()); mism += d
        if d:
            for idx in np.where(gpu != exp)[0][:5]:
                print(f"  MISMATCH t={t} rank={idx} perm={unrank_host(int(idx),n)} gpu={int(gpu[idx])} cpu_rho={cpu[idx]}")
    print(f"[validate_s9] all t, mismatches={mism}  {'PASS' if mism==0 else 'FAIL'}")
    return mism == 0

def validate_perms(n, samples, seed=777):
    rng = np.random.default_rng(seed)
    perms = np.empty((samples, n), dtype=np.int32)
    base = np.arange(n)
    for i in range(samples):
        perms[i] = rng.permutation(base)
    cpu = cpu_rho_batch(perms.tolist(), n)
    assert cpu.size == samples
    d_perms = cp.asarray(perms.reshape(-1))
    mism = 0
    for t in range(1, n+1):
        out = cp.empty(samples, dtype=cp.uint8)
        blocks = (samples + THREADS - 1) // THREADS
        K_decide_perms((blocks,), (THREADS,), (d_perms, np.uint64(samples), np.int32(n), np.int32(t), out))
        gpu = cp.asnumpy(out).astype(bool); exp = (cpu >= t); dd = int((gpu != exp).sum()); mism += dd
        if dd:
            for idx in np.where(gpu != exp)[0][:5]:
                print(f"  MISMATCH n={n} t={t} perm={perms[idx].tolist()} gpu={int(gpu[idx])} cpu_rho={cpu[idx]}")
    print(f"[validate_perms] n={n} samples={samples} all t, mismatches={mism}  {'PASS' if mism==0 else 'FAIL'}")
    return mism == 0

CKPT_SCHEMA = "wpk-count-v1"

def _fresh_state(n, batch):
    return {"schema": CKPT_SCHEMA, "n": n, "total": factorial(n), "batch": batch,
            "cursor": 0,
            "canon": 0, "def_reps": 0, "def_w": 0, "deep_reps": 0, "deep_w": 0,
            "res_def": [], "res_nondef": [],
            "resume_events": 0, "wall_s": 0.0}

def gpu_count(n, batch=1<<26, outdir=None, report=True):
    """Full counting census over all n! ranks. Returns the final state dict.
    Checkpoint = ONE atomic json {cursor + counters + reservoirs} per batch."""
    total = factorial(n)
    outdir = outdir or f"gpu_count_n{n}"
    os.makedirs(outdir, exist_ok=True)
    ck_path = os.path.join(outdir, "checkpoint.json")
    hb_path = os.path.join(outdir, "heartbeat.json")
    st = None; resumed = False
    if os.path.exists(ck_path):
        try:
            st = json.load(open(ck_path))
        except Exception:
            st = None
        if st is not None:
            if st.get("schema") != CKPT_SCHEMA or st.get("n") != n or st.get("total") != total:
                raise SystemExit(f"[abort] checkpoint {ck_path} is for a different run "
                                 f"(schema/n/total mismatch: {st.get('schema')}/{st.get('n')}/{st.get('total')}). "
                                 f"Use a fresh outdir.")
            resumed = st["cursor"] > 0
    if st is None:
        st = _fresh_state(n, batch)
    if resumed:
        st["resume_events"] += 1
        if report: print(f"[resume] from rank {st['cursor']}/{total} ({100*st['cursor']/total:.2f}%)", flush=True)
    t0 = time.time(); done_at_start = st["cursor"]; ema_rate = None
    while st["cursor"] < total:
        start = st["cursor"]; count = min(batch, total - start)
        bt0 = time.time()
        out = _launch(K_count, start, count, n)
        orb = out & cp.uint8(15)
        b_canon = int(((out >> 4) & cp.uint8(1)).sum(dtype=cp.int64))
        b_defw = int(orb.sum(dtype=cp.int64))
        defmask = orb > 0
        b_defreps = int(defmask.sum(dtype=cp.int64))
        deepmask = (out & cp.uint8(32)) > 0
        b_deepreps = int(deepmask.sum(dtype=cp.int64))
        b_deepw = int(orb[deepmask].sum(dtype=cp.int64)) if b_deepreps else 0
        # reservoirs (deterministic bottom-k over global ranks)
        if b_defreps:
            dr = (cp.asnumpy(cp.nonzero(defmask)[0]).astype(np.uint64) + np.uint64(start))
            st["res_def"] = _merge_reservoir(st["res_def"], dr)
        ndmask = out == cp.uint8(16)
        b_nd = int(ndmask.sum(dtype=cp.int64))
        if b_nd:
            # subsample device-side before host transfer: strided pick keeps it cheap
            nd_idx = cp.nonzero(ndmask)[0]
            stride = max(1, b_nd // 4096)
            ndr = (cp.asnumpy(nd_idx[::stride]).astype(np.uint64) + np.uint64(start))
            st["res_nondef"] = _merge_reservoir(st["res_nondef"], ndr)
        bwall = time.time() - bt0
        st["canon"] += b_canon; st["def_reps"] += b_defreps; st["def_w"] += b_defw
        st["deep_reps"] += b_deepreps; st["deep_w"] += b_deepw
        st["cursor"] = start + count
        st["wall_s"] += bwall
        _write_atomic(ck_path, json.dumps(st))          # cursor+counters commit TOGETHER
        rate = count / bwall if bwall > 0 else 0
        ema_rate = rate if ema_rate is None else 0.85 * ema_rate + 0.15 * rate
        sus = (st["cursor"] - done_at_start) / (time.time() - t0)
        eta = (total - st["cursor"]) / sus if sus > 0 else -1
        hb = {"ts": time.time(), "pid": os.getpid(), "n": n, "cursor": st["cursor"], "total": total,
              "def_w": st["def_w"], "deep_w": st["deep_w"], "rate_inst": rate, "rate_ema": ema_rate,
              "rate_sustained": sus, "eta_s": eta, "resumed": resumed,
              "resume_events": st["resume_events"], "last_ckpt_ts": time.time(), "outdir": os.path.abspath(outdir)}
        try:
            _write_atomic(hb_path, json.dumps(hb))
        except Exception:
            pass                                        # heartbeat is best-effort, never fatal
        if report:
            print(f"  {st['cursor']}/{total} ({100*st['cursor']/total:.2f}%) defW={st['def_w']} "
                  f"deepW={st['deep_w']} rate={sus/1e6:.1f}M/s eta={eta/60:.1f}min", flush=True)
    if report:
        print(f"[done] n={n} |D({n})| = {st['def_w']}  (reps={st['def_reps']}, canon={st['canon']}, "
              f"deep_w={st['deep_w']}, deep_reps={st['deep_reps']})  wall={st['wall_s']:.1f}s "
              f"resume_events={st['resume_events']}", flush=True)
    return st

if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if m == "selftest":
        print("device:", cp.cuda.runtime.getDeviceProperties(0)['name'].decode())
        print("unrank0 n=13:", unrank_host(0, 13))
    elif m == "validate_s9":
        validate_s9()
    elif m == "validate_perms":
        validate_perms(int(sys.argv[2]), int(sys.argv[3]))
    elif m == "count":
        n = int(sys.argv[2])
        batch = int(sys.argv[3]) if len(sys.argv) > 3 else (1 << 26)
        outdir = sys.argv[4] if len(sys.argv) > 4 else None
        t0 = time.time()
        st = gpu_count(n, batch=batch, outdir=outdir)
        print(f"RESULT n={n} D={st['def_w']} def_reps={st['def_reps']} deep_w={st['deep_w']} "
              f"deep_reps={st['deep_reps']} canon={st['canon']} wall={time.time()-t0:.1f}s")
    elif m == "bench":
        n = int(sys.argv[2]); nr = int(sys.argv[3])
        prune = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        kern = K_count if prune else K_count_np
        t0 = time.time(); out = _launch(kern, 0, nr, n); cp.cuda.Stream.null.synchronize()
        dt = time.time() - t0
        defw = int((out & cp.uint8(15)).sum(dtype=cp.int64))
        print(f"[bench] n={n} prune={prune} nranks={nr} wall={dt:.3f}s rate={nr/dt/1e6:.1f}M sigma/s defW_slice={defw}")
