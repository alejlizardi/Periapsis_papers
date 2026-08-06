"""T6 independent certification of census witnesses / near-min survivors.
For each perm: certify rho with (a) cpsat_rho.py (CP-SAT arc-MTZ MILP — different
formalism), and (b) pc.py (state-dedup DFS reference). Report agreement vs the
census's claimed rho. Zero tolerance.

Usage: py -3 cert_witnesses.py <claimed_rho> <perm1> <perm2> ...
   or:  py -3 cert_witnesses.py --file <perms_file> [claimed_rho_col]
perm = comma-separated 0-based. Prints PASS/FAIL per perm and a summary.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pc
import cpsat_rho

def certify(perm, claimed):
    s = tuple(perm)
    t0=time.time(); r_pc = pc.rho(s); t_pc=time.time()-t0
    t0=time.time(); st, r_cp = cpsat_rho.rho_cpsat(list(perm), timeout=600.0); t_cp=time.time()-t0
    ok = (r_pc == claimed) and (r_cp == claimed) and (st == "OPTIMAL")
    print(f"  perm={','.join(map(str,perm))} claimed={claimed} pc.py={r_pc}({t_pc:.1f}s) "
          f"cpsat={r_cp}[{st}]({t_cp:.1f}s)  {'PASS' if ok else 'FAIL'}")
    return ok

def main():
    args=sys.argv[1:]
    perms=[]  # (perm, claimed)
    if args and args[0]=="--file":
        path=args[1]
        for line in open(path):
            line=line.strip()
            if not line or line.startswith("#"): continue
            # accept "rho perm" or "SURV rho b orbit [..]" or just "perm"
            parts=line.replace("["," ").replace("]"," ").split()
            if parts[0]=="SURV":
                claimed=int(parts[1]); perm=[int(x) for x in parts[4].split(",")]
            else:
                claimed=int(parts[0]); perm=[int(x) for x in parts[1].split(",")]
            perms.append((perm,claimed))
    else:
        claimed=int(args[0])
        for a in args[1:]:
            perms.append(([int(x) for x in a.split(",")], claimed))
    allok=True
    for perm,claimed in perms:
        allok &= certify(perm,claimed)
    print("CERT SUMMARY:", "ALL PASS" if allok else "SOME FAIL", f"({len(perms)} perms)")

if __name__=="__main__": main()
