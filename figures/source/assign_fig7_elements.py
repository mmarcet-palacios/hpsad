import csv, re, json, itertools
D="/mnt/user-data/uploads/Promoter/HPSAD_v7/datasets/motifs"
RUNS=[("C1 Pol II","C1_pol2_d"),("C1 coding","C1_protein_coding_d"),("C1 lncRNA","C1_non_coding_d"),("C1 small RNA","C1_small_rna_d"),
      ("C1 pseudo","C1_pseudogene_d"),("C2","C2_d"),("C2b","C2b_d"),("C3","C3_d"),("C4 Pol II","C4_pol2_d"),("C4 small RNA","C4_small_rna_d")]
IU={"A":"A","C":"C","G":"G","T":"T","R":"[AGR]","Y":"[CTY]","S":"[CGS]","W":"[ATW]","K":"[GTK]","M":"[ACM]","N":"."}
def rc(s): return s[::-1].translate(str.maketrans("ACGTRYKMSWN","TGCAYRMKSWN"))
VOC=[("TATA box",["TATAAA","WWWAAA","ATAAA"],(-40,-20)),
     ("Initiator",["CAG","TCA","CAT","CAST"],(-8,2)),
     ("YY1",["AAGATGG","ATGGCGGC","GCCATNTT"],(-20,25)),
     ("CCAAT box",["CCAAT","ATTGG"],(-200,-30)),
     ("GC box",["GGGCGG","CCGCCC","GGGGCG"],(-200,-1)),
     ("ETS",["CCGGAA","TTCCGG","GGAAGT","ACTTCC"],(-200,-1)),
     ("CRE",["TGACG","CGTCA"],(-200,-1)),
     ("lncRNA CRE-like",["TGTGA.NTCA","TGA.RTCA","CACAATG","CATTGTGA"],(-200,-1)),
     ("NRF1",["GCGCA.GCGC","TGCGCA"],(-200,-1)),
     ("E-box",["CACGTG"],(-200,-1)),
     ("THAP11/Ronin",["ACTACA","TGTAGT"],(-200,-1)),
     ("ZBTB33",["TCTCGCGAGA","CTCGCGAG"],(-200,-1)),
     ("CCTTTAAGG (ZNF536-like)",["CCTTTAAG","CTTAAAGG"],(-200,-1)),
     ("Start codon context",["CCATGG","ACCATG","CATGG"],(5,100)),
     ("5' splice donor",["GGTRAG","GTAAGT","GGTAAG","AGGTA","GTGAG"],(20,100)),
     ("tRNA A box",["TAGC.CAG","AG..GG","AGT.G"],(1,25)),
     ("tRNA B box",["GGTTCGA","GTTC.A..C"],(35,70)),
     ("C/D box C (snoRNA)",["TGATGA"],(1,40))]
SETS={"A":"A","C":"C","G":"G","T":"T","R":"AG","Y":"CT","S":"CG","W":"AT","K":"GT","M":"AC","N":"ACGT"}
def match(cons,pats):
    for p in pats:
        for c in (cons,rc(cons)):
            for i in range(len(c)-len(p)+1):
                ok=True
                for a,b in zip(c[i:i+len(p)],p):
                    if b==".": continue
                    sa=SETS[a]; sb=SETS[b]
                    if len(sa)>2 or not (set(sa)&set(sb)): ok=False;break
                if ok: return True
    return False
out={}
copyl=[]
for lab,run in RUNS:
    rows=[r for r in csv.DictReader(open(f"{D}/hpsad_motifs_{run}_zones_f.tsv"),delimiter="\t") if r["pass_fdr"]=="True"]
    out[lab]={"n_pass":len(rows),"elements":{},"unassigned":[]}
    for r in rows:
        pos=int(r["position"]); n=int(r["instances"]); ratio=float(r["ratio"]); cons=r["consensus"]
        hit=None
        if r["zone"]=="core":
            for nm,pats,(lo,hi) in VOC:
                if lo<=pos<=hi and match(cons,pats):
                    if nm=="Initiator" and len(cons.strip("N"))>11: continue
                    if nm in("tRNA A box","tRNA B box","C/D box C (snoRNA)") and "small" not in lab: continue
                    if nm in("Start codon context","5' splice donor") and "small" in lab: continue
                    hit=nm;break
        if hit:
            e=out[lab]["elements"].setdefault(hit,[])
            e.append(dict(pos=pos,n=n,ratio=ratio,cons=cons,score=float(r["score"])))
        else:
            out[lab]["unassigned"].append(dict(zone=r["zone"],pos=pos,n=n,ratio=ratio,cons=cons))
            copyl.append((lab,r["zone"],pos,n,ratio,cons))
json.dump(out,open("fig7_f.json","w"),indent=1)
names=[v[0] for v in VOC]
print(f"{'':26}"+"".join(f"{l[:9]:>10}" for l,_ in RUNS))
for nm in names:
    line=f"{nm[:25]:26}"
    for l,_ in RUNS:
        e=out[l]["elements"].get(nm)
        if e:
            b=max(e,key=lambda x:x['score']); line+=f"{b['pos']:>+6}({len(e)})"
        else: line+="         ."
    print(line)
print(f"{'passing':26}"+"".join(f"{out[l]['n_pass']:>10}" for l,_ in RUNS))
print(f"{'unassigned':26}"+"".join(f"{len(out[l]['unassigned']):>10}" for l,_ in RUNS))
print(f"{'  of which O/E>20,n<150':26}"+"".join(f"{sum(1 for u in out[l]['unassigned'] if u['ratio']>20 and u['n']<150):>10}" for l,_ in RUNS))
for l,_ in RUNS:
    print("\n##",l, "unassigned:", "; ".join(f"{u['cons']} {u['pos']:+d} x{u['ratio']:.0f} n{u['n']}" for u in out[l]["unassigned"]))
