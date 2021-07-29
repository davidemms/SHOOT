import sys
import glob

import ete3

def resolve_shoot_polytomies(d_db):
    """
    Resolve all the bifurcating trees in a shoot DB that EPA would have to use
    Args:
        d_db - The shoot database directory
    Post-condition:
        For each tree with a polytomy create a new file ORIG.bifur.tre for directories:
        - Gene_Trees/*_tree.txt
        - subtrees_*/sub/
        - subtrees_*/super/
    """
    if not d_db.endswith("/"):
        d_db += "/"
    targets = ["Gene_Trees/*_tree.txt", "Gene_Trees/subtrees*/sub/*.tre"]
    exclude = [".tre.unroot.tre", ".without.tre", ".super.tre.unroot.tre"]
    n_edit = 0
    try:
        for t in targets:
            for fn in glob.glob(d_db + t):
                if any(e in fn for e in exclude):
                    continue
                q_updated = False
                t = ete3.Tree(fn)
                for n in t.traverse():
                    if (not n.is_root()) and len(n.children) > 2:
                        q_updated = True
                        n.resolve_polytomy()
                if q_updated:
                    n_edit += 1
                    t.write(outfile = fn + ".bifur.tre")
    except:
        print(fn)
    print("%d trees updated" % n_edit)
    
if __name__ == "__main__":
    d_shoot = sys.argv[1]
    resolve_shoot_polytomies(d_shoot)