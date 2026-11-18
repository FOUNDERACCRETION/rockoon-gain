"""Checks before upload: every \\R{key} in the paper exists, every citation resolves,
results.tex is not stale, and the source has no leftover notes."""
import json, os, re, sys
import rockoon_model as rm
ROOT = os.path.dirname(rm.HERE)
tex = open(os.path.join(ROOT, 'paper', 'paper.tex')).read()
res = open(os.path.join(ROOT, 'paper', 'results.tex')).read()
keys = set(re.findall(r'res@(\w+)\\endcsname', res))
used = set(re.findall(r'\\R\{(\w+)\}', tex))
bad = sorted(used - keys)
cites = set(c.strip() for grp in re.findall(r'\\cite\{([^}]*)\}', tex) for c in grp.split(','))
bibs = set(re.findall(r'\\bibitem\{(\w+)\}', tex))
labels = set(re.findall(r'\\label\{([^}]*)\}', tex)); refs = set(re.findall(r'\\ref\{([^}]*)\}', tex))
residue = [w for w in ('TODO', 'FIXME', 'XXX', 'TBD', 'VERIFY', 'placeholder', 'illustrative', 'fill in', 'Here is', 'As an AI', "I've", 'Let me')
           if re.search(w, tex)]
comments = [l for l in tex.splitlines() if re.search(r'(?<!\\)%', l)]
R = json.load(open(os.path.join(ROOT, 'results', 'results.json')))
stale = os.path.getmtime(os.path.join(ROOT, 'results', 'results.json')) < max(
    os.path.getmtime(os.path.join(rm.HERE, f)) for f in ('rockoon_model.py', 'run_analysis.py'))
ok = True
for name, items in (('undefined result keys', bad), ('cited but not in bibliography', sorted(cites - bibs)),
                    ('in bibliography but never cited', sorted(bibs - cites)), ('undefined cross-references', sorted(refs - labels)),
                    ('residue words', residue), ('LaTeX comments', comments)):
    print('%-34s %s' % (name + ':', items if items else 'none'))
    ok = ok and not items
print('%-34s %s' % ('results older than the code:', 'YES - rerun run_analysis.py' if stale else 'no'))
print('%d result values used, %d defined, %d references' % (len(used), len(keys), len(bibs)))
sys.exit(0 if ok and not stale else 1)
