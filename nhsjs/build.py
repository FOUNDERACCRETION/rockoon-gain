"""Build the two NHSJS versions: standard (superscript numbers) and online (full citation in double parentheses)."""
import re, subprocess, shutil, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
t = open('nhsjs.tex').read()
defined = set(re.findall(r"\\defref\{(\w+)\}", open('refs.tex').read()))
order = []
for grp in re.findall(r"\\cite\{([^}]*)\}", t):
    for k in grp.split(','):
        k = k.strip()
        assert k in defined, 'no NHSJS-format reference for ' + k
        if k not in order:
            order.append(k)
open('bib_order.tex', 'w').write('\n'.join('\\bibitem{%s} \\showref{%s}' % (k, k) for k in order) + '\n')
for mode in ('standard', 'online'):
    open('mode.tex', 'w').write('\\onlinetrue\n' if mode == 'online' else '\\onlinefalse\n')
    for f in ('nhsjs.aux', 'nhsjs.out'):
        if os.path.exists(f):
            os.remove(f)
    for _ in range(3):
        subprocess.run(['pdflatex', '-interaction=nonstopmode', 'nhsjs.tex'], stdout=subprocess.DEVNULL)
    shutil.copy('nhsjs.pdf', 'nhsjs_%s.pdf' % mode)
    pages = int(re.search(r'Pages:\s+(\d+)', subprocess.run(['pdfinfo', 'nhsjs.pdf'], capture_output=True, text=True).stdout).group(1))
    print('%s: %d pages%s' % (mode, pages, '  ** OVER THE NHSJS 20-PAGE LIMIT **' if pages > 20 else ' (limit 20)'))
    print(mode, 'built;', len(order), 'references cited')
for f in ('nhsjs.aux', 'nhsjs.log', 'nhsjs.out', 'nhsjs.pdf'):
    if os.path.exists(f):
        os.remove(f)
