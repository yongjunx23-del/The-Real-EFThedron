"""Assemble a standalone editor-compatible appendix and manuscript fragment."""
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'research/figure_repro_20261008'
OUT = ROOT / 'output/FIGURE_REPRO_20261008'

PREAMBLE = r'''\documentclass[11pt]{article}
\usepackage[margin=0.9in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,longtable,array,booktabs,tabularx,url}
\usepackage[hidelinks]{hyperref}
\setlength{\emergencystretch}{3em}
\setlength{\parskip}{0.35em}
\setlength{\parindent}{0pt}
\allowdisplaybreaks
\title{Certified $g_0$--$g_2$ bounds in four dimensions\\
\large Constraints, optimization, certification, and reproducibility appendix}
\author{}
\date{8 October 2026}
\begin{document}
\maketitle
\begin{abstract}
This appendix specifies the inequalities underlying all certified points and the
certified polygon in the final comparison with the reference primal and dual
boundaries. It gives the complete two retained null-row banks, the witness-search
restrictions and approximations, the continuum certification methods, and the
exact-data format. The accompanying package preserves all selected multipliers,
proof dependencies and replay sources. There are 38 numerical inequalities and
one analytic inequality; 22 distinct certificates are marked in the lower panel.
\end{abstract}
\tableofcontents
\clearpage
\appendix
'''

def markdown_from_tex(s):
    s = re.sub(r'\\label\{[^}]+\}', '', s)
    s = re.sub(r'\\section\{([^}]+)\}', r'## \1', s)
    s = re.sub(r'\\subsection\{([^}]+)\}', r'### \1', s)
    s = re.sub(r'\\texttt\{([^}]+)\}', lambda m:'`'+m[1].replace(r'\_','_')+'`', s)
    s = s.replace(r'Eq.~\eqref{eq:coercivity}', 'the negative-V support inequality')
    s = s.replace(r'\begin{itemize}', '').replace(r'\end{itemize}', '')
    s = re.sub(r'^\s*\\item ', '- ', s, flags=re.M)
    s = s.replace(r'\(', '$').replace(r'\)', '$')
    s = s.replace(r'\[', '$$').replace(r'\]', '$$')
    s = s.replace(r'\begin{align}', '$$\n'+r'\begin{aligned}')
    s = s.replace(r'\end{align}', r'\end{aligned}'+'\n$$')
    s = s.replace(r'\nonumber', '')
    return s

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    names = ['constraints_appendix.tex','certification_appendix.tex',
             'geometry_appendix.tex','selected_bound_table.tex']
    body = '\n\n'.join((SRC/n).read_text(encoding='utf-8') for n in names)
    (OUT/'appendix.tex').write_text(PREAMBLE+body+'\n\\end{document}\n', encoding='utf-8')
    (OUT/'appendix_fragment.tex').write_text('% Requires amsmath, amssymb, longtable, array, booktabs.\n'
        '% Add \\appendix once in the host manuscript before this fragment.\n'+body, encoding='utf-8')
    md = '# Certified g₀–g₂ bounds: complete appendix\n\n'
    md += 'Scope: all certified points and the certified polygon in the final dual/primal comparison. '
    md += 'The exact coefficients and proof records accompany this document.\n\n'
    md += (SRC/'constraints_appendix.md').read_text(encoding='utf-8')+'\n\n'
    md += markdown_from_tex((SRC/'certification_appendix.tex').read_text(encoding='utf-8'))+'\n\n'
    md += markdown_from_tex((SRC/'geometry_appendix.tex').read_text(encoding='utf-8'))+'\n\n'
    md += (SRC/'selected_bound_table.md').read_text(encoding='utf-8')
    (OUT/'appendix.md').write_text(md, encoding='utf-8')
    for n in ['optimization_records.json','inventory_verify.py']:
        shutil.copy2(SRC/n, OUT/n)
    print('Built standalone appendix.tex, manuscript fragment, Markdown appendix and verification entry point.')

if __name__ == '__main__':
    main()
