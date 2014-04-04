(TeX-add-style-hook "manual"
 (lambda ()
    (LaTeX-add-bibliographies
     "/home/ghislain/Documents/Ghislain-CIRAD/Biblio/Biblio-These-GV")
    (LaTeX-add-labels
     "eq:binomial"
     "eq:bernoulli"
     "fig:Altitude"
     "fig:theta-binomial"
     "fig:observations-binomial"
     "fig:mcmc-binomial"
     "fig:predictions-binomial"
     "fig:pred-obs-binomial"
     "eq:siteocc")
    (TeX-add-symbols
     '("keywords" 1)
     '("mymulticolumn" 3)
     '("SetRowColor" 1)
     '("bs" 1)
     "logit"
     "p"
     "sizeBigTable"
     "newline")
    (TeX-run-style-hooks
     "colortbl"
     "xcolor"
     "booktabs"
     "longtable"
     "lineno"
     "hyperref"
     "citecolor=blue"
     "colorlinks=true"
     "setspace"
     "array"
     "inputenc"
     "utf8"
     "babel"
     "english"
     "francais"
     "amsfonts"
     "amsmath"
     "amssymb"
     "verbatim"
     "natbib"
     "sort"
     "round"
     "graphicx"
     "a4wide"
     ""
     "latex2e"
     "art12"
     "article"
     "leqno"
     "12pt"
     "a4paper")))

