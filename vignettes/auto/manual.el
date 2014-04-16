(TeX-add-style-hook "manual"
 (lambda ()
    (LaTeX-add-bibliographies
     "/home/ghislain/Documents/Ghislain-CIRAD/Biblio/Biblio-These-GV")
    (LaTeX-add-labels
     "sec:binomial"
     "eq:binomial"
     "fig:Altitude"
     "eq:bernoulli"
     "fig:theta-binomial"
     "fig:observations-binomial"
     "fig:mcmc-binomial"
     "fig:predictions-binomial"
     "fig:pred-obs-binomial"
     "eq:siteocc"
     "fig:mcmc-siteocc"
     "fig:predictions-siteocc"
     "fig:predictions-siteocc-glm")
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

