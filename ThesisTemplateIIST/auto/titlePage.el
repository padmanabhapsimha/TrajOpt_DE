(TeX-add-style-hook
 "titlePage"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("report" "10pt" "a4paper" "final")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("natbib" "numbers" "sort" "authoryear") ("geometry" "margin=1in") ("hyperref" "hidelinks")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "report"
    "rep10"
    "graphicx"
    "natbib"
    "geometry"
    "hyperref"
    "titling"
    "datetime"
    "setspace"
    "nomencl"
    "times")
   (TeX-add-symbols
    "nidtitle"
    "nidname"
    "nidscode"
    "nidguide"))
 :latex)

