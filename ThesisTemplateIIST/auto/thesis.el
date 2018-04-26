(TeX-add-style-hook
 "thesis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("Classes/PhDThesisPSnPDF" "a4paper" "12pt" "times" "numbered" "print" "authoryear")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "Preamble/preamble"
    "thesis-info"
    "Certificate/certificate"
    "Declaration/declaration"
    "Acknowledgement/acknowledgement"
    "Abstract/abstract"
    "Chapter1/chapter1"
    "Chapter2/chapter2"
    "Appendix1/appendix1"
    "Classes/PhDThesisPSnPDF"
    "Classes/PhDThesisPSnPDF12")
   (LaTeX-add-bibliographies
    "References/references"))
 :latex)

