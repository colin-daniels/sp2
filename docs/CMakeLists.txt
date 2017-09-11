# Sphinx targets are:
#  html        to make standalone HTML files
#  dirhtml     to make HTML files named index.html in directories
#  singlehtml  to make a single large HTML file
#  pickle      to make pickle files
#  json        to make JSON files
#  htmlhelp    to make HTML files and an HTML help project
#  qthelp      to make HTML files and a qthelp project
#  devhelp     to make HTML files and a Devhelp project
#  epub        to make an epub
#  latex       to make LaTeX files, you can set PAPER=a4 or PAPER=letter
#  latexpdf    to make LaTeX and PDF files (default pdflatex)
#  latexpdfja  to make LaTeX files and run them through platex/dvipdfmx
#  text        to make text files
#  man         to make manual pages
#  texinfo     to make Texinfo files
#  info        to make Texinfo files and run them through makeinfo
#  gettext     to make PO message catalogs
#  changes     to make an overview of all changed/added/deprecated items
#  xml         to make Docutils-native XML files
#  pseudoxml   to make pseudoxml-XML files for display purposes
#  linkcheck   to check all external links for integrity
#  doctest     to run all doctests embedded in the documentation (if enabled)
#  coverage    to run coverage check of the documentation (if enabled)
set(SPHINX_TARGET html)

set(SPHINX_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/source)
set(SPHINX_BUILD_DIR  ${CMAKE_CURRENT_BINARY_DIR}/build)
set(SPHINX_OPTIONS    "")

add_custom_target(docs
  COMMAND sphinx-build -M ${SPHINX_TARGET} "${SPHINX_SOURCE_DIR}" "${SPHINX_BUILD_DIR}" ${SPHINX_OPTIONS}
  COMMENT "Generating ${PROJECT_NAME} documentation with sphinx-build" VERBATIM)
