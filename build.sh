pandoc index.md --from markdown --to rst -s -o index.rst
make html
open _build/html/index.html

