pandoc index.md --mathjax --from markdown --to rst -s -o index_tmp.rst
tail -n +4 index_tmp.rst | sed '1s/^/:tocdepth: 1/' >index.rst
make html
rm index_tmp.rst
#open _build/html/index.html

