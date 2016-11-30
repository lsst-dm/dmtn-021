tail -n +4 index.md > index_tmp.md
pandoc index_tmp.md --mathjax --from markdown --to rst -s -o index_tmp.rst
tail -n +4 index_tmp.rst | sed '1s/^/:tocdepth: 2/' >index.rst
rpl -e ':alt:' ':name:' index.rst
make html
rm index_tmp.rst index_tmp.md
#open _build/html/index.html

