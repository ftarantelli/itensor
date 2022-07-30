#git checkout cluster
git checkout home
git pull

aux=''
ARGV="$@"
num="$#"

fname="$ARGV"
aux=`echo "${aux} ${fname}"`

git checkout origin/cluster -- ${PWD}/${aux}
git commit -m "update ${aux} from cluster"
git push -u origin home
