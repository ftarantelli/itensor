git pull
#git checkout cluster
#git checkout home

aux=''
ARGV="$@"
num="$#"

fname="$ARGV"
aux=`echo "${aux} ${fname}"`

git checkout cluster -- ${PWD}/${aux}
git commit -m "update ${aux} from cluster"
git push -u origin home
