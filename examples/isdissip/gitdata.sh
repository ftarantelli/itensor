act='home'
fin='cluster'

#git checkout $fin
#git checkout $act
git pull

aux=''
ARGV="$@"
num="$#"

fname="$ARGV"
aux=`echo "${aux} ${fname}"`

git checkout origin/$fin -- ${PWD}/${aux}
git commit -m "update ${aux} from ${fin}"
git push -u origin $act
