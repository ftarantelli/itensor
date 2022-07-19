aux=''
ARGV="$@"

fname="$ARGV"
aux=`echo "${aux} ${fname}"`

#echo "${aux}"

git add ${aux}
git commit -m "update files"
git push
