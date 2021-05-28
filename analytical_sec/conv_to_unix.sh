
OIFS="$IFS"
IFS=$'\n'

for file in `find ./data -type f -name "*.csv"`
do
  echo ${file}
  dos2unix ${file}
done

IFS="$OIFS"

