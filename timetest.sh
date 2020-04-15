#!/usr/bin/env bash

tester=./serial/a.out
while getopts "ht:c:su" opt; do
  case "$opt" in
    t) tester="$OPTARG";;
    h) printf "Usage: %s: [-t TESTER] [-c CANONICAL] [-su] FILES...\n" "$0"
       exit 0;;
    ?) printf "Usage: %s: [-t TESTER] [-c CANONICAL] [-su] FILES...\n" "$0"
       exit -1;;
  esac
done
shift $((OPTIND-1))

etimes=""
for f in $@; do
  echo -n ">>>$f "
  ftime=$(/usr/bin/time --quiet -f'%U' $tester < $f 2>&1 >/dev/null)
  etimes="$ftime $etimes"
  echo "$ftime"
done

echo -e "max\tmin\tmean\tmedian"
echo "$etimes" | xargs | tr ' ' '\n' | datamash max 1 min 1 mean 1 median 1
# vim: ft=sh
