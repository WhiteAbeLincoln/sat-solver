#!/usr/bin/env bash

tester=./serial/a.out
canon=minisat
mode=
while getopts "ht:c:su" opt; do
  case "$opt" in
    t) tester="$OPTARG";;
    c) canon="$OPTARG";;
    s) mode="SAT";;
    u) mode="UNSAT";;
    h) printf "Usage: %s: [-t TESTER] [-c CANONICAL] [-su] FILES...\n" "$0"
       exit 0;;
    ?) printf "Usage: %s: [-t TESTER] [-c CANONICAL] [-su] FILES...\n" "$0"
       exit -1;;
  esac
done
shift $((OPTIND-1))

pcount=0
fcount=0
failures=""
for f in $@; do
  failed=
  echo ">>>$f"
  $tester < $f
  tr="$?"
  if [[ "$mode" == "SAT" ]]; then
    if [[ "$tr" -ne "10" ]]; then
      echo "Failure, expected SAT, got $tr"
      failed="$f	$tr"
    fi
  elif [[ "$mode" == "UNSAT" ]]; then
    if [[ "$tr" -ne "20" ]]; then
      echo "Failure, expected UNSAT, got $tr"
      failed="$f	$tr"
    fi
  else
    # minisat doesn't like the trailing %\n0 that our tests have
    $canon <(sed 's/^%$\|^0$//' "$f")
    cr="$?"
    if [[ "$tr" -ne "$cr" ]]; then
      echo "Failure, expected $cr, got $tr"
      failed="$f	$tr"
    fi
  fi
  if [[ -z "$failed" ]]; then
    pcount=$((pcount+1))
  else
    fcount=$((fcount+1))
    failures="$failures:$failed"
  fi
done

echo
echo "$pcount of $# passed"
echo "$fcount of $# failed"
echo "Filename	Exit Code$failures" | tr ':' '\n'
if [[ "$fcount" -ne 0  ]]; then
  exit -1
fi
# vim: ft=sh
