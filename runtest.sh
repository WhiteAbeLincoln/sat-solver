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
    h) printf "Usage: %s: [-t TESTER] [-c CANONICAL] [-su] FILES...";;
    ?) printf "Usage: %s: [-t TESTER] [-c CANONICAL] [-su] FILES..."
       exit -1;;
  esac
done
shift $((OPTIND-1))

for f in $@; do
  echo ">>>$f"
  $tester < $f
  tr="$?"
  if [[ "$mode" == "SAT" ]]; then
    if [[ "$tr" -ne "10" ]]; then
      echo "Failure, expected SAT, got $tr"
      exit -1
    fi
  elif [[ "$mode" == "UNSAT" ]]; then
    if [[ "$tr" -ne "20" ]]; then
      echo "Failure, expected UNSAT, got $tr"
      exit -1
    fi
  else
    # minisat doesn't like the trailing %\n0 that our tests have
    $canon <(sed 's/^%$\|^0$//' "$f")
    cr="$?"
    if [[ "$tr" -ne "$cr" ]]; then
      echo "Failure, expected $cr, got $tr"
      exit -1
    fi
  fi
done

# vim: ft=sh
