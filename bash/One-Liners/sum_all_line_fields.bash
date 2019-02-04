#!/bin/bash
#
awk '{ for (i = 1; i &lt;= NF; i++) s = s+$i }; END { print s+0 }' "$@"
