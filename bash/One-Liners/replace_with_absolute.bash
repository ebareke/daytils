#!/bin/bash
#
awk '{ for (i = 1; i &lt;= NF; i++) if ($i &lt; 0) $i = -$i; print }' "$@"
