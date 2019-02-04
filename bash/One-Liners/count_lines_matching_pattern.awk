#!/bin/bash
#
awk '/GENE/ { n++ }; END { print n+0 }' "$@"
