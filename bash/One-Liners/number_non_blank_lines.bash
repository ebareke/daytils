#!/bin/bash
#
awk 'NF { $0=++a " :" $0 }; { print }'
