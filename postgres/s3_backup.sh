#!/bin/bash
# first argument is the s3 object to write
set -e -o pipefail
pg_dump -h localhost -p 5001 -U postgres -Fc --compress=0 | lz4 --fast | aws s3 cp  --expected-size 800000000000  - $1
