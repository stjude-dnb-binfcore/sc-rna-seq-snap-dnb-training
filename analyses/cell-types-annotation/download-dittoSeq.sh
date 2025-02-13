#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Download R scripts for dittoSeq library

cd util
git clone https://github.com/dtm2451/dittoSeq.git
