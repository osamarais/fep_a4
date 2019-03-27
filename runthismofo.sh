#!/bin/bash -e
./setup.sh
./build.sh main
./build/main reorder_a.dmg reorder_a.smb
./build/main reorder_b.dmg reorder_b.smb
./build/main reorder_c.dmg reorder_c.smb
