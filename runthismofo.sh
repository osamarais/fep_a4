#!/bin/bash -e
./setup.sh
./build.sh main
./build/main reorder_a.dmg reorder_a.smb
./build/main reorder_b.dmg reorder_b.smb
./build/main reorder_c.dmg reorder_c.smb
#./build/main wehshi_mesh.prt wehshi_mesh.sms
#paraview ./reorder_a.smb/reorder_a.smb.pvtu ./reorder_b.smb/reorder_b.smb.pvtu ./reorder_c.smb/reorder_c.smb.pvtu
