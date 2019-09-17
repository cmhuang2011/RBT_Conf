#!/bin/sh

# run the relaxation procedure

./oxDNA input_weak_constant
./oxDNA input_weak2
./oxDNA input_medium_constant
./oxDNA input_strong_constant

mv last_conf.dat relaxed_go.conf

./oxDNA input_strong_constantMajorMinor1
./oxDNA input_strong_constantMajorMinor2
./oxDNA input_strong_constantMajorMinor3


mv last_conf.dat relaxed_go2.conf

./oxDNA input_go_dS
mv last_conf.dat relaxed_go3.conf


./oxDNA input_goOXDNA2


#mv last_conf.dat relaxed.conf
#mv last_conf.dat relaxed_weak.conf
