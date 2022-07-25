#!/bin/bash

for i in {8..50..2}
do
	echo "====> python functions_clconf_SNR.py 10 config_felix.ini signal $i $((i+2))"
	python functions_clconf_SNR.py 1000 config_felix.ini signal $i $((i+2))

	echo "====> python functions_clconf_SNR.py 10 config_felix.ini noise $i $((i+2))"
	python functions_clconf_SNR.py 1000 config_felix.ini noise $i $((i+2))

	mv signal_add.hdf5 signal_${i}_$((i+2)).hdf5
	echo "------> file changed in signal_${i}_$((i+2)).hdf5"
	mv noise_tests.hdf5 noise_${i}_$((i+2)).hdf5
done

for i in {8..50..2}
do
	echo "====> python snr_dependency.py signal_${i}_$((i+2)).hdf5 noise_${i}_$((i+2)).hdf5 $i $((i+2)) 1000"
	python snr_dependency.py signal_${i}_$((i+2)).hdf5 noise_${i}_$((i+2)).hdf5 $i $((i+2)) 1000 100
done
