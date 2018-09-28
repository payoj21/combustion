# combustion
This repo provides a framework for statistical inverse and forward problems for
combustion kinetics with model inadequacy.

Dependencies: Queso, Antioch, MPI, Boost, GSL

Add the following directories to maintain current structure: bin, build,
time-series-data

Generate data with binary gen_data.

Run the inverse and forward problem with binary hydrogen_ip. This will create a
directory ``outputData''. Subsequent runs will append to this directory, so
rename or remove this directory if you want a clean run.


