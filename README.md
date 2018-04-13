# parallel-max-p-region
This is the parallel implementation of max-p region algorithm

In order to use this code, you will need pysal library. You can install pysal using the following command on windows command line.
pip install pysal

You can run this code using the following sample code max_p_region_run_example.py

Change lattice size as per your requirement w = pysal.lat2W(20,20).
Change threshold value as per your requirement floor = 100.

Since, this code uses multiprocessing pool, it can only be run from command line.
Download and save all the codes in this project into a folder.
Go to the folder where all the codes are saved from command line and use the following command to run the code.
python max_p_region_run_example.py
