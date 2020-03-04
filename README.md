To build code on UCR HPCC clusters "cluster.hpcc.ucr.edu": 
(1) git clone your repository 
(2) login to an available gpu: srun -p gpu -c 4 --gres=gpu:1 --pty bash -l 
(3) cd into Collagen_Elastin 
(4) run build_bend_model_hpcc.sh to run cmake and put files in the 'build' folder using the command ./build_model_hpcc.sh build

To run the code: 
(0) Stay logged into the gpu. 
(1) First load the modules used in build_bend_model_hpcc.sh 
(2) Run the exacutable: ./build/bend_model -dt=0.0001 data.xml
The xml file contains initial conditions for the model. 

To run from restart:
(0) use the flag 'reset' and use the xml file to initiate the model followed by the current state file. 
(1) e.g. ./build/bend_model reset -dt=0.0001 data.xml State_1.0.sta

To submit a job, use the SBATCH.sh file and run the command: 
(1) sbatch -p gpu --gres=gpu:1 --mem=10g --time=250:00:00 SBATCH.sh 
This file runs the SBATCH.sh file with a request for 1 gpu with 10gb memory for 250 hours. Note: You must be logged out of the gpu for job scheduling.

To change clusters, you will need to change the modules in build_model files. 
There are a few set up but the will go out of date soon. 
I prefer the dockerhub build since your image will always be the same.
 To use the CRC you can log into a gpu using the command "qrsh -q gpu -l gpu_card=1", and build then run. 
To use on a container based cluster, you must choose your own image (in yaml file) and then log in for building and log out to submit jobs.