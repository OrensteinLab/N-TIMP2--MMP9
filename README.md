# Predicting the affinity landscape of N-TIMP2/MMP9<sub>CAT</sub> by combining deep neural networks and deep mutational scans

# Introduction

We present a novel approach using deep neural networks (DNN) trained on High-Throughput Sequencing (HTS) data for protein mutagenesis library affinity screens. 
In this study, we have focused on the experimental raw data from the N-terminal domain of the tissue inhibitor of metalloproteinases 2 (N-TIMP2) with the catalytic domain 
of matrix metalloproteinase 9 (MMP9<sub>CAT</sub>). Our goal is to comprehensively measure Protein-Protein Interactions (PPIs) and accurately predict unobserved affinity-enhancing 
or affinity-reducing variants.

# Purpose
This repository includes the code related to our research project aimed at predicting the affinity landscape of N-TIMP2/MMP9<sub>CAT</sub> by combining deep neural networks
and deep mutational scans

# Data
The experimental raw data utilized in this study is derived from the N-TIMP2/MMP9<sub>CAT</sub> complex. High-Throughput Sequencing (HTS) data has been used to train deep neural networks,
enabling us to accurately predict unobserved affinity-enhancing or affinity-reducing variants.

# Setup Environment Instructions (Windows)

Before you proceed with the setup, make sure to have Python and Anaconda installed on your system.


1. **Download the Code Repository:**
   - Visit the GitHub repository: [https://github.com/OrensteinLab/N-TIMP2--MMP9/tree/main/Code](https://github.com/OrensteinLab/N-TIMP2--MMP9/tree/main/Code)
   - Download the contents of the "Code" folder.

2. **Inside the "Code" Folder, Add Your Raw Data:**
   - For code validation, include the datasets named "All_variant_ala.csv" and "All_variant_no_ala.csv."

3. **Create a Virtual Conda Environment:**
   - Open a command prompt.
   - Navigate to the directory where you downloaded the "Code" repository.
   - Run the following command to create a virtual conda environment named "my_env" with Python 3.9.16 and the required modules:
     ```bash
     conda create --name my_env --file requirements.txt python=3.9.16
     ```

4. **Activate the New Environment and Run the Script:**
   - Activate the environment using the following command:
     ```bash
     conda activate my_env
     ```
   - Run the scripts according to the provided usage instructions.


# Usage
### 1.	Run the pre-processing script:
Execute the script `pre-proccesing_to_raw_data.py` in the presence of the raw data files located in the data folder. Upon completion, the script will generate 
two files in the folder:

•	all_variant_ala.csv: Contains all variants containing Alanine (Ala)

•	all_variant_no_ala.csv: Contains all variants without Alanine (Ala)
### 2.	Train models script:
While running the `train_model.py` script, you will be prompted to enter input representing the action you want to perform. The code supports the following three options:

1-	Use all the data to train the model without making any predictions.

2-	Split the data into three sets: a training set (80%), a validation set (10%), and a test set (10%). Train the model using the training set and make predictions 
on the validation set.

3-	Split the data into a training set (90%) and a test set (10%). Train the model using the training set and make predictions on the validation set.

•	If you choose Option 1, the trained models will be saved in the folder.

•	If you choose Option 2 or 3, prediction files will be generated at the end of the script execution.
### 3.	Get predictions:
There are two scripts available for performing predictions:

1.	Independent Dataset Prediction:
`train_inference_ki.py` script removes variants from the training dataset matching those in the external dataset. Then, trains the models and makes predictions for
these variants.

2.	All Single Mutations Prediction:
The script `inference_heatmap.py` performs predictions for all single mutations.
