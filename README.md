# comp_animal_disease


## Installation

You can install the comp_animal_disease with following command.
	
	git clone https://github.com/Gyungbu/comp_animal_disease.git
 
The list of required packages for `script` is shown in the `requirements.txt` file. When you run the command below, the required packages will be downloaded. (version : `python 3.9.7`)
	
	conda create -n env_comp_animal_disease
	conda activate env_comp_animal_disease
	conda install pip  
	conda install python=3.9.7
	pip install -r ./comp_animal_disease/requirements.txt 

# comp_update_mrs : (Optional) Update Reference Data
## How to use

### 1. Prepare Merged Proportion File
Place the csv file of your Merged Proportion File in the `./comp_animal_disease/input/`?folder.

Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.
3. In the case of dog, the name of the Merged Proportion File should start with 'PD', and in the case of cat, the name of the Merged Proportion File should start with 'PC'.

### 2. Run comp_update_mrs
To run comp_update_mrs,
 
Run the command below:
  
    python ./comp_animal_disease/comp_update_mrs.py {path_exp}
    ### ex) python comp_update_mrs.py "/home/kbkim/comp_animal_disease/input/PDmirror_output_dog_1629.csv"
    ### ex) python comp_update_mrs.py "/home/kbkim/comp_animal_disease/input/PCmirror_output_cat_1520.csv"  
   
    
    
When comp_update_mrs is executed as above, the file?`db_abundance_{self.species}.xlsx`, `comp_mrs_db_{self.species}.xlsx`, `comp_percentile_rank_db_{self.species}.csv`?will be created or modified in the?`./comp_animal_disease/input/`?folder (where, {self.species} : dog or cat).
And the file `mrs_hist_{self.species}g.png` will be created or modified in the?`./comp_animal_disease/output/`?folder (where, {self.species} : dog or cat).


# comp_percentile_rank : Calculate MRS, Dysbiosis, HealthyDistance and percentile rank of samples, Benefical & Harmful microbiome abuncdance.
## How to use

### 1. Prepare Input data
Place the csv file of your proportion file in the `./comp_animal_disease/input/`?folder.
Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.
3. In the case of dog, the name of the Merged Proportion File should start with 'PD', and in the case of cat, the name of the Merged Proportion File should start with 'PC'.
4. In comp_percentile_rank_forRPA.py, under if __name__ == '__main__': enter the path of the proportion file you want to analyze in the 'path_exp' value and save it.

### 2. Run comp_percentile_rank
To run comp_percentile_rank,
 
Run the command below:

    python ./comp_animal_disease/comp_percentile_rank_forRPA.py.py 
    

When comp_percentile_rank is executed as above, the file `comp_eval_{self.species}.csv`, `comp_percentile_rank_{self.species}.csv` and `comp_scatterplot_dog.png` will be created or modified in the?`./comp_animal_disease/output/` folder.


