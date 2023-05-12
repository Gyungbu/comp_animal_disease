##<Usage: python Script.py {path_exp}>
### ex) python comp_percentile_rank.py "/home/kbkim/comp_animal_disease/input/mirror_output_dog_1340.csv" 
### ex) python comp_percentile_rank.py "/home/kbkim/comp_animal_disease/input/mirror_output_cat_1409.csv"

import os, datetime
import pandas as pd
from scipy.stats import percentileofscore, pearsonr
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from skbio.stats.composition import multiplicative_replacement, clr

# path_exp : Path of Merged Proportion file to analyze
#path_exp = sys.argv[1] 

#-------------------------------------------------------
# 공통 함수
#-------------------------------------------------------
def WriteLog(functionname, msg, type='INFO', fplog=None):
    #strmsg = "[%s][%s][%s] %s\n" % (datetime.datetime.now(), type, functionname, msg)
    #if( DEBUGMODE ): 
    #    print(strmsg)
    #else:
    #    if( fplog != None ):
    #        fplog.write(strmsg)
    
    head = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    writestr = f"[{head}][{functionname}] {msg}\n"
    #if( DEBUGMODE ):
    if( True ):
        #print(message)
        writestr = f"[{functionname}] {msg}\n"
        print(writestr)
        
    if( fplog != None ):
        fplog.write(writestr)
        fplog.flush()

###################################
# MainClass
###################################
class CompAnimalDisease:
    def __init__(self, path_exp, outdir=None, fplog=None):
        """
        Initializes a CompAnimalDisease object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.
        """
        self.path_exp = path_exp
        self.species = None
        self.__fplog=fplog
        if( os.path.basename(self.path_exp).startswith('PD') ):
            self.species = 'dog'
        elif( os.path.basename(self.path_exp).startswith('PC') ):
            self.species = 'cat'

        if( self.species is None ):
            print("The species should be dog[PD] or cat[PC]")
            print("Please check the path_exp")


        curdir = os.path.dirname(os.path.abspath(__file__))
        self.path_beta = f"{curdir}/input/phenotype_microbiome_{self.species}.xlsx"
        self.path_healthy = f"{curdir}/input/healthy_profile_{self.species}.xlsx"
        self.path_mrs_db = f"{curdir}/input/comp_mrs_{self.species}.xlsx"
        
        
        ###output
        if( outdir is not None ):
            self.outdir = outdir
        else:
            self.outdir = f"{curdir}/output"

        #self.path_comp_percentile_rank_output = f"{self.outdir}/{os.path.basename(self.path_exp).replace('_report.txt','')}.csv"
        #self.path_comp_eval_output = f"{self.outdir}/{os.path.basename(self.path_exp).replace('_report.txt','_eval')}.csv"
        
        self.path_comp_percentile_rank_output = f"{curdir}/output/comp_percentile_rank_{self.species}.csv"
        self.path_comp_eval_output = f"{curdir}/output/comp_eval_{self.species}.csv"

        ##ReadDB  에서 읽어들인데이타
        self.df_beta = None
        self.df_healthy = None
        self.df_exp = None
        self.df_mrs_db = None
        self.df_exp_healthy = None
        
        self.df_mrs = None
        self.df_percentile_rank = None
        self.df_eval = None
        
        self.li_diversity = None
        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_microbiome = None


    
    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_exp : Data frame of Experimental result information - Abundance    
    def ReadDB(self):
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:
            self.df_beta = pd.read_excel(self.path_beta)
            self.df_healthy = pd.read_excel(self.path_healthy)
            self.df_exp = pd.read_csv(self.path_exp)
            self.df_mrs_db = pd.read_excel(self.path_mrs_db, index_col=0) 
            self.df_exp_healthy = pd.read_csv(self.path_exp)

            self.df_beta.rename(columns = {"Disease": "phenotype", "NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_beta = self.df_beta[["phenotype", "ncbi_name", "microbiome", "beta", "microbiome_subtract"]]
            self.df_beta['beta'] = self.df_beta['beta'].replace({'유해': 1, '유익': -1})

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            
        return rv, rvmsg

    def SubtractAbundance(self): 
        """
        Subtract the abundance for each microbiome in the df_exp.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            print(self.df_exp)
            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']):
                self.li_diversity = list(self.df_exp.iloc[0,1:]) # li_diversity : Alpha-Diversity list 
                self.df_exp = self.df_exp.iloc[2:,:]
                self.df_exp_healthy = self.df_exp_healthy.iloc[2:,:]
                            
            # li_new_sample_name : Sample name list 
            # li_phenotype : Phenotype list 
            self.li_new_sample_name = list(self.df_exp.columns)[1:]  
            self.li_phenotype = list(dict.fromkeys(self.df_beta['phenotype']))
            
            print(self.df_beta)
            
            # Subtract the abundance - df_exp
            for idx_beta, row_beta in self.df_beta.iterrows(): 
                li_micro_sub = []

                if pd.isna(row_beta['microbiome_subtract']) is False:
                    li_micro_sub = row_beta['microbiome_subtract'].split('\n')

                    for micro_sub in li_micro_sub:
                        condition = (self.df_exp.taxa == row_beta['microbiome'])
                        condition_sub = (self.df_exp.taxa == micro_sub)

                        if len(self.df_exp[condition_sub]) > 0:

                            for sample_name in self.li_new_sample_name:
                                self.df_exp.loc[condition, sample_name] -= self.df_exp[condition_sub][sample_name].values[0]
           
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Check the diversity & observed rows in the exp file or db file")
    
        return rv, rvmsg

    def CalculateMRS(self): 
        """
        Calculate the MRS (Microbiome Risk Score).

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:                
            # df_mrs : Data frame of MRS corresponding to specific phenotype and sample
            self.df_mrs = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            self.df_mrs = self.df_mrs.fillna(0) 

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    condition_phen = (self.df_beta.phenotype == self.li_phenotype[j])   
                    mrs = 0

                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows():
                        condition_micro = (self.df_exp.taxa == row_beta['microbiome'])

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance = self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0] 
                            mrs += row_beta['beta'] * math.log10(100*abundance + 1) 

                    mrs /= len(self.df_beta[condition_phen])       
                    self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = -mrs

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the CalculateMRS process")
    
        return rv, rvmsg

    def CalculateDysbiosis(self): 
        """
        Calculate the Dysbiosis.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_mrs['Dysbiosis'] = 0
            self.li_microbiome = list(dict.fromkeys(self.df_beta['microbiome']))
            
            for i in range(len(self.li_new_sample_name)):
                dysbiosis = 0
                
                for j in range(len(self.li_microbiome)):
                    condition_harmful = (self.df_beta.microbiome == self.li_microbiome[j]) & (self.df_beta.beta == 1) 
                    condition_beneficial = (self.df_beta.microbiome == self.li_microbiome[j]) & (self.df_beta.beta == -1) 
                    
                    if (len(self.df_beta[condition_harmful]) >= 1) & (len(self.df_beta[condition_beneficial]) == 0):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance = self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0] 
                            dysbiosis += math.log10(100*abundance + 1)            
                            
                    elif (len(self.df_beta[condition_harmful]) == 0) & (len(self.df_beta[condition_beneficial]) >= 1):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance = self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]  
                            dysbiosis -= math.log10(100*abundance + 1)      
                            
                self.df_mrs.loc[self.li_new_sample_name[i], 'Dysbiosis'] = -dysbiosis
                 
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the CalculateDysbiosis process")
    
        return rv, rvmsg
     
    def CalculateHealthyDistance(self): 
        """
        Calculate the Healthy Distance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """           
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_mrs['HealthyDistance'] = 0     
            self.df_mrs['Diversity'] = self.li_diversity

            np_abundance = np.array([], dtype=np.float64).reshape(0,len(self.li_new_sample_name))
            np_abundance_others = np.ones((1,len(self.li_new_sample_name)), dtype=float)
            
            # Subtract the abundance - df_exp_healthy
            for idx_healthy, row_healthy in self.df_healthy.iterrows(): 
                li_micro_sub = []
                li_micro = row_healthy['microbiome'].split('\n')
                np_abundance_temp = np.zeros((1,len(self.li_new_sample_name)), dtype=float)
                
                for micro in li_micro:
                    condition_append = (self.df_exp_healthy.taxa == micro)

                    if len(self.df_exp_healthy[condition_append]) > 0:
                        np_abundance_temp += self.df_exp_healthy[condition_append].to_numpy()[:,1:].astype(np.float64)
                        np_abundance_others -= self.df_exp_healthy[condition_append].to_numpy()[:,1:].astype(np.float64)
                
                if pd.isna(row_healthy['microbiome_subtract']) is False:
                    li_micro_sub = row_healthy['microbiome_subtract'].split('\n')

                    for micro_sub in li_micro_sub:
                        condition_sub = (self.df_exp_healthy.taxa == micro_sub)
            
                        if len(self.df_exp_healthy[condition_sub]) > 0:
                            np_abundance_temp -= self.df_exp_healthy[condition_sub].to_numpy()[:,1:].astype(np.float64)
                            np_abundance_others += self.df_exp_healthy[condition_sub].to_numpy()[:,1:].astype(np.float64)
                            
                np_abundance = np.concatenate((np_abundance,np_abundance_temp),axis=0)

            np_abundance = np.concatenate((np_abundance,np_abundance_others),axis=0)
            np_abundance = np_abundance.transpose()
            
            # Apply multiplicative replacement and CLR transformations
            np_abundance = multiplicative_replacement(np_abundance)
            np_abundance = clr(np_abundance)   

            np_healthy_abundance = self.df_healthy['RA'].to_numpy()
            np_healthy_abundance = np.append(np_healthy_abundance, 100-np_healthy_abundance.sum())

            np_healthy_abundance = clr(np_healthy_abundance)
            
            # Calculate healthy distance for each new sample
            for idx in range(len(self.li_new_sample_name)):          
                healthy_dist = np.linalg.norm(np_abundance[idx] - np_healthy_abundance)            
                self.df_mrs.loc[self.li_new_sample_name[idx], 'HealthyDistance'] = -healthy_dist
            
            # Calculate the TotalScore
            self.df_mrs['TotalScore'] = self.df_mrs['Dysbiosis'] + self.df_mrs['HealthyDistance'] + self.df_mrs['Diversity']
 
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the CalculateHealthyDistance process")
            sys.exit()
    
        return rv, rvmsg          
    
    def CalculatePercentileRank(self):
        """
        Calculate the Percentile Rank and Save the Percentile Rank data as an Csv file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:      
            # Append the Dysbiosis, HealthyDistance, Diversity, TotalRiskScore to phenotype list
            self.li_phenotype += ['Dysbiosis', 'HealthyDistance', 'Diversity', 'TotalScore']

            # Create an empty data frame with the same index and columns as the df_mrs data frame
            self.df_percentile_rank = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            # Loop through all samples and phenotypes and calculate the percentile rank
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = (percentileofscore(list(self.df_mrs_db[self.li_phenotype[j]]), self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]], kind='mean')).round(1)
                 
            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_phenotype)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]<=5, self.li_phenotype[i]] = 5.0
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]>=95, self.li_phenotype[i]] = 95.0        
        
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')

            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_csv(self.path_comp_percentile_rank_output, encoding="utf-8-sig")

            #print('Analysis Complete')         
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the CalculatePercentileRank process")
            sys.exit()
    
        return rv, rvmsg

    def EvaluatePercentileRank(self):
        """
        Evaluate based on percentile rank value and Save the Evaluation data as an Csv file

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:                  
            # Define the conditions and corresponding values
            conditions = [
                self.df_percentile_rank > 80,
                (self.df_percentile_rank > 60) & (self.df_percentile_rank <= 80),
                (self.df_percentile_rank > 40) & (self.df_percentile_rank <= 60),
                (self.df_percentile_rank > 10) & (self.df_percentile_rank <= 40),
                self.df_percentile_rank <= 10
            ]
            values = ['excellent', 'good', 'average', 'caution', 'bad']

            # Apply the conditions and values using np.where()
            self.df_eval = pd.DataFrame(np.where(conditions[0], values[0],
                                np.where(conditions[1], values[1],
                                np.where(conditions[2], values[2],
                                np.where(conditions[3], values[3], values[4])))),
                                index=self.df_percentile_rank.index, columns=self.df_percentile_rank.columns)

            # Save the output file - df_eval
            self.df_eval.to_csv(self.path_comp_eval_output, encoding="utf-8-sig")
            
            print('Analysis Complete')         
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the EvaluatePercentileRank process")
            sys.exit()
    
        return rv, rvmsg
    
####################################
# main
####################################
if __name__ == '__main__':
    
    path_exp = 'input/PDmirror_output_dog.csv'
    #path_exp = 'input/PCmirror_output_cat_1409.csv'
    
    companimal = CompAnimalDisease(path_exp)
    companimal.ReadDB()
    companimal.SubtractAbundance()
    companimal.CalculateMRS()    
    companimal.CalculateDysbiosis()    
    companimal.CalculateHealthyDistance()
    companimal.CalculatePercentileRank()
    companimal.EvaluatePercentileRank()
    