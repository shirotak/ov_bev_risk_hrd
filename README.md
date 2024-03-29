# Time-dependent changes in risk of progression during use of bevacizumab for ovarian cancer

  This page provides processed data and analysis codes to reproduce the main results of the paper.   
  Other codes for preprocessing public or restricted-access data are available from the author upon reasonable request.

---
## Analysis environment
### Python 3.8.15
- lifelines 0.26.3
- matplotlib 3.4.3
- numpy 1.20.3
- pandas 1.3.4
- scikit-learn 1.0.1
- scipy 1.7.2
### R 4.1.2
- survRM2  1.0-4
- survival 3.4-0
- RISCA 1.0.3
- dplyr 1.0.10 

---
## Contents
### Python scripts  
- ```ICON7A_analysis.ipynb```
- ```Kaplan_Meier_image_analysis.ipynb```
### R scripts
- ```RMST.R```  
R script for RMST analysis
- ```ARMST.R```  
R script for adjusted RMST analysis
- ```ARMST_function.R```
- input
- R_results
- Refs
### data
ICON7-A cohort data are summarized in ```ICON7A_clinical_data.txt```.  
Other files contain coordinates of Kaplan-Meier curves extracted from published papers.
### results　　　
Output of analysis results

---
## Citation
JAMA Netw Open. 2023 Aug 1;6(8):e2326834.  
Time-Dependent Changes in Risk of Progression During Use of Bevacizumab for Ovarian Cancer  
Shiro Takamatsu, Hidekatsu Nakai, Ken Yamaguchi, Junzo Hamanishi, Masaki Mandai, Noriomi Matsumura  
PMID: 37531111  
DOI: 10.1001/jamanetworkopen.2023.26834
