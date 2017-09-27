# MixtureModel
This series of scripts accompanies the paper "When Averaging Goes Wrong: The Case for Mixture Model Estimation in Psychological Science".

Files included:
- Prior_Prob.R
Documents illustrative examples and graphs presented in the paper. To run mixture modeling on your data, use MixtureModel.R

- Mixture_Model.R
The workhorse of the paper (may take a while to run). Presents all functions and analyses used to estimate the number of mixture components and posterior probabilities. Calls the function Boot_Comp.R. To reproduce the cognitive training example presented in the paper, download PPS_MelbyLervag.csv (included here).

- Boot_Comp.R
Bootstrap function to determine the number of components of the mixture. Automatically called via the source function in Mixture_Model.R

- Plug_In.R
User-friendly version of Mixture_Model, to allow pluging-in one's data. Each step is detailed and annotated. 

Please cite as:
Moreau, D., & Corballis, M. C. (2018). When Averaging Goes Wrong: The Case for Mixture Model Estimation in Psychological Science

Suggestions or comments? Please send an email to d.moreau@auckland.ac.nz
