## Code for "Spatially Varying Deep Functional Neural Network: Application in Large-Scale Crop Yield Prediction"

## Data
- `loc.info.RData`: County location information  
- `regdat.loc.RData`: Regional location data  
- `precipitation.RData`: Precipitation information  

Due to the large volume of corn yield data and daily temperature trajectories, please request the data link by email:  
yeonjoo.park [at] utsa.edu

### Source Code
- `FNN_source/`: Source code for Functional Neural Network (FNN) models - borrowed from: Wang, S., & Cao, G. (2025). Deep Neural Network for Functional Graphical Models Structure Learning. Journal of Computational and Graphical Statistics, 1–12. https://doi.org/10.1080/10618600.2025.2551268
- `functions_source/`: Supporting functions used across models and experiments  

### Main Script
- `CropYield_Prediction_SVD_funNet.R`: Script for running 10-fold cross-validation experiments for crop yield prediction using SVD-funNet
