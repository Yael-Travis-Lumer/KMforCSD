# KMforCSD

KMforCSD is an R package containing the algorithm and the simulated data examples from the paper Kernel Machines for Current Status Data. 
The package also contains the artificially censored data from Section 6.2.

## Installation
You can install the package from its [GitHub repository](https://github.com/Yael-Travis-Lumer/KMforCSD/). You first need to install the [devtools](https://github.com/r-lib/devtools) package.
```{r}
install.packages("devtools")
```
Then install KMforCSD using the `install_github` function in the [devtools](https://github.com/r-lib/devtools) package.
```{r}
library(devtools)
install_github("Yael-Travis-Lumer/KMforCSD")
```

## Example
1. Generate simulated current status data from weibull distribution (Setting 2 in paper):
```{r}
library(KMforCSD)
n=500 #size of dataset
data_list = weibull_data(n=500)
sim_data = data_list[[1]]
train_ind = sample(seq_len(n), size = 0.8*n)
train_data = sim_data[train_ind, ]
test_data = sim_data[-train_ind, ]
```
2. Train the KM-CSD:
```{r}
sol = KMforCSD(train_data)
```
3. Predicted values on test set:
```{r}
pred = predict(sol, test_data) #prediction
decision_function = pred$response #estimated conditional expectation
```

A line I wrote on my local computer
