# KMforCSD

R packgae KMforCSD containing the algorithm and the simulated data examples from the paper Kernel Machines for Current Status Data. 
The package also contains the artificially censored data from Section 6.2.

## Example
1. Generate simulated current status data from weibull distribution (Setting 2 in paper):
```{r}
n=500 #size of dataset
data_list = weibull_data(n=500)
data = data_list[[1]]
train_ind = sample(seq_len(n), size = 0.8*n)
train_data = data[train_ind, ]
test_data = data[-train_ind, ]
```
2. Train the KM-CSD:
```{r}
sol = KMforCSD(train_data)
```
3. Predicted values on test set:
```{r}
pred = predict(sol, test_data) #prediction
decision_function = pred$prediction #estimated conditional expectation
```

