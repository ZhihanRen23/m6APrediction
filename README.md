
# m6APrediction

**m6APrediction** is an R package that predicts **m6A methylation sites** based on genomic feature data using a pre-trained **Random Forest** model.  
This package was developed as part of the BIO215 module (Practical 7–8) at Xi’an Jiaotong-Liverpool University.

##  Installation

You can install the package directly from GitHub using the **devtools** or **remotes** package:

```r
# If you don't have devtools installed yet:
install.packages("devtools")

# Install from your own GitHub repository:
devtools::install_github("ZhihanRen23/m6APrediction")

# Load the package
library(m6APrediction)

If you cannot access Github, you may also install from a local aechive file:

install.packages("m6APrediction_1.0.0.tar.gz", repos = NULL, type = "source")

##  Usage Example

library(m6APrediction)

# Load example model and data
miniDB <- readRDS(system.file("extdata", "m6A_model.rds", package = "m6APrediction"))
m6A_model <- read.csv(system.file("extdata", "mini_db.csv", package = "m6APrediction"))

# Run batch prediction
res <- prediction_multiple(m6A_model, head(mini_db, 10))
head(res)

# Single prediction example
prediction_single(m6A_model, mini_db[1, , drop = FALSE])

#Author
ZhihanRen23(zhihan.ren23@student.xjtlu.edu.cn)
