# Dan install packages ----

# What packages are needed
packages <- c("tidyverse", 
              "here")

# Packages that are already installed
installed_packages <- packages %in% rownames(installed.packages())

# Install all missing packages
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Above is just a fancy way to see if you have these packages are installed
# and installing the ones that aren't already there

# # This is basically equivalent to the following but it avoids installing it if 
# # it is already there 
# install.packages("tidyverse")
# install.packages("here")


