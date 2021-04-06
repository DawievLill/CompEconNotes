

# Quick script to show some of the functionality of RCall. We could do something similar with Python, Matlab and C++. Julia has packages that interact with all these languages (and more).

using RCall
using Gadfly

# This initial portion is all R code that is called by Julia 

R"""
library(ggplot2)
library(dplyr)

df <- ToothGrowth %>%
  group_by(supp, dose) %>%
  summarise(se = sd(len) / n(), len = mean(len), n = n())

ggplot(df, aes(x = dose, y = len, group = supp, color = supp)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = len - se, ymax = len + se), width = 0.1, alpha = 0.5,
  position = position_dodge(0.005)) +  
  scale_color_manual(values = c(VC = "skyblue", OJ = "orange")) + 
  labs(x = "Dose", y = "Tooth Length") 
"""

# Retrieve dataframe from R to Julia workspace

@rget df


# Use Julia commands on the dataframe as retrieved from R. 
df[:ymin] = df[:len] - df[:se]
df[:ymax] = df[:len] + df[:se]
Gadfly.plot(df, 
    x = :dose, 
    y = :len, 
    color = :supp, 
    Geom.point,
    Guide.xlabel("Dose"), 
    Guide.ylabel("Tooth Length"), 
    Guide.xticks(ticks = [0.5, 1.0, 1.5, 2.0]),
    Geom.line, 
    Geom.errorbar, 
    ymin = :ymin, 
    ymax = :ymax, 
    Scale.color_discrete_manual("orange", "skyblue"))


