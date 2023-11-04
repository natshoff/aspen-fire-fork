


# Functions for Aim 2
calculate_ipw <- function(df,y,X) {
 # Calculate the numerator (expected distribution of aspen %)
 model_num <- lm(y ~ 1, data = df)
 num <- dnorm(y, predict(model_num), sd(model_num$residuals))
 # Now calculate the denominator (exp. distribution regressed against confounders)
 model_den <- lm(y ~ X, data=df)
 den <- dnorm(y, predict(model_den), sd(model_den$residuals))
 # Join back to the data frame
 return(df %>% mutate(ipw = num / den))
 rm(model_num,num,model_den,den)
}