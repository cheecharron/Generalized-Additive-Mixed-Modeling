library(haven)
library(dplyr)
library(sqldf)
library(ggplot2)
library(mgcv)
library(MASS)
library(jtools)

setwd('R:/DDN_QUERI_IRB_Exempt/ACDC/Statisticians/SAS Data Sets/Prelim A1c/')

a1c = read_sas('a1cdata.sas7bdat') %>%
  mutate(shifttimezeropre=dayssinceenroll+365) %>%
  mutate(shifttimezeropresq=shifttimezeropre^2) %>%
  mutate(dayssinceenrollsq=dayssinceenroll^2) %>%
  mutate(shifttimezerocub=shifttimezeropresq*shifttimezeropre) %>%
  mutate(dayssinceenrollcub=dayssinceenrollsq*dayssinceenroll) %>%
  mutate(dayssincetrtend=dayssinceenroll-180) %>%
  mutate(post=ifelse(pre_post==2,1,0)) %>%
  mutate(trtend=ifelse(dayssinceenroll<180,0,1))

a1c$PatientSSN = factor(a1c$PatientSSN)
a1c$post = factor(a1c$post)

table(a1c$post, a1c$pre_post)
summary(a1c$dayssinceenroll)

# Plotting
ggplot(data=a1c, aes(x = post, y = labValue, group = DDN, color = DDN)) +
  geom_line() +
  geom_point() +
  labs(title = "A1C Values Pre and Post Treatment",
       x = "Period",
       y = "A1C Value") +
  theme_minimal()

ggplot(data=a1c, aes(x = dayssinceenroll, y = labValue)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +  # Line indicating the start of intervention
  labs(title = "A1C Values Pre and Post Treatment",
       x = "Day (Negative values before intervention; Positive values after)",
       y = "A1C Value") +
  theme_minimal()

### GAM MODEL ###

gam_model = gam(labValue ~ s(dayssinceenroll, by=post, k=10) + s(PatientSSN, bs="re"), data=a1c, method='REML')
summary(gam_model)

## need to generate estimates from the model for plotting ##

# Create a reference grid with all days of interest
new_data = sqldf('select distinct dayssinceenroll, post
                 from a1c') 
new_data$PatientSSN = rep(a1c$PatientSSN[1],nrow(new_data))
new_pre = new_data %>% filter(post==0)
new_post = new_data %>% filter(post==1)

####### PRE #########

# Predict fitted values
pred <- predict(gam_model, new_pre, se.fit = TRUE) %>% data.frame()

### generate robust sandwich standard errors ###

# Helper function to get the linear predictor matrix
get_lp <- function(model, new_data) {
  X <- predict(model, new_data, type = "lpmatrix")
  return(X)
}

# Obtain the linear predictor matrix
X <- get_lp(gam_model, new_pre)
residuals <- resid(gam_model)
sigma <- sqrt(sum(residuals^2) / gam_model$df.residual)

# Calculate the variance-covariance matrix for coefficients
# Create the robust covariance matrix using the generalized inverse
robust_vcov <- sigma^2 * ginv(t(X) %*% X)

# Calculate robust standard errors
pred$robust_se <- sqrt(diag(X %*% robust_vcov %*% t(X)))
pred$ci_lower <- pred$fit - 1.96 * pred$robust_se
pred$ci_upper <- pred$fit + 1.96 * pred$robust_se

pred = cbind(new_pre,pred)

pred_pre = pred

####### POST #########

# Predict fitted values
pred <- predict(gam_model, new_post, se.fit = TRUE) %>% data.frame()

### generate robust sandwich standard errors ###

# Obtain the linear predictor matrix
X <- get_lp(gam_model, new_post)
residuals <- resid(gam_model)
sigma <- sqrt(sum(residuals^2) / gam_model$df.residual)

# Calculate the variance-covariance matrix for coefficients
# Create the robust covariance matrix using the generalized inverse
robust_vcov <- sigma^2 * ginv(t(X) %*% X)

# Calculate robust standard errors
pred$robust_se <- sqrt(diag(X %*% robust_vcov %*% t(X)))
pred$ci_lower <- pred$fit - 1.96 * pred$robust_se
pred$ci_upper <- pred$fit + 1.96 * pred$robust_se

pred = cbind(new_post,pred)

pred_post = pred

# combine pred_pre and pred_post
pred = rbind(pred_pre,pred_post)

#dev.off()
ggplot(data=pred,aes(x=dayssinceenroll,y=fit,color=post)) +
  geom_line() +
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper, fill=post), alpha=0.2) +
  scale_x_continuous(breaks = c(-360,-270,-180,-90,0,90,180,270,360,450), 
                     labels=c(-360,-270,-180,-90,0,90,180,270,360,450)) +
  xlab("Days (before and after start of intervention)") + 
  ylab("A1C Levels") +
  theme_apa()  +
  theme(legend.position = "none")

ggsave("R:/DDN_QUERI_IRB_Exempt/ACDC/Statisticians/SAS Output/Modeled A1C values.jpeg", width=8, height=4, dpi=600)

