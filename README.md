Timing of breeding reveals a trade-off between constitutive immune
investment and life history in a migratory bird
================
David Chang van Oordt, Conor C. Taff, Thomas A. Ryan, Maren N. Vitousek

-   [Method assessment](#method-assessment)
    -   [Assay repeatability](#assay-repeatability)
    -   [Treament type comparisons](#treament-type-comparisons)
-   [Results](#results)
    -   [BKA in our population](#bka-in-our-population)
    -   [BKA and Reproductive Effort](#bka-and-reproductive-effort)
        -   [General trend evaluation](#general-trend-evaluation)
            -   [BKA and age](#bka-and-age)
            -   [BKA and Lay Date](#bka-and-lay-date)
        -   [Clutch Size](#clutch-size)
            -   [Clutch size data
                distribution](#clutch-size-data-distribution)
            -   [Clutch size models](#clutch-size-models)
        -   [Nestling Feeding Rate](#nestling-feeding-rate)
            -   [Feeding Rate variation](#feeding-rate-variation)
            -   [Feeding Rate models](#feeding-rate-models)
            -   [Including Mate
                provisioning](#including-mate-provisioning)
        -   [Mass Loss](#mass-loss)
            -   [Mass loss distribution](#mass-loss-distribution)
            -   [Mass loss models](#mass-loss-models)
    -   [BKA and reproductive success](#bka-and-reproductive-success)
        -   [Fledging success](#fledging-success)
            -   [Fleding success models](#fleding-success-models)
        -   [Number of nestlings fledged](#number-of-nestlings-fledged)
            -   [Number fledged
                distribution](#number-fledged-distribution)
            -   [Number fledged models](#number-fledged-models)
        -   [Nestling mass](#nestling-mass)
            -   [Nestling mass
                distribution](#nestling-mass-distribution)
            -   [Nestling mass models](#nestling-mass-models)
    -   [Survival](#survival)
        -   [2020 Return Rate](#2020-return-rate)
            -   [Return rate models](#return-rate-models)
-   [References](#references)

# Method assessment

## Assay repeatability

``` r
# Repeatability by plate

pseurep <- bk %>%
  ggplot() + geom_point(aes(x=rep1, y=rep2)) +
  geom_abline(slope =1) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1),
        legend.position = "none") +
  labs(title="Pseureplicates - Plate effect",
       x = "No. of CFUs in Plate 1",
       y = "No. of CFUs in Plate 2") +
  scale_y_continuous(expand = c(0,0), limits = c(0,700)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,700)) +
  annotate("text", x = 500, y = 600, 
           label = paste("r","=",
                         round(cor(cbind(bk$rep1, bk$rep2),
                                   method = "pearson")[1,2], 3) ) )

# Repeatability by reaction

reps <- inner_join(filter(bk, rxn.rep == 1) %>% select(Individual_Band, bkc),
           filter(bk, rxn.rep == 2) %>% select(Individual_Band, bkc),
           by = "Individual_Band")

rxnrep <- reps %>%
  ggplot() + geom_point(aes(x=bkc.x, y=bkc.y)) +
  geom_abline(slope =1) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1),
        legend.position = "none") +
  labs(title="Reaction repeatibility",
       x = "Prop. Bacteria Killed in Challenge 1",
       y = "Prop. Bacteria Killed in Challenge 2") +
  annotate("text", x = 0.5, y = 0.95, 
           label = paste("r =",round(cor(reps[,-1])[1,2], 3) ) ) 

grid.arrange(pseurep, rxnrep, ncol = 2)
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/repeatability-1.png" style="display: block; margin: auto;" />

## Treament type comparisons

``` r
# Boxplot
df %>% 
  mutate(timing = as.factor(
    ifelse(Experiment == "Water", "30 minutes", "<3 minutes")
  )) %>% 
  ggplot() +
  geom_boxplot(aes(x = Treatment2, y = Bacteria_Killing_Assay),
               size = 1, width = 0.25, fill = "grey") +
  geom_point(aes(x = Treatment2, y = Bacteria_Killing_Assay), 
             position = position_jitter(width = 0.1)) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1),
        legend.position = "none",
        legend.title.align = 0.5,
        legend.background = element_rect(fill=NA),
        axis.text.x=element_text(angle=20,hjust=1)) +
  scale_x_discrete(labels = c("Predator-Dull", "Predator-Control", "Long-term Control", "Control-Control", "Control-Dull")) +
  labs(x = "Timing of Blood Sample",
       y = "Prop. of Bacteria Killed")
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/treatment comparison-1.png" style="display: block; margin: auto;" />

``` r
# T-test
bka_mod1 <- lm(bk ~ Treatment2, data = df)
summary(aov(bk ~ Treatment2, data = df))
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Treatment2   4  0.158 0.03944   0.314  0.868
    ## Residuals   55  6.917 0.12576

# Results

## BKA in our population

``` r
bk %>%
  group_by(Individual_Band) %>% 
  summarise(bkc = mean(bkc)) %>%
  ggplot() + geom_histogram(aes(bkc), bins = 15, colour = "white") +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1)) +
  scale_y_continuous(limits = c(0, 10), breaks=seq(0,10,by=2),
                     expand = c(0,0)) +
  geom_vline(xintercept=0, linetype = "dashed", colour = "grey") +
  labs(x = "Prop. Bacteria Killed", y = "No. of individuals")
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/bka summary-1.png" style="display: block; margin: auto;" />

``` r
summ_bk <- group_by(bk, Individual_Band) %>% summarise(bkc = mean(bkc))
summ <- matrix(0, nrow = 1, ncol = 6)
colnames(summ) <- c("Mean", "S.D.", "2.5%", "97.5%", "Min.", "Max.")
summ[1,1] <- mean(summ_bk$bkc, na.rm = TRUE)
summ[1,2] <- sd(summ_bk$bkc, na.rm = TRUE)
summ[1,3] <- summ[1,1] - 1.96*summ[1,2]/sqrt(nrow(summ_bk))
summ[1,4] <- summ[1,1] + 1.96*summ[1,2]/sqrt(nrow(summ_bk))
summ[1,5] <- min(summ_bk)
summ[1,6] <- max(summ_bk)
kable(summ)
```

|      Mean |      S.D. |      2.5% |     97.5% |       Min. |      Max. |
|----------:|----------:|----------:|----------:|-----------:|----------:|
| 0.4850821 | 0.3462755 | 0.3974623 | 0.5727019 | -0.1838843 | 281128707 |

## BKA and Reproductive Effort

### General trend evaluation

#### BKA and age

``` r
summary(glm(bk ~ Age,  family=Gamma(link="inverse"),
            data = df))
```

    ## 
    ## Call:
    ## glm(formula = bk ~ Age, family = Gamma(link = "inverse"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.7340  -0.7319  -0.1405   0.3950   1.0823  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   1.7168     0.2433   7.056 2.36e-09 ***
    ## AgeSY         0.3976     0.3392   1.172    0.246    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.4620108)
    ## 
    ##     Null deviance: 46.235  on 59  degrees of freedom
    ## Residual deviance: 45.611  on 58  degrees of freedom
    ## AIC: 41.409
    ## 
    ## Number of Fisher Scoring iterations: 6

#### BKA and Lay Date

``` r
summary(glm(bk ~ ylaydate, family = Gamma(link="inverse"),
                  data = df))
```

    ## 
    ## Call:
    ## glm(formula = bk ~ ylaydate, family = Gamma(link = "inverse"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.7590  -0.6939  -0.1986   0.4916   1.1956  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) -5.48379    5.57535  -0.984    0.329
    ## ylaydate     0.05364    0.04041   1.328    0.190
    ## 
    ## (Dispersion parameter for Gamma family taken to be 0.4847828)
    ## 
    ##     Null deviance: 46.235  on 59  degrees of freedom
    ## Residual deviance: 45.375  on 58  degrees of freedom
    ## AIC: 41.06
    ## 
    ## Number of Fisher Scoring iterations: 6

### Clutch Size

#### Clutch size data distribution

``` r
ggplot(df) +
  geom_histogram(aes(x=Clutch_Size), colour = "white", bins = 15) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Clutch size",
       y = "No. of individuals") +
  scale_x_continuous(limits = c(0,8),
                                    breaks = 1:8)
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/clutch distribution-1.png" style="display: block; margin: auto;" />

#### Clutch size models

##### Null model

``` r
clutch_m1 <- glm(Clutch_Size ~ Age + ylaydate,
                family = quasipoisson(link = "log"), 
                data = df)
summary(clutch_m1)
```

    ## 
    ## Call:
    ## glm(formula = Clutch_Size ~ Age + ylaydate, family = quasipoisson(link = "log"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.49967  -0.20721   0.01626   0.19048   0.79329  
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.692395   0.610564   6.048  1.2e-07 ***
    ## AgeSY       -0.091341   0.039065  -2.338  0.02291 *  
    ## ylaydate    -0.014077   0.004417  -3.187  0.00234 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.1175148)
    ## 
    ##     Null deviance: 9.2551  on 59  degrees of freedom
    ## Residual deviance: 7.1867  on 57  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

##### BKA model

``` r
clutch_m2 <- glm(Clutch_Size ~ Bacteria_Killing_Assay + Age + ylaydate,
                family = quasipoisson(link = "log"), 
                data = df)
summary(clutch_m2)
```

    ## 
    ## Call:
    ## glm(formula = Clutch_Size ~ Bacteria_Killing_Assay + Age + ylaydate, 
    ##     family = quasipoisson(link = "log"), data = df)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.44530  -0.19636   0.02834   0.16755   0.76566  
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)             3.60270    0.61444   5.863 2.53e-07 ***
    ## Bacteria_Killing_Assay -0.06489    0.05703  -1.138   0.2600    
    ## AgeSY                  -0.08525    0.03941  -2.163   0.0348 *  
    ## ylaydate               -0.01323    0.00447  -2.961   0.0045 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.1173973)
    ## 
    ##     Null deviance: 9.2551  on 59  degrees of freedom
    ## Residual deviance: 7.0351  on 56  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

##### Interaction model

``` r
clutch_m3 <- glm(Clutch_Size ~ Bacteria_Killing_Assay*ylaydate + Age,
                 family = quasipoisson(link = "log"), 
                 data = df)
summary(clutch_m3)
```

    ## 
    ## Call:
    ## glm(formula = Clutch_Size ~ Bacteria_Killing_Assay * ylaydate + 
    ##     Age, family = quasipoisson(link = "log"), data = df)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -1.18629  -0.20471   0.01383   0.17487   0.74587  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                      1.984922   0.938139   2.116   0.0389 *
    ## Bacteria_Killing_Assay           3.663823   1.665632   2.200   0.0321 *
    ## ylaydate                        -0.001609   0.006769  -0.238   0.8130  
    ## AgeSY                           -0.081385   0.038846  -2.095   0.0408 *
    ## Bacteria_Killing_Assay:ylaydate -0.026740   0.011934  -2.241   0.0291 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.1133617)
    ## 
    ##     Null deviance: 9.2551  on 59  degrees of freedom
    ## Residual deviance: 6.4620  on 55  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

##### Model comparison

``` r
# redefine models to poisson models
clutch_m1p <- glm(Clutch_Size ~ Age + ylaydate,
                family = poisson(link = "log"), 
                data = df)
clutch_m2p <- glm(Clutch_Size ~ Bacteria_Killing_Assay + Age + ylaydate,
                family = poisson(link = "log"), 
                data = df)
clutch_m3p <- glm(Clutch_Size ~ Bacteria_Killing_Assay*ylaydate + Age,
                 family = poisson(link = "log"), 
                 data = df)

ICtab(clutch_m1p, clutch_m2p, clutch_m3p, dispersion = dfun(clutch_m1p), type = "qAICc")
```

    ##            dqAICc df
    ## clutch_m3p 0.0    5 
    ## clutch_m1p 1.3    3 
    ## clutch_m2p 2.4    4

### Nestling Feeding Rate

#### Feeding Rate variation

``` r
# Note that Offset refers to Nestling Age.

ggplot(rfid_d) + 
  geom_point(aes(x= Offset, y = FemFeed, group = f_rfid), alpha = 0.1, show.legend = F) + 
  geom_line(aes(x= Offset, y = FemFeed, group = f_rfid), alpha = 0.1, show.legend = F) +
  geom_boxplot(aes(x = Offset, y = FemFeed, group = as.factor(Offset), fill = as.factor(Offset)),
               width = 0.3) +
  scale_x_continuous(expand = c(0,0), limits = c(0,20)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1),
        legend.position = "none") +
  labs(title="Daily provisioning visits per female",
       x = "Nestling Age (days)",
       y = "Number of visits per day")
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/visits variation-1.png" style="display: block; margin: auto;" />

#### Feeding Rate models

``` r
scaled_pc <- pc %>%
  mutate(scBrood_Size_Hatching = scale(Brood_Size_Hatching),
         scOffset = scale(Offset),
         scHour = scale(Hour),
         scjlaydate = scale(jlaydate),
         sctemp = scale(temp),
         scMalFeed = scale(MalFeed),
         scbkc = scale(bkc),
         Individual_Band = as.factor(Individual_Band))
```

##### For all nests

###### Null model

These models are modified from (Vitousek et al. 2018)

``` r
provisioning_m1 <- glmer(data = scaled_pc, FemFeed ~ scBrood_Size_Hatching + 
                           I(scBrood_Size_Hatching^2) + scOffset + I(scOffset^2) + 
                           scHour + I(scHour^2) + scjlaydate + sctemp + 
                           Treatment + scBrood_Size_Hatching*scOffset + 
                           (1|Individual_Band),
                         family=poisson(link="log")); 
summary(provisioning_m1)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: FemFeed ~ scBrood_Size_Hatching + I(scBrood_Size_Hatching^2) +  
    ##     scOffset + I(scOffset^2) + scHour + I(scHour^2) + scjlaydate +  
    ##     sctemp + Treatment + scBrood_Size_Hatching * scOffset + (1 |  
    ##     Individual_Band)
    ##    Data: scaled_pc
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  48307.2  48403.3 -24139.6  48279.2     7027 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1653 -1.0546 -0.0786  0.9313 11.1717 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Individual_Band (Intercept) 0.01315  0.1147  
    ## Number of obs: 7041, groups:  Individual_Band, 30
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     2.524134   0.057815  43.659  < 2e-16 ***
    ## scBrood_Size_Hatching           0.047065   0.028375   1.659   0.0972 .  
    ## I(scBrood_Size_Hatching^2)     -0.021022   0.019618  -1.072   0.2839    
    ## scOffset                       -0.031015   0.005228  -5.932 2.98e-09 ***
    ## I(scOffset^2)                  -0.245156   0.005854 -41.875  < 2e-16 ***
    ## scHour                          0.052005   0.004048  12.848  < 2e-16 ***
    ## I(scHour^2)                     0.026744   0.004907   5.450 5.04e-08 ***
    ## scjlaydate                     -0.056493   0.025524  -2.213   0.0269 *  
    ## sctemp                          0.131264   0.005103  25.722  < 2e-16 ***
    ## TreatmentControl_Control       -0.019011   0.069585  -0.273   0.7847    
    ## TreatmentPredator_Control       0.053191   0.063692   0.835   0.4036    
    ## TreatmentControl_Dull           0.143297   0.077548   1.848   0.0646 .  
    ## scBrood_Size_Hatching:scOffset -0.037268   0.004550  -8.191 2.58e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 13 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

###### BKA model

``` r
# bkc model
provisioning_m2 <- glmer(data = scaled_pc, FemFeed ~ scBrood_Size_Hatching + 
                           I(scBrood_Size_Hatching^2) + scOffset + I(scOffset^2) + 
                           scHour + I(scHour^2) + scjlaydate + sctemp + 
                           Treatment + scBrood_Size_Hatching*scOffset + 
                           (1|Individual_Band) + scbkc,
                         family=poisson(link="log"))
summary(provisioning_m2)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: FemFeed ~ scBrood_Size_Hatching + I(scBrood_Size_Hatching^2) +  
    ##     scOffset + I(scOffset^2) + scHour + I(scHour^2) + scjlaydate +  
    ##     sctemp + Treatment + scBrood_Size_Hatching * scOffset + (1 |  
    ##     Individual_Band) + scbkc
    ##    Data: scaled_pc
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  48303.7  48406.6 -24136.9  48273.7     7026 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1613 -1.0541 -0.0798  0.9329 11.1854 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Individual_Band (Intercept) 0.01077  0.1038  
    ## Number of obs: 7041, groups:  Individual_Band, 30
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     2.564641   0.055088  46.556  < 2e-16 ***
    ## scBrood_Size_Hatching           0.064507   0.026742   2.412   0.0159 *  
    ## I(scBrood_Size_Hatching^2)     -0.030332   0.018217  -1.665   0.0959 .  
    ## scOffset                       -0.031120   0.005225  -5.956 2.59e-09 ***
    ## I(scOffset^2)                  -0.245177   0.005851 -41.904  < 2e-16 ***
    ## scHour                          0.052067   0.004048  12.863  < 2e-16 ***
    ## I(scHour^2)                     0.026681   0.004907   5.438 5.40e-08 ***
    ## scjlaydate                     -0.032876   0.025109  -1.309   0.1904    
    ## sctemp                          0.131151   0.005102  25.704  < 2e-16 ***
    ## TreatmentControl_Control       -0.031373   0.063541  -0.494   0.6215    
    ## TreatmentPredator_Control       0.012457   0.060240   0.207   0.8362    
    ## TreatmentControl_Dull           0.104120   0.072188   1.442   0.1492    
    ## scbkc                          -0.059141   0.023948  -2.470   0.0135 *  
    ## scBrood_Size_Hatching:scOffset -0.037252   0.004548  -8.191 2.60e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 14 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

###### Interaction model

``` r
provisioning_m3 <- glmer(data = scaled_pc, FemFeed ~ scBrood_Size_Hatching + 
                            I(scBrood_Size_Hatching^2) + scOffset + I(scOffset^2) + 
                            scHour + I(scHour^2) + sctemp + 
                            Treatment + scBrood_Size_Hatching*scOffset + 
                            (1|Individual_Band) + scbkc*scjlaydate,
                          family=poisson(link="log"))
summary(provisioning_m3)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: FemFeed ~ scBrood_Size_Hatching + I(scBrood_Size_Hatching^2) +  
    ##     scOffset + I(scOffset^2) + scHour + I(scHour^2) + sctemp +  
    ##     Treatment + scBrood_Size_Hatching * scOffset + (1 | Individual_Band) +  
    ##     scbkc * scjlaydate
    ##    Data: scaled_pc
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  48305.3  48415.0 -24136.6  48273.3     7025 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1614 -1.0554 -0.0808  0.9332 11.1836 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Individual_Band (Intercept) 0.01056  0.1028  
    ## Number of obs: 7041, groups:  Individual_Band, 30
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     2.578049   0.057955  44.483  < 2e-16 ***
    ## scBrood_Size_Hatching           0.075762   0.031161   2.431  0.01504 *  
    ## I(scBrood_Size_Hatching^2)     -0.035262   0.019424  -1.815  0.06947 .  
    ## scOffset                       -0.031102   0.005225  -5.953 2.64e-09 ***
    ## I(scOffset^2)                  -0.245120   0.005851 -41.894  < 2e-16 ***
    ## scHour                          0.052085   0.004048  12.867  < 2e-16 ***
    ## I(scHour^2)                     0.026652   0.004907   5.432 5.59e-08 ***
    ## sctemp                          0.131100   0.005103  25.692  < 2e-16 ***
    ## TreatmentControl_Control       -0.042493   0.065028  -0.653  0.51346    
    ## TreatmentPredator_Control      -0.004650   0.064672  -0.072  0.94268    
    ## TreatmentControl_Dull           0.075564   0.082707   0.914  0.36091    
    ## scbkc                          -0.065171   0.025295  -2.576  0.00998 ** 
    ## scjlaydate                     -0.030236   0.025176  -1.201  0.22977    
    ## scBrood_Size_Hatching:scOffset -0.037235   0.004548  -8.187 2.67e-16 ***
    ## scbkc:scjlaydate                0.019781   0.028796   0.687  0.49212    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 15 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

###### Model comparison

``` r
ICtab(provisioning_m1, provisioning_m2, provisioning_m3, type = "AICc")
```

    ##                 dAICc df
    ## provisioning_m2  0.0  15
    ## provisioning_m3  1.5  16
    ## provisioning_m1  3.5  14

#### Including Mate provisioning

##### Null model

``` r
provisioning_m4 <- glmer(data = scaled_pc, FemFeed ~ scBrood_Size_Hatching + 
                           I(scBrood_Size_Hatching^2) + scOffset + I(scOffset^2) + 
                           scHour + I(scHour^2) + scjlaydate + sctemp + 
                           Treatment + scMalFeed + scBrood_Size_Hatching*scOffset + 
                           (1|Individual_Band),
                         family=poisson(link="log"),
                         control = glmerControl(optimizer = "Nelder_Mead")); 
summary(provisioning_m4)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: FemFeed ~ scBrood_Size_Hatching + I(scBrood_Size_Hatching^2) +  
    ##     scOffset + I(scOffset^2) + scHour + I(scHour^2) + scjlaydate +  
    ##     sctemp + Treatment + scMalFeed + scBrood_Size_Hatching *  
    ##     scOffset + (1 | Individual_Band)
    ##    Data: scaled_pc
    ## Control: glmerControl(optimizer = "Nelder_Mead")
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  28888.3  28984.7 -14429.2  28858.3     4563 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.4972 -0.9298 -0.0864  0.7842 10.6191 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Individual_Band (Intercept) 0.04206  0.2051  
    ## Number of obs: 4578, groups:  Individual_Band, 24
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     2.200613   0.157950  13.932  < 2e-16 ***
    ## scBrood_Size_Hatching          -0.032378   0.055083  -0.588   0.5567    
    ## I(scBrood_Size_Hatching^2)      0.032780   0.042568   0.770   0.4413    
    ## scOffset                       -0.049269   0.006660  -7.398 1.39e-13 ***
    ## I(scOffset^2)                  -0.197602   0.007677 -25.741  < 2e-16 ***
    ## scHour                          0.040522   0.005094   7.955 1.79e-15 ***
    ## I(scHour^2)                     0.040509   0.006509   6.224 4.85e-10 ***
    ## scjlaydate                     -0.063380   0.052483  -1.208   0.2272    
    ## sctemp                          0.059061   0.006804   8.680  < 2e-16 ***
    ## TreatmentControl_Control        0.202587   0.165512   1.224   0.2210    
    ## TreatmentPredator_Control       0.309934   0.148464   2.088   0.0368 *  
    ## TreatmentControl_Dull           0.419452   0.181846   2.307   0.0211 *  
    ## scMalFeed                       0.209892   0.004665  44.990  < 2e-16 ***
    ## scBrood_Size_Hatching:scOffset -0.049505   0.005655  -8.754  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 14 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

##### BKA model

``` r
provisioning_m5 <- glmer(data = scaled_pc, FemFeed ~ scBrood_Size_Hatching + 
                           I(scBrood_Size_Hatching^2) + scOffset + I(scOffset^2) + 
                           scHour + I(scHour^2) + scjlaydate + sctemp + 
                           Treatment + scMalFeed + scBrood_Size_Hatching*scOffset + 
                           (1|Individual_Band) + scbkc,
                         family=poisson(link="log"))
summary(provisioning_m5)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: FemFeed ~ scBrood_Size_Hatching + I(scBrood_Size_Hatching^2) +  
    ##     scOffset + I(scOffset^2) + scHour + I(scHour^2) + scjlaydate +  
    ##     sctemp + Treatment + scMalFeed + scBrood_Size_Hatching *  
    ##     scOffset + (1 | Individual_Band) + scbkc
    ##    Data: scaled_pc
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  28889.4  28992.2 -14428.7  28857.4     4562 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.4959 -0.9303 -0.0859  0.7825 10.6228 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Individual_Band (Intercept) 0.0414   0.2035  
    ## Number of obs: 4578, groups:  Individual_Band, 24
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     2.220637   0.155582  14.273  < 2e-16 ***
    ## scBrood_Size_Hatching          -0.017043   0.056300  -0.303   0.7621    
    ## I(scBrood_Size_Hatching^2)      0.029265   0.042065   0.696   0.4866    
    ## scOffset                       -0.049326   0.006660  -7.406 1.30e-13 ***
    ## I(scOffset^2)                  -0.197694   0.007677 -25.752  < 2e-16 ***
    ## scHour                          0.040540   0.005094   7.959 1.74e-15 ***
    ## I(scHour^2)                     0.040528   0.006509   6.227 4.76e-10 ***
    ## scjlaydate                     -0.035004   0.061328  -0.571   0.5682    
    ## sctemp                          0.059073   0.006804   8.682  < 2e-16 ***
    ## TreatmentControl_Control        0.204987   0.164892   1.243   0.2138    
    ## TreatmentPredator_Control       0.271618   0.149405   1.818   0.0691 .  
    ## TreatmentControl_Dull           0.382574   0.180656   2.118   0.0342 *  
    ## scMalFeed                       0.209838   0.004665  44.980  < 2e-16 ***
    ## scbkc                          -0.053107   0.055736  -0.953   0.3407    
    ## scBrood_Size_Hatching:scOffset -0.049535   0.005655  -8.759  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 15 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

##### Interaction model

``` r
provisioning_m6 <- glmer(data = scaled_pc, FemFeed ~ scBrood_Size_Hatching + 
                           I(scBrood_Size_Hatching^2) + scOffset + I(scOffset^2) + 
                           scHour + I(scHour^2) + sctemp + 
                           Treatment + scBrood_Size_Hatching*scOffset + scMalFeed +
                           (1|Individual_Band) + scbkc*scjlaydate,
                         family=poisson(link="log"))
summary(provisioning_m6)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: FemFeed ~ scBrood_Size_Hatching + I(scBrood_Size_Hatching^2) +  
    ##     scOffset + I(scOffset^2) + scHour + I(scHour^2) + sctemp +  
    ##     Treatment + scBrood_Size_Hatching * scOffset + scMalFeed +  
    ##     (1 | Individual_Band) + scbkc * scjlaydate
    ##    Data: scaled_pc
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  28890.9  29000.2 -14428.4  28856.9     4561 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.4960 -0.9297 -0.0853  0.7839 10.6226 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Individual_Band (Intercept) 0.03445  0.1856  
    ## Number of obs: 4578, groups:  Individual_Band, 24
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                     2.212220   0.153355  14.425  < 2e-16 ***
    ## scBrood_Size_Hatching          -0.055482   0.074168  -0.748   0.4544    
    ## I(scBrood_Size_Hatching^2)      0.033762   0.040546   0.833   0.4050    
    ## scOffset                       -0.049339   0.006659  -7.410 1.27e-13 ***
    ## I(scOffset^2)                  -0.197590   0.007679 -25.732  < 2e-16 ***
    ## scHour                          0.040513   0.005094   7.953 1.81e-15 ***
    ## I(scHour^2)                     0.040544   0.006509   6.229 4.69e-10 ***
    ## sctemp                          0.059079   0.006804   8.683  < 2e-16 ***
    ## TreatmentControl_Control        0.213560   0.159295   1.341   0.1800    
    ## TreatmentPredator_Control       0.323801   0.160657   2.015   0.0439 *  
    ## TreatmentControl_Dull           0.472971   0.216946   2.180   0.0292 *  
    ## scMalFeed                       0.209876   0.004665  44.990  < 2e-16 ***
    ## scbkc                          -0.015699   0.070870  -0.222   0.8247    
    ## scjlaydate                     -0.067626   0.070531  -0.959   0.3377    
    ## scBrood_Size_Hatching:scOffset -0.049562   0.005655  -8.765  < 2e-16 ***
    ## scbkc:scjlaydate               -0.055049   0.074802  -0.736   0.4618    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## 
    ## Correlation matrix not shown by default, as p = 16 > 12.
    ## Use print(x, correlation=TRUE)  or
    ##     vcov(x)        if you need it

###### Model Comparison

``` r
ICtab(provisioning_m4, provisioning_m5, provisioning_m6, type = "AIC")
```

    ##                 dAIC df
    ## provisioning_m4  0.0 15
    ## provisioning_m5  1.0 16
    ## provisioning_m6  2.6 17

### Mass Loss

#### Mass loss distribution

``` r
ggplot(df) +
  geom_histogram(aes(x=prop.ml),
                 colour = "white", bins = 15) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Mass lost (g)",
       y = "No. of individuals")
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/mass loss histogram-1.png" style="display: block; margin: auto;" />

#### Mass loss models

##### Null model

``` r
prop_mass_loss_m1 <- glm(prop.ml ~ ylaydate + Brood_Size_Hatching + Treatment2, 
                         family = gaussian(link = "identity"),
                         data = df)
summary(prop_mass_loss_m1)
```

    ## 
    ## Call:
    ## glm(formula = prop.ml ~ ylaydate + Brood_Size_Hatching + Treatment2, 
    ##     family = gaussian(link = "identity"), data = df)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -0.113319  -0.033811   0.002857   0.024224   0.167666  
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                -0.2608728  0.3981631  -0.655    0.518
    ## ylaydate                    0.0008605  0.0027244   0.316    0.755
    ## Brood_Size_Hatching        -0.0048227  0.0127106  -0.379    0.707
    ## Treatment2Predator_Control  0.0073017  0.0312655   0.234    0.817
    ## Treatment2Control_Control   0.0387403  0.0331994   1.167    0.254
    ## Treatment2Control_Dull     -0.0549700  0.0380404  -1.445    0.160
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.004192246)
    ## 
    ##     Null deviance: 0.13787  on 31  degrees of freedom
    ## Residual deviance: 0.10900  on 26  degrees of freedom
    ##   (28 observations deleted due to missingness)
    ## AIC: -77.017
    ## 
    ## Number of Fisher Scoring iterations: 2

##### BKA model

``` r
prop_mass_loss_m2 <- glm(prop.ml ~ Bacteria_Killing_Assay + ylaydate + Brood_Size_Hatching + Treatment2, 
                         family = gaussian(link = "identity"),
                         data = df)
summary(prop_mass_loss_m2)
```

    ## 
    ## Call:
    ## glm(formula = prop.ml ~ Bacteria_Killing_Assay + ylaydate + Brood_Size_Hatching + 
    ##     Treatment2, family = gaussian(link = "identity"), data = df)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -0.11390  -0.03648   0.01116   0.02438   0.14134  
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                -0.415451   0.397010  -1.046    0.305
    ## Bacteria_Killing_Assay     -0.055859   0.033967  -1.644    0.113
    ## ylaydate                    0.002187   0.002760   0.792    0.436
    ## Brood_Size_Hatching        -0.004754   0.012314  -0.386    0.703
    ## Treatment2Predator_Control -0.001435   0.030751  -0.047    0.963
    ## Treatment2Control_Control   0.044602   0.032359   1.378    0.180
    ## Treatment2Control_Dull     -0.057877   0.036894  -1.569    0.129
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.003934348)
    ## 
    ##     Null deviance: 0.137868  on 31  degrees of freedom
    ## Residual deviance: 0.098359  on 25  degrees of freedom
    ##   (28 observations deleted due to missingness)
    ## AIC: -78.304
    ## 
    ## Number of Fisher Scoring iterations: 2

##### Interaction model

``` r
prop_mass_loss_m3 <- glm(prop.ml ~ Bacteria_Killing_Assay*ylaydate + Brood_Size_Hatching + Treatment2, 
                         family = gaussian(link = "identity"),
                         data = df)
summary(prop_mass_loss_m3)
```

    ## 
    ## Call:
    ## glm(formula = prop.ml ~ Bacteria_Killing_Assay * ylaydate + Brood_Size_Hatching + 
    ##     Treatment2, family = gaussian(link = "identity"), data = df)
    ## 
    ## Deviance Residuals: 
    ##       Min         1Q     Median         3Q        Max  
    ## -0.111031  -0.038201   0.008797   0.029682   0.119833  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                     -0.852603   0.599812  -1.421    0.168
    ## Bacteria_Killing_Assay           1.231646   1.323565   0.931    0.361
    ## ylaydate                         0.005564   0.004436   1.254    0.222
    ## Brood_Size_Hatching             -0.012628   0.014745  -0.856    0.400
    ## Treatment2Predator_Control       0.011692   0.033610   0.348    0.731
    ## Treatment2Control_Control        0.050271   0.032913   1.527    0.140
    ## Treatment2Control_Dull          -0.039845   0.041321  -0.964    0.345
    ## Bacteria_Killing_Assay:ylaydate -0.009272   0.009529  -0.973    0.340
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.003942726)
    ## 
    ##     Null deviance: 0.137868  on 31  degrees of freedom
    ## Residual deviance: 0.094625  on 24  degrees of freedom
    ##   (28 observations deleted due to missingness)
    ## AIC: -77.542
    ## 
    ## Number of Fisher Scoring iterations: 2

##### Model comparison

``` r
ICtab(prop_mass_loss_m1, prop_mass_loss_m2, prop_mass_loss_m3, type = "AICc")
```

    ##                   dAICc df
    ## prop_mass_loss_m1 0.0   7 
    ## prop_mass_loss_m2 0.3   8 
    ## prop_mass_loss_m3 3.0   9

## BKA and reproductive success

### Fledging success

#### Fleding success models

##### Null model

``` r
success_m1 <- glm(success ~ Brood_Size_Hatching + ylaydate + Treatment2, 
                  family = binomial(link = "logit"), 
                  data = df )
summary(success_m1)
```

    ## 
    ## Call:
    ## glm(formula = success ~ Brood_Size_Hatching + ylaydate + Treatment2, 
    ##     family = binomial(link = "logit"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8744  -1.0280   0.6686   0.8348   1.3522  
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)                 3.982e-01  1.164e+01   0.034    0.973
    ## Brood_Size_Hatching        -5.352e-02  3.159e-01  -0.169    0.865
    ## ylaydate                    2.454e-03  8.129e-02   0.030    0.976
    ## Treatment2Predator_Control -8.252e-01  9.444e-01  -0.874    0.382
    ## Treatment2Water             9.152e-01  9.546e-01   0.959    0.338
    ## Treatment2Control_Control   3.814e-01  1.038e+00   0.368    0.713
    ## Treatment2Control_Dull      1.710e+01  1.769e+03   0.010    0.992
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 68.021  on 54  degrees of freedom
    ## Residual deviance: 59.086  on 48  degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## AIC: 73.086
    ## 
    ## Number of Fisher Scoring iterations: 16

##### BKA model

``` r
success_m2 <- glm(success ~ Bacteria_Killing_Assay + Brood_Size_Hatching + 
                    ylaydate + Treatment2, family = binomial(link = "logit"), 
                  data = df )
summary(success_m2)
```

    ## 
    ## Call:
    ## glm(formula = success ~ Bacteria_Killing_Assay + Brood_Size_Hatching + 
    ##     ylaydate + Treatment2, family = binomial(link = "logit"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8208  -0.9654   0.5924   0.8283   1.4999  
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)                  -1.11606   11.87824  -0.094    0.925
    ## Bacteria_Killing_Assay       -0.90123    1.02150  -0.882    0.378
    ## Brood_Size_Hatching          -0.03889    0.31814  -0.122    0.903
    ## ylaydate                      0.01649    0.08358   0.197    0.844
    ## Treatment2Predator_Control   -0.94767    0.97092  -0.976    0.329
    ## Treatment2Water               0.83899    0.96750   0.867    0.386
    ## Treatment2Control_Control     0.39138    1.05041   0.373    0.709
    ## Treatment2Control_Dull       17.04807 1745.28325   0.010    0.992
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 68.021  on 54  degrees of freedom
    ## Residual deviance: 58.281  on 47  degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## AIC: 74.281
    ## 
    ## Number of Fisher Scoring iterations: 16

##### Interaction model

``` r
success_m3 <- glm(Success ~ Bacteria_Killing_Assay*ylaydate + 
                    Brood_Size_Hatching + Treatment2, 
                  family = binomial(link = "logit"), 
                  data = df %>%
                    mutate(Success = ifelse(Nest_Fate == "Fledged", 1, 0)))
summary(success_m3)
```

    ## 
    ## Call:
    ## glm(formula = Success ~ Bacteria_Killing_Assay * ylaydate + Brood_Size_Hatching + 
    ##     Treatment2, family = binomial(link = "logit"), data = df %>% 
    ##     mutate(Success = ifelse(Nest_Fate == "Fledged", 1, 0)))
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.8209  -0.9651   0.5924   0.8285   1.4993  
    ## 
    ## Coefficients:
    ##                                   Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)                     -1.053e+00  1.990e+01  -0.053    0.958
    ## Bacteria_Killing_Assay          -1.060e+00  4.049e+01  -0.026    0.979
    ## ylaydate                         1.602e-02  1.468e-01   0.109    0.913
    ## Brood_Size_Hatching             -3.812e-02  3.729e-01  -0.102    0.919
    ## Treatment2Predator_Control      -9.489e-01  1.017e+00  -0.933    0.351
    ## Treatment2Water                  8.380e-01  1.001e+00   0.837    0.402
    ## Treatment2Control_Control        3.908e-01  1.062e+00   0.368    0.713
    ## Treatment2Control_Dull           1.705e+01  1.745e+03   0.010    0.992
    ## Bacteria_Killing_Assay:ylaydate  1.143e-03  2.915e-01   0.004    0.997
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 68.021  on 54  degrees of freedom
    ## Residual deviance: 58.281  on 46  degrees of freedom
    ##   (5 observations deleted due to missingness)
    ## AIC: 76.281
    ## 
    ## Number of Fisher Scoring iterations: 16

##### Model comparison

``` r
ICtab(success_m1, success_m2, success_m3, type = "AICc")
```

    ##            dAICc df
    ## success_m1 0.0   7 
    ## success_m2 1.9   8 
    ## success_m3 4.8   9

### Number of nestlings fledged

#### Number fledged distribution

``` r
ggplot(filter(df, .fledged > 0 )) +
  geom_histogram(aes(x=.fledged),
                 colour = "white", bins = 15) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Number of nestlings fledged",
       y = "No. of individuals")
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/fledged histogram-1.png" style="display: block; margin: auto;" />

#### Number fledged models

##### Null model

``` r
fledged_m1 <- glm(.fledged ~ Brood_Size_Hatching + ylaydate +  Treatment2, 
                  family = quasipoisson(link = "log"), 
                  data = filter(df, .fledged > 0 ) )
summary(fledged_m1)
```

    ## 
    ## Call:
    ## glm(formula = .fledged ~ Brood_Size_Hatching + ylaydate + Treatment2, 
    ##     family = quasipoisson(link = "log"), data = filter(df, .fledged > 
    ##         0))
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5411  -0.2850   0.1599   0.4132   0.7797  
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                -2.38563    1.90587  -1.252   0.2197  
    ## Brood_Size_Hatching         0.14183    0.05698   2.489   0.0182 *
    ## ylaydate                    0.02180    0.01326   1.645   0.1098  
    ## Treatment2Predator_Control  0.06574    0.21351   0.308   0.7602  
    ## Treatment2Water             0.08246    0.20080   0.411   0.6841  
    ## Treatment2Control_Control   0.07124    0.21436   0.332   0.7418  
    ## Treatment2Control_Dull     -0.11480    0.23288  -0.493   0.6254  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.4003495)
    ## 
    ##     Null deviance: 19.605  on 38  degrees of freedom
    ## Residual deviance: 14.525  on 32  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

##### BKA model

``` r
fledged_m2 <- glm(.fledged ~ Bacteria_Killing_Assay + Brood_Size_Hatching + 
                    ylaydate + Treatment2, family = quasipoisson(link = "log"), 
                  data = filter(df, .fledged > 0 ) )
summary(fledged_m2)
```

    ## 
    ## Call:
    ## glm(formula = .fledged ~ Bacteria_Killing_Assay + Brood_Size_Hatching + 
    ##     ylaydate + Treatment2, family = quasipoisson(link = "log"), 
    ##     data = filter(df, .fledged > 0))
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.6369  -0.2827   0.1686   0.3872   0.8323  
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                -2.741354   1.938648  -1.414   0.1673  
    ## Bacteria_Killing_Assay     -0.153011   0.160231  -0.955   0.3470  
    ## Brood_Size_Hatching         0.147631   0.057386   2.573   0.0151 *
    ## ylaydate                    0.024865   0.013609   1.827   0.0773 .
    ## Treatment2Predator_Control  0.007871   0.222796   0.035   0.9720  
    ## Treatment2Water             0.043265   0.205975   0.210   0.8350  
    ## Treatment2Control_Control   0.066107   0.215119   0.307   0.7607  
    ## Treatment2Control_Dull     -0.141857   0.235896  -0.601   0.5520  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.4011574)
    ## 
    ##     Null deviance: 19.605  on 38  degrees of freedom
    ## Residual deviance: 14.159  on 31  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

##### Interaction model

``` r
fledged_m3 <- glm(.fledged ~ Bacteria_Killing_Assay*ylaydate + 
                    Brood_Size_Hatching + Treatment2, 
                  family = quasipoisson(link = "log"), 
                  data = filter(df, .fledged > 0 ) )
summary(fledged_m3)
```

    ## 
    ## Call:
    ## glm(formula = .fledged ~ Bacteria_Killing_Assay * ylaydate + 
    ##     Brood_Size_Hatching + Treatment2, family = quasipoisson(link = "log"), 
    ##     data = filter(df, .fledged > 0))
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5004  -0.3345   0.1914   0.3913   0.7919  
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                     -4.60109    2.97687  -1.546   0.1327  
    ## Bacteria_Killing_Assay           4.82774    5.92646   0.815   0.4217  
    ## ylaydate                         0.03865    0.02153   1.795   0.0827 .
    ## Brood_Size_Hatching              0.13058    0.06126   2.132   0.0413 *
    ## Treatment2Predator_Control       0.04842    0.23090   0.210   0.8353  
    ## Treatment2Water                  0.07023    0.20965   0.335   0.7400  
    ## Treatment2Control_Control        0.07624    0.21742   0.351   0.7283  
    ## Treatment2Control_Dull          -0.07905    0.24872  -0.318   0.7528  
    ## Bacteria_Killing_Assay:ylaydate -0.03572    0.04248  -0.841   0.4071  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for quasipoisson family taken to be 0.4100625)
    ## 
    ##     Null deviance: 19.605  on 38  degrees of freedom
    ## Residual deviance: 13.869  on 30  degrees of freedom
    ## AIC: NA
    ## 
    ## Number of Fisher Scoring iterations: 4

##### Model comparison

``` r
# redefine models to poisson
fledged_m1p <- glm(.fledged ~ Brood_Size_Hatching + ylaydate +  Treatment, 
                  family = poisson(link = "log"), 
                  data = filter(df, .fledged > 0 ) )
fledged_m2p <- glm(.fledged ~ Bacteria_Killing_Assay + Brood_Size_Hatching + 
                    ylaydate + Treatment, family = poisson(link = "log"), 
                  data = filter(df, .fledged > 0 ) )
fledged_m3p <- glm(.fledged ~ Bacteria_Killing_Assay*ylaydate + 
                    Brood_Size_Hatching + Treatment, 
                  family = poisson(link = "log"), 
                  data = filter(df, .fledged > 0 ) )

ICtab(fledged_m1p, fledged_m2p, fledged_m3p, dispersion = dfun(fledged_m1p), type = "qAICc")
```

    ##             dqAICc df
    ## fledged_m1p 0.0    5 
    ## fledged_m2p 2.2    6 
    ## fledged_m3p 4.1    7

### Nestling mass

#### Nestling mass distribution

``` r
ggplot(filter(nd, Nestling_Fate == "Fledged")) +
  geom_histogram(aes(x = Nestling_Mass),
                 colour = "white", bins = 15) +
  theme(panel.background = element_rect(fill=NA),
        axis.line = element_line(size=1)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x= "Nestling mass at day 12 (g)",
       y = "No. of individuals")
```

<img src="BKA_laydate_tradeoffs_files/figure-gfm/nestling mass histogram-1.png" style="display: block; margin: auto;" />

#### Nestling mass models

##### Null model

``` r
nestling_m1 <- lmer(Nestling_Mass ~ Brood_Size_Hatching + 
                     ylaydate + Treatment2 + temp + (1|Unit_Box), 
                   data = filter(nd, 
                                 Nestling_Fate == "Fledged",
                                 Bacteria_Killing_Assay != "NA") )
summary(nestling_m1)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Nestling_Mass ~ Brood_Size_Hatching + ylaydate + Treatment2 +  
    ##     temp + (1 | Unit_Box)
    ##    Data: filter(nd, Nestling_Fate == "Fledged", Bacteria_Killing_Assay !=  
    ##     "NA")
    ## 
    ## REML criterion at convergence: 562.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.90108 -0.50326  0.06341  0.49698  3.05004 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Unit_Box (Intercept) 2.146    1.465   
    ##  Residual             3.289    1.813   
    ## Number of obs: 132, groups:  Unit_Box, 34
    ## 
    ## Fixed effects:
    ##                            Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)                18.04627   10.41281 23.96477   1.733 0.095932 .  
    ## Brood_Size_Hatching        -1.34529    0.34944 22.64691  -3.850 0.000834 ***
    ## ylaydate                   -0.09624    0.11263 24.17820  -0.854 0.401246    
    ## Treatment2Control_Dull     -0.91266    1.17312 23.19759  -0.778 0.444448    
    ## Treatment2Predator_Control -2.53446    1.36576 22.66705  -1.856 0.076543 .  
    ## Treatment2Predator_Dull    -1.56940    1.37326 22.26054  -1.143 0.265255    
    ## Treatment2Water            -0.64936    0.95434 20.90728  -0.680 0.503700    
    ## temp                        1.21669    0.57544 21.26182   2.114 0.046464 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Br_S_H ylaydt Tr2C_D Tr2P_C Tr2P_D Trtm2W
    ## Brd_Sz_Htch -0.159                                          
    ## ylaydate    -0.744 -0.156                                   
    ## Trtmnt2Cn_D  0.171  0.149 -0.231                            
    ## Trtmnt2Pr_C  0.139  0.268 -0.009  0.463                     
    ## Trtmnt2Pr_D  0.074  0.451 -0.172  0.485  0.512              
    ## Tretmnt2Wtr  0.236  0.082 -0.365  0.630  0.492  0.541       
    ## temp         0.132  0.225 -0.748  0.098 -0.235  0.059  0.241

##### BKA model

``` r
nestling_m2 <- lmer(Nestling_Mass ~ Bacteria_Killing_Assay + 
                     Brood_Size_Hatching + temp +
                      ylaydate + Treatment2 + (1|Unit_Box), 
                   data = filter(nd, 
                                 Nestling_Fate == "Fledged",
                                 Bacteria_Killing_Assay != "NA") )
summary(nestling_m2)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Nestling_Mass ~ Bacteria_Killing_Assay + Brood_Size_Hatching +  
    ##     temp + ylaydate + Treatment2 + (1 | Unit_Box)
    ##    Data: filter(nd, Nestling_Fate == "Fledged", Bacteria_Killing_Assay !=  
    ##     "NA")
    ## 
    ## REML criterion at convergence: 560.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.90469 -0.49278  0.07571  0.48610  3.01971 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Unit_Box (Intercept) 2.275    1.508   
    ##  Residual             3.285    1.813   
    ## Number of obs: 132, groups:  Unit_Box, 34
    ## 
    ## Fixed effects:
    ##                            Estimate Std. Error      df t value Pr(>|t|)   
    ## (Intercept)                 17.3125    10.9901 22.1720   1.575  0.12935   
    ## Bacteria_Killing_Assay      -0.2499     0.9515 19.8695  -0.263  0.79556   
    ## Brood_Size_Hatching         -1.3362     0.3581 21.8515  -3.731  0.00117 **
    ## temp                         1.1968     0.5928 20.0308   2.019  0.05707 . 
    ## ylaydate                    -0.0876     0.1197 21.8447  -0.732  0.47189   
    ## Treatment2Control_Dull      -0.9558     1.2078 21.9654  -0.791  0.43718   
    ## Treatment2Predator_Control  -2.6066     1.4202 21.9134  -1.835  0.08006 . 
    ## Treatment2Predator_Dull     -1.6029     1.4053 21.1208  -1.141  0.26679   
    ## Treatment2Water             -0.7296     1.0179 19.8849  -0.717  0.48185   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Bc_K_A Br_S_H temp   ylaydt Tr2C_D Tr2P_C Tr2P_D
    ## Bctr_Klln_A  0.258                                                 
    ## Brd_Sz_Htch -0.179 -0.095                                          
    ## temp         0.160  0.133  0.209                                   
    ## ylaydate    -0.762 -0.281 -0.122 -0.749                            
    ## Trtmnt2Cn_D  0.197  0.136  0.134  0.114 -0.256                     
    ## Trtmnt2Pr_C  0.181  0.193  0.244 -0.203 -0.063  0.476              
    ## Trtmnt2Pr_D  0.090  0.074  0.442  0.068 -0.186  0.490  0.516       
    ## Tretmnt2Wtr  0.292  0.289  0.051  0.267 -0.416  0.637  0.519  0.539

##### Interaction model

``` r
nestling_m3 <- lmer(Nestling_Mass ~ Bacteria_Killing_Assay*ylaydate + 
                     Brood_Size_Hatching + temp +
                      Treatment2 + (1|Unit_Box), 
                   data = filter(nd, 
                                 Nestling_Fate == "Fledged",
                                 Bacteria_Killing_Assay != "NA") )
summary(nestling_m3)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: 
    ## Nestling_Mass ~ Bacteria_Killing_Assay * ylaydate + Brood_Size_Hatching +  
    ##     temp + Treatment2 + (1 | Unit_Box)
    ##    Data: filter(nd, Nestling_Fate == "Fledged", Bacteria_Killing_Assay !=  
    ##     "NA")
    ## 
    ## REML criterion at convergence: 559.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8637 -0.5641  0.0054  0.4731  2.9983 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Unit_Box (Intercept) 2.040    1.428   
    ##  Residual             3.318    1.822   
    ## Number of obs: 132, groups:  Unit_Box, 34
    ## 
    ## Fixed effects:
    ##                                 Estimate Std. Error      df t value Pr(>|t|)
    ## (Intercept)                       0.9654    15.9020 19.9640   0.061 0.952196
    ## Bacteria_Killing_Assay           42.6338    31.2011 19.1392   1.366 0.187647
    ## ylaydate                          0.0301     0.1434 20.3838   0.210 0.835812
    ## Brood_Size_Hatching              -1.4985     0.3651 19.9866  -4.104 0.000552
    ## temp                              1.2382     0.5717 17.2152   2.166 0.044635
    ## Treatment2Control_Dull           -0.6239     1.1901 18.5057  -0.524 0.606310
    ## Treatment2Predator_Control       -2.3377     1.3844 18.7620  -1.689 0.107854
    ## Treatment2Predator_Dull          -1.7192     1.3579 18.4830  -1.266 0.221209
    ## Treatment2Water                  -0.6210     0.9831 17.1603  -0.632 0.535919
    ## Bacteria_Killing_Assay:ylaydate  -0.3080     0.2239 19.1018  -1.375 0.184928
    ##                                    
    ## (Intercept)                        
    ## Bacteria_Killing_Assay             
    ## ylaydate                           
    ## Brood_Size_Hatching             ***
    ## temp                            *  
    ## Treatment2Control_Dull             
    ## Treatment2Predator_Control         
    ## Treatment2Predator_Dull            
    ## Treatment2Water                    
    ## Bacteria_Killing_Assay:ylaydate    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Bc_K_A ylaydt Br_S_H temp   Tr2C_D Tr2P_C Tr2P_D Trtm2W
    ## Bctr_Klln_A -0.740                                                        
    ## ylaydate    -0.852  0.587                                                 
    ## Brd_Sz_Htch  0.130 -0.326 -0.286                                          
    ## temp         0.066  0.059 -0.569  0.180                                   
    ## Trtmnt2Cn_D -0.021  0.207 -0.083  0.059  0.124                            
    ## Trtmnt2Pr_C  0.013  0.149  0.035  0.181 -0.193  0.490                     
    ## Trtmnt2Pr_D  0.111 -0.067 -0.189  0.438  0.063  0.464  0.497              
    ## Tretmnt2Wtr  0.136  0.087 -0.287  0.022  0.270  0.637  0.522  0.528       
    ## Bctr_Kll_A:  0.745 -1.000 -0.593  0.324 -0.056 -0.203 -0.143  0.069 -0.078

##### Model comparison

``` r
ICtab(nestling_m1, nestling_m2, nestling_m3, type = "AICc")
```

    ##             dAICc df
    ## nestling_m1  0.0  10
    ## nestling_m2  0.6  11
    ## nestling_m3  2.3  12

## Survival

### 2020 Return Rate

#### Return rate models

##### Null model

``` r
return_m1 <- glm(return2020 ~ success + ylaydate + Age + Treatment2,
              family = binomial(link = "logit"),
              data = df)
summary(return_m1)
```

    ## 
    ## Call:
    ## glm(formula = return2020 ~ success + ylaydate + Age + Treatment2, 
    ##     family = binomial(link = "logit"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.7360  -0.7646  -0.3535   0.6929   2.5041  
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                -0.27000   10.86830  -0.025   0.9802  
    ## successTRUE                 2.30432    0.94444   2.440   0.0147 *
    ## ylaydate                   -0.00490    0.07785  -0.063   0.9498  
    ## AgeSY                      -0.09686    0.70803  -0.137   0.8912  
    ## Treatment2Predator_Control -2.02818    1.16917  -1.735   0.0828 .
    ## Treatment2Water            -2.45812    1.06319  -2.312   0.0208 *
    ## Treatment2Control_Control  -1.70600    1.14393  -1.491   0.1359  
    ## Treatment2Control_Dull     -1.04942    1.20966  -0.868   0.3857  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 71.529  on 59  degrees of freedom
    ## Residual deviance: 57.260  on 52  degrees of freedom
    ## AIC: 73.26
    ## 
    ## Number of Fisher Scoring iterations: 5

##### BKA model

``` r
return_m2 <- glm(return2020 ~ success + ylaydate + Age + Treatment2 + Bacteria_Killing_Assay,
                 family = binomial(link = "logit"),
                 data = df)
summary(return_m2)
```

    ## 
    ## Call:
    ## glm(formula = return2020 ~ success + ylaydate + Age + Treatment2 + 
    ##     Bacteria_Killing_Assay, family = binomial(link = "logit"), 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5959  -0.7622  -0.3458   0.6404   2.6438  
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                -1.32885   10.88434  -0.122   0.9028  
    ## successTRUE                 2.34739    0.95850   2.449   0.0143 *
    ## ylaydate                    0.00547    0.07820   0.070   0.9442  
    ## AgeSY                      -0.04504    0.71261  -0.063   0.9496  
    ## Treatment2Predator_Control -2.20528    1.21557  -1.814   0.0696 .
    ## Treatment2Water            -2.58708    1.08950  -2.375   0.0176 *
    ## Treatment2Control_Control  -1.70156    1.14104  -1.491   0.1359  
    ## Treatment2Control_Dull     -1.17745    1.21207  -0.971   0.3313  
    ## Bacteria_Killing_Assay     -0.79191    0.97686  -0.811   0.4176  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 71.529  on 59  degrees of freedom
    ## Residual deviance: 56.596  on 51  degrees of freedom
    ## AIC: 74.596
    ## 
    ## Number of Fisher Scoring iterations: 5

##### Interaction model

``` r
return_m3 <- glm(return2020 ~ success + Age + Treatment2 + Bacteria_Killing_Assay*ylaydate,
                 family = binomial(link = "logit"),
                 data = df)
summary(return_m3)
```

    ## 
    ## Call:
    ## glm(formula = return2020 ~ success + Age + Treatment2 + Bacteria_Killing_Assay * 
    ##     ylaydate, family = binomial(link = "logit"), data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5669  -0.7417  -0.3440   0.5552   2.6309  
    ## 
    ## Coefficients:
    ##                                  Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                      12.36830   17.50447   0.707   0.4798  
    ## successTRUE                       2.25325    0.96608   2.332   0.0197 *
    ## AgeSY                            -0.11367    0.71732  -0.158   0.8741  
    ## Treatment2Predator_Control       -2.55984    1.29036  -1.984   0.0473 *
    ## Treatment2Water                  -2.64115    1.10780  -2.384   0.0171 *
    ## Treatment2Control_Control        -1.69134    1.14468  -1.478   0.1395  
    ## Treatment2Control_Dull           -1.36529    1.27866  -1.068   0.2856  
    ## Bacteria_Killing_Assay          -36.03552   32.29362  -1.116   0.2645  
    ## ylaydate                         -0.09158    0.12522  -0.731   0.4645  
    ## Bacteria_Killing_Assay:ylaydate   0.25277    0.23150   1.092   0.2749  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 71.529  on 59  degrees of freedom
    ## Residual deviance: 55.347  on 50  degrees of freedom
    ## AIC: 75.347
    ## 
    ## Number of Fisher Scoring iterations: 5

##### Model comparison

``` r
ICtab(return_m1, return_m2, return_m3, type = "AICc")
```

    ##           dAICc df
    ## return_m1  0.0  8 
    ## return_m2  2.1  9 
    ## return_m3  3.8  10

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-vitousekLingeringImpactStress2018" class="csl-entry">

Vitousek, Maren N., Conor C. Taff, Daniel R. Ardia, Jocelyn M. Stedman,
Cedric Zimmer, Timothy C. Salzman, and David W. Winkler. 2018. “The
Lingering Impact of Stress: Brief Acute Glucocorticoid Exposure Has
Sustained, Dose-Dependent Effects on Reproduction.” *Proceedings of the
Royal Society B: Biological Sciences* 285 (1882): 20180722–22.
<https://doi.org/10.1098/rspb.2018.0722>.

</div>

</div>
