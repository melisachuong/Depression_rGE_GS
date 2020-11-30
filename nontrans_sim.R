
```{r}
sym_1= c()
sym_2= c()
sym_3= c()
sym_4= c()
sym_5= c()
sym_6= c()
sym_7= c()

for(j in seq(0,1, by=0.1)){ #j being V_noise; the PRS noise
  for(k in seq(0.1,1, by=0.1)){ #k being heritability
    for(i in 1:20){ #i being the iteration/no. of replications
        
        iteration <- i
      
        
        h2 = k
        V_P = 1
        V_A = h2*V_P
        V_Trans = 0.5*V_A #You will ALWAYS transmit 50% of the genetic variance
        V_nTrans = V_A - V_Trans
        V_E = (1-h2)*V_P
        n = 10000
        
        #Creating a genetic variable that captures ALL transmitted genetic variance
        maternal_genetic_trans = rnorm(n, mean=0, sqrt(V_Trans))
        paternal_genetic_trans = rnorm(n, mean=0, sqrt(V_Trans))
        
        
        #Creating a genetic variable that captures ALL non-transmitted genetic variance
        maternal_genetic_ntrans = rnorm(n, mean=0, sqrt(V_nTrans))
        paternal_genetic_ntrans = rnorm(n, mean=0, sqrt(V_nTrans))
        
        
        #Creating a variable that accounts for all genetic variance, sum of transmitted and non-transmitted
        maternal_genetic = maternal_genetic_trans + maternal_genetic_ntrans
        paternal_genetic = paternal_genetic_trans + paternal_genetic_ntrans
        
        
        #Creating a genetic variable that accounts for all environmental variance
        maternal_environment = rnorm(n, mean=0, sd = sqrt(V_E))
        paternal_environment = rnorm(n, mean=0, sd = sqrt(V_E))
        
        #Creating a variable that accounts for all environmental variance
        maternal_phenotype = maternal_genetic + maternal_environment
        paternal_phenotype = paternal_genetic + paternal_environment
        
        
        #Simulating offspring segregation
        offspring_segregation = rnorm(n, 0, sd=sqrt(V_A/2))
        
        #Creating an offspring genetic variable
        offspring_genetic =  (maternal_genetic_trans +  paternal_genetic_trans) + offspring_segregation
        
        offspring_environment = rnorm(n, mean=0, sd = sqrt(V_E))
        
        offspring_phenotype = offspring_genetic + offspring_environment
        
        #Specifying a noise variance variable that is between 0-V_P
        V_noise = j
        
        
        #Simulating parental PRSs where PRSs are derived of parental genetic variables + noise derived using the noise variance variable created in the previous step
        
        maternal_tnoise = rnorm(n, 0, sd=sqrt(V_noise/2))
        maternal_ntnoise = rnorm(n, 0, sd=sqrt(V_noise/2))
        
        maternal_tprs = maternal_genetic_trans + maternal_tnoise
        maternal_ntprs = maternal_genetic_ntrans + maternal_ntnoise
        maternal_prs = maternal_tprs + maternal_ntprs
        
        
        paternal_tnoise = rnorm(n, 0, sd=sqrt(V_noise/2))
        paternal_ntnoise = rnorm(n, 0, sd=sqrt(V_noise/2))
        
        paternal_tprs = paternal_genetic_trans + paternal_tnoise
        paternal_ntprs = paternal_genetic_ntrans + paternal_ntnoise
        paternal_prs = paternal_tprs + paternal_ntprs
        
        #Simulating offspring PRSs derived from offspring genetic variable + noise derived from the average parental noise and offspring noise
        #offspring_prs = offspring_genetic + rnorm(n, 0, sd=sqrt(V_noise/2)) + (maternal_tnoise + paternal_tnoise)/2
        offspring_prs = maternal_tprs + paternal_tprs
        
        #QUICK CHECKS
        cor(maternal_prs, offspring_prs)
        cor(paternal_prs, offspring_prs)
        
        cor((maternal_tprs+paternal_tprs),offspring_prs)
        
        
        #Creating a data.frame with all important variables included
        cohort <- data.frame(maternal_genetic, maternal_environment, maternal_phenotype, maternal_prs, maternal_tprs, maternal_ntprs,
                             paternal_genetic, paternal_environment, paternal_phenotype, paternal_prs, paternal_tprs, paternal_ntprs,
                             offspring_genetic, offspring_environment, offspring_phenotype, offspring_prs)}}}
        
        

```    
