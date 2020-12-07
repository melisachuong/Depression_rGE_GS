
sym_1= c()
sym_2= c()
sym_3= c()
sym_4= c()
sym_5= c()
sym_6= c()
sym_7= c()
sym_8= c()
sym_9= c()

for(j in seq(0,1, by=0.1)){ #j being V_noise; the PRS noise
  for(k in seq(0.1,1, by=0.1)){ #k being heritability
    for(i in 1:20){ #i being the iteration/no. of replications
      for(l in seq(0.1,1, by=0.1)){ #l being the tagged genetic variance
        
        k=0.8
        l=0.45
        j=0.5
        i = 1
      iteration <- i
      
      h2 <- k
      V_P <- 1
      V_A <- h2*V_P
      V_Tagged <- l*V_A #the value of l should be SNP heritability i.e. the proportion of the genetic variance that is actually tagged
      V_nTagged <- V_A - V_Tagged
      V_E <- (1-h2)*V_P
      n <- 10000



#Creating a parent genetic variable that captures ALL tagged genetic variance
maternal_genetic_tagged <- rnorm(n, mean=0, sqrt(V_Tagged))
paternal_genetic_tagged <- rnorm(n, mean=0, sqrt(V_Tagged))


#Creating a parent genetic variable that captures ALL non-tagged genetic variance
maternal_genetic_ntagged <- rnorm(n, mean=0, sqrt(V_nTagged))
paternal_genetic_ntagged <- rnorm(n, mean=0, sqrt(V_nTagged))


#Creating a parent variable that accounts for all genetic variance
maternal_genetic <- maternal_genetic_tagged + maternal_genetic_ntagged
paternal_genetic <- paternal_genetic_tagged + paternal_genetic_ntagged


#Creating a parent variable that accounts for all environmental variance
maternal_environment <- rnorm(n, mean=0, sd = sqrt(V_E))
paternal_environment <- rnorm(n, mean=0, sd = sqrt(V_E))


#Creating a parent phenotypic variable that is partially genetic and partially environmental
maternal_phenotype <- maternal_genetic + maternal_environment
paternal_phenotype <- paternal_genetic + paternal_environment


#Simulating offspring segregation for tagged and non-tagged genetic variance from each parent
offspring_segregation_tagged_mat <- rnorm(n, 0, sd=sqrt(V_Tagged/4)) #the /4 value comes from quantitative genetics theory
offspring_segregation_ntagged_mat <- rnorm(n, 0, sd=sqrt(V_nTagged/4))
offspring_segregation_tagged_pat <- rnorm(n, 0, sd=sqrt(V_Tagged/4))
offspring_segregation_ntagged_pat <- rnorm(n, 0, sd=sqrt(V_nTagged/4))


#Creating offspring genetic, environment and phenotypic variables
offspring_genetic_tagged_mat <- maternal_genetic_tagged/2 + offspring_segregation_tagged_mat
offspring_genetic_ntagged_mat <- maternal_genetic_ntagged/2 + offspring_segregation_ntagged_mat
offspring_genetic_tagged_pat <- paternal_genetic_tagged/2 + offspring_segregation_tagged_pat
offspring_genetic_ntagged_pat <- paternal_genetic_ntagged/2 + offspring_segregation_ntagged_pat

offspring_genetic <- offspring_genetic_tagged_mat + offspring_genetic_ntagged_mat + offspring_genetic_tagged_pat + offspring_genetic_ntagged_pat

offspring_environment <- rnorm(n, mean=0, sd = sqrt(V_E))

offspring_phenotype <- offspring_genetic + offspring_environment


#Create parent transmitted and non-transmitted genetic variables 
maternal_genetic_tagged_ntrans <- maternal_genetic_tagged - offspring_genetic_tagged_mat
maternal_genetic_tagged_trans <- offspring_genetic_tagged_mat
paternal_genetic_tagged_ntrans <- paternal_genetic_tagged - offspring_genetic_tagged_pat
paternal_genetic_tagged_trans <- offspring_genetic_tagged_pat

maternal_genetic_ntagged_ntrans <- maternal_genetic_ntagged - offspring_genetic_ntagged_mat
maternal_genetic_ntagged_trans <- offspring_genetic_ntagged_mat
paternal_genetic_ntagged_ntrans <- paternal_genetic_ntagged - offspring_genetic_ntagged_pat
paternal_genetic_ntagged_trans <- offspring_genetic_ntagged_pat

maternal_genetic_trans  <- maternal_genetic_ntagged_trans + maternal_genetic_tagged_trans
maternal_genetic_ntrans <- maternal_genetic_tagged_ntrans + maternal_genetic_ntagged_ntrans
paternal_genetic_trans  <- paternal_genetic_ntagged_trans + paternal_genetic_tagged_trans
paternal_genetic_ntrans <- paternal_genetic_tagged_ntrans + paternal_genetic_ntagged_ntrans


#Simulating parental PRSs where PRSs are derived of parental genetic variables + noise derived using the noise variance variable created in the previous step

#First we need to create noise terms to ensure we can play around with the PRS accuracy
V_noise = j #j represents the noisyness of the PRS, as j becomes larger, the PRS becomes more noisy

maternal_tagged_trans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
maternal_ntagged_trans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
maternal_tagged_ntrans_noise<- rnorm(n, 0, sd=sqrt(V_noise/4))
maternal_ntagged_ntrans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
maternal_noise <- rnorm(n, 0, sd=sqrt(V_noise))

paternal_tagged_trans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
paternal_ntagged_trans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
paternal_tagged_ntrans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
paternal_ntagged_ntrans_noise <- rnorm(n, 0, sd=sqrt(V_noise/4))
paternal_noise <- rnorm(n, 0, sd=sqrt(V_noise))

#We will separately create parental polygenic risk scores 
maternal_trans_prs <- maternal_genetic_tagged_trans + maternal_genetic_ntagged_trans + maternal_tagged_trans_noise + maternal_ntagged_trans_noise
paternal_trans_prs <- paternal_genetic_tagged_trans + paternal_genetic_ntagged_trans + paternal_tagged_trans_noise + paternal_ntagged_trans_noise

maternal_ntrans_prs <- maternal_genetic_tagged_ntrans + maternal_genetic_ntagged_ntrans + maternal_tagged_ntrans_noise + maternal_ntagged_ntrans_noise
paternal_ntrans_prs <- paternal_genetic_tagged_ntrans + paternal_genetic_ntagged_ntrans + paternal_tagged_ntrans_noise + paternal_ntagged_ntrans_noise

maternal_prs <- maternal_genetic + maternal_noise
paternal_prs <- paternal_genetic + paternal_noise

offspring_prs <- offspring_genetic + rnorm(n, 0, sd=sqrt(V_noise/2)) + (maternal_tagged_trans_noise + maternal_ntagged_trans_noise + 
                 maternal_tagged_ntrans_noise + maternal_ntagged_ntrans_noise + paternal_tagged_trans_noise + paternal_ntagged_trans_noise +
                 paternal_tagged_ntrans_noise + paternal_ntagged_ntrans_noise)/2


#Creating a data.frame with all important variables included
cohort <- data.frame(maternal_genetic, maternal_environment, maternal_phenotype, maternal_trans_prs, maternal_ntrans_prs,
                     paternal_genetic, paternal_environment, paternal_phenotype, paternal_trans_prs, paternal_ntrans_prs,
                     offspring_genetic, offspring_environment, offspring_phenotype, offspring_prs)


regression <- summary(lm(offspring_phenotype ~ offspring_prs + maternal_trans_prs + paternal_trans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_1 <- rbind(sym_1,reginfo)
write.table(sym_1, "sim_outcome_offprs_parenttrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + maternal_ntrans_prs + paternal_ntrans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_2 <- rbind(sym_2,reginfo)
write.table(sym_2, "sim_outcome_offprs_parentnontrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + maternal_trans_prs + paternal_trans_prs + maternal_ntrans_prs + paternal_ntrans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_3 <- rbind(sym_3,reginfo)
write.table(sym_3, "sim_outcome_offprs_parenttrans_ntrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + maternal_trans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_4 <- rbind(sym_4,reginfo)
write.table(sym_4, "sim_outcome_offprs_mattrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + maternal_ntrans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_5 <- rbind(sym_5,reginfo)
write.table(sym_5, "sim_outcome_offprs_matntrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + maternal_trans_prs + maternal_ntrans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_6 <- rbind(sym_6,reginfo)
write.table(sym_6, "sim_outcome_offprs_matrans_ntrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + paternal_trans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_7 <- rbind(sym_7,reginfo)
write.table(sym_7, "sim_outcome_offprs_pattrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + paternal_ntrans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_8 <- rbind(sym_8,reginfo)
write.table(sym_8, "sim_outcome_offprs_patntrans", row.names=T, quote=F)

regression <- summary(lm(offspring_phenotype ~ offspring_prs + paternal_trans_prs + paternal_ntrans_prs, data=cohort))
reginfo <-  cbind(regression$coefficients, regression$adj.r.squared, V_noise, iteration, k,l)  
sym_9 <- rbind(sym_9,reginfo)
write.table(sym_9, "sim_outcome_offprs_pattrans_ntrans", row.names=T, quote=F)}}}}



















# Checks the variances and correlations
var(offspring_genetic)
var(offspring_phenotype)
var(offspring_genetic)/var(offspring_phenotype)
var(maternal_genetic)/var(maternal_phenotype)
var(paternal_genetic)/var(paternal_phenotype)

var(maternal_genetic_trans)
var(maternal_genetic_ntrans)
var(paternal_genetic_trans)
var(paternal_genetic_ntrans)

summary(lm(offspring_genetic~maternal_genetic))
summary(lm(offspring_genetic~paternal_genetic))
summary(lm(offspring_genetic~maternal_genetic_ntrans))
summary(lm(offspring_genetic~paternal_genetic_ntrans))
summary(lm(offspring_genetic~maternal_genetic_trans))
summary(lm(offspring_genetic~paternal_genetic_trans))

cor(maternal_prs, offspring_prs)
cor(paternal_prs, offspring_prs)

cor((maternal_trans_prs+maternal_ntrans_prs), offspring_prs)
cor((paternal_trans_prs+paternal_ntrans_prs), offspring_prs) #we expect these results to be similar to the previous 2 correlations

cor((maternal_prs-maternal_trans_prs), maternal_ntrans_prs)
cor((paternal_prs-paternal_trans_prs), paternal_ntrans_prs)
cor((maternal_trans_prs+paternal_trans_prs), offspring_prs)



