# CALCULATE eFAST INDICES ####
# ---------------------------------------------------------------------------- #

#All outputs
sim.results.eFAST <- as.matrix(y)

for (i in (1:length(sim.results.eFAST[,1])))
{
  tell(Experience, sim.results.eFAST[i,])
}

plot(Experience)

Variance <- Experience$V # Total variance
names(Variance) <- GSA.parms

Done <- Experience$D1
names(Done) <- GSA.parms

Dt <- Experience$Dt
names(Dt) <- GSA.parms

first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
names(first_order) <- GSA.parms

total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
names(total_order) <- GSA.parms

lowry_data <- data.frame(
  Parameter = GSA.parms,
  Main.Effect = first_order,
  Interaction = total_order) %>% 
  mutate(
    Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
    Interaction = Interaction/sum(Interaction)  # to normalise
  ) 

write.csv(lowry_data, file = here("eFASTresults.csv"), row.names = FALSE)

# CP
CP.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CP") &
                                            !str_detect(rownames(sim.results.eFAST), "TT"), ]

for (i in (1:length(CP.sim.results.eFAST[,1])))
{
  tell(Experience, CP.sim.results.eFAST[i,])
}

plot(Experience)

Variance <- Experience$V # Total variance
names(Variance) <- GSA.parms

Done <- Experience$D1
names(Done) <- GSA.parms

Dt <- Experience$Dt
names(Dt) <- GSA.parms

first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
names(first_order) <- GSA.parms

total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
names(total_order) <- GSA.parms

lowry_data <- data.frame(
  Parameter = GSA.parms,
  Main.Effect = first_order,
  Interaction = total_order) %>% 
  mutate(
    Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
    Interaction = Interaction/sum(Interaction)  # to normalise
  ) 

write.csv(lowry_data, file = here("CP.eFASTresults.csv"), row.names = FALSE)

#CPTT
CPTT.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CPTT"), ]

for (i in (1:length(CPTT.sim.results.eFAST[,1])))
{
  tell(Experience, CPTT.sim.results.eFAST[i,])
}

plot(Experience)

Variance <- Experience$V # Total variance
names(Variance) <- GSA.parms

Done <- Experience$D1
names(Done) <- GSA.parms

Dt <- Experience$Dt
names(Dt) <- GSA.parms

first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
names(first_order) <- GSA.parms

total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
names(total_order) <- GSA.parms

lowry_data <- data.frame(
  Parameter = GSA.parms,
  Main.Effect = first_order,
  Interaction = total_order) %>% 
  mutate(
    Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
    Interaction = Interaction/sum(Interaction)  # to normalise
  ) 

write.csv(lowry_data, file = here("CPTT.eFASTresults.csv"), row.names = FALSE)

#CPTL
CPTL.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CPTL"), ]

for (i in (1:length(CPTL.sim.results.eFAST[,1])))
{
  tell(Experience, CPTL.sim.results.eFAST[i,])
}

plot(Experience)

Variance <- Experience$V # Total variance
names(Variance) <- GSA.parms

Done <- Experience$D1
names(Done) <- GSA.parms

Dt <- Experience$Dt
names(Dt) <- GSA.parms

first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
names(first_order) <- GSA.parms

total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
names(total_order) <- GSA.parms

lowry_data <- data.frame(
  Parameter = GSA.parms,
  Main.Effect = first_order,
  Interaction = total_order) %>% 
  mutate(
    Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
    Interaction = Interaction/sum(Interaction)  # to normalise
  ) 

write.csv(lowry_data, file = here("CPTL.eFASTresults.csv"), row.names = FALSE)

#CL_ic
CL_ic.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CL_ic"), ]

for (i in (1:length(CL_ic.sim.results.eFAST[,1])))
{
  tell(Experience, CL_ic.sim.results.eFAST[i,])
}

plot(Experience)

Variance <- Experience$V # Total variance
names(Variance) <- GSA.parms

Done <- Experience$D1
names(Done) <- GSA.parms

Dt <- Experience$Dt
names(Dt) <- GSA.parms

first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
names(first_order) <- GSA.parms

total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
names(total_order) <- GSA.parms

lowry_data <- data.frame(
  Parameter = GSA.parms,
  Main.Effect = first_order,
  Interaction = total_order) %>% 
  mutate(
    Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
    Interaction = Interaction/sum(Interaction)  # to normalise
  ) 

write.csv(lowry_data, file = here("CL_ic.eFASTresults.csv"), row.names = FALSE)





