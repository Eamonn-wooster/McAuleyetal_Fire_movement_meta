## Testing Decline Effect on Early vs Late Studies
# Workflow from ChatGPT

#prepare and use es from script 06

es$Publication_Year <- as.numeric(as.character(data$Publication_Year))

median_year <- median(es$Publication_Year)

early <- es[es$Publication_Year <= median_year, ]
late  <- es[es$Publication_Yea >  median_year, ]

early.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = early) 

late.vcv <- vcalc(vi, cluster = Study_ID, obs = Obs_ID, rho = 0.5, # rho is usually 0.5 or 0.8
                   data = late) 


res_early1 <- rma.mv(yi = yi, V = early.vcv,
                      random = list(~1 | Study_ID,
                                    ~1 | Species_tree, # phylo effect 
                                    ~1 | Species2, # non-phylo effect 
                                    ~1 | Obs_ID), 
                      data =  early,
                      control = list(optimizer="BFGS"),
                      test = "t",
                      sparse = TRUE,
                      R = list(Species_tree = cor1)
)

res_late1 <- rma.mv(yi = yi, V = late.vcv,
                     random = list(~1 | Study_ID,
                                   ~1 | Species_tree, # phylo effect 
                                   ~1 | Species2, # non-phylo effect 
                                   ~1 | Obs_ID), 
                     data =  late,
                     control = list(optimizer="BFGS"),
                     test = "t",
                     sparse = TRUE,
                     R = list(Species_tree = cor1)
)

summary(res_early1)
summary(res_late1)
