# Group numbers/ranges in strings into matched categories
dat0[,PR.fac:=case_when(
  grepl("^<1|^0$", PR.stain) ~ "<1%", 
  grepl("^>1$|^1-5$|^2$|^<10$|^<5$", PR.stain) ~ "1-10%", 
  grepl("^>10$|^15$|^45$|^50$", PR.stain) ~ "11-50%",  
  grepl("^75$|^95$|^>50$", PR.stain) ~ "51-100%",
  grepl("^unk", PR.stain) ~ NA,
  TRUE ~ PR.stain)]
