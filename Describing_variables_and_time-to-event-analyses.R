### Codes for "Construction of a machine learning model predicting febrile neutropenia in patients with AML undergoing chemotherapy" ###

# ---------
# Figure 1A
# ---------

g <- ggplot(df, 
              aes(
                x = days, 
                y = b_neut,
                color = cycle                                                                        
                 )) +
  geom_hline(yintercept = 0.5, 
             linetype = "dashed", 
             color = "black", 
             linewidth = 0.5, 
             alpha = 0.35) +
  geom_smooth(method = "loess",
              se=TRUE,          
              alpha = 0.25
  ) +
  scale_color_manual(values = c("IND" = "#E41A1C", "CONS" = "#377EB8")) +
  guides(color = guide_legend(title = "Mean [95% CI]"))  + 
  labs(x = "Days", 
       y = "ANC 10⁹/L") + 
  scale_y_continuous(breaks = seq(0,4, by = 0.5)) +
  scale_x_continuous(breaks = seq(1,28, by = 3)) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# ---------
# Figure 1B
# ---------

g <- ggplot(df, 
                aes(
                  x = days, 
                  y = b_neut_trend_from_3_d,
                  color = cycle  
                )) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black", 
             linewidth = 0.5, 
             alpha = 0.35) +
  geom_smooth(method = "loess", 
              se=TRUE, 
              alpha = 0.25) +
  scale_color_manual(values = c("IND" = "#E41A1C", "CONS" = "#377EB8")) +
  guides(color = guide_legend(title = "Mean [95% CI]"))  +  
  labs(x = "Days", 
       y = "ANC 10⁹/L Δ3 days") +
  scale_y_continuous(breaks = seq(-3, 3, by = 0.5)
  ) +
  scale_x_continuous(breaks = seq(1,28, by = 3)) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# ---------
# Figure 1C
# ---------

df_final <- df %>%                                                                             
  dplyr::select(id, cycle, cycle_nro, !!df1, !!df_3_min, !!df_3_max) %>%                                           
  reshape2::melt(id.vars =  c("id", "cycle", "cycle_nro"), variable.name = "labtest", value.name = "result") %>%    
  dplyr::group_by(labtest, cycle, cycle_nro) %>%                                                                               
  summarise(result = median(result, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(labtest) %>%                                                                                               
  dplyr::mutate(result_norm_scale = scales::rescale(result, to = c(0,1))) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(cycle_name = ifelse(
    sykli == "2", "IND", "CONS"),                                                                            
    combo_name = paste0(cycle_name, cycle_nro),
    combo_name = factor(
      combo_name, 
      levels = c("IND1", "IND2", "CONS1", "CONS2", "CONS3")
    ), 
    subgroup_name = factor(                                                                                     
      case_when(
        labtest %in% c(cbc, cbc_3_min, cbc_3_max) ~ "CBC",
        labtest %in% c(tls, tls_3_min, tls_3_max) ~ "TLS",             
        labtest == "p_crp" | labtest == "p_crp_trend_from_3_d_min" | labtest == "p_crp_trend_from_3_d_max" ~ "Inflammatory", 
        labtest %in% c(liver, liver_3_min, liver_3_max) ~ "Liver function",     
        TRUE ~ labtest),
      levels = c("Liver function", "TLS", "Inflammatory", "CBC")
    ),
    subgroup = as.integer(subgroup_name),                                                                       
    type_name = factor(                                                                                         
      case_when(
        labtest %in% !!df1 ~ "Day 1",
        labtest %in% !!df_3_min ~ "Steepest 3-day\n↓ trend",
        labtest %in% !!df_3_max ~ "Steepest 3-day\n↑ trend",
        TRUE ~ labtest),
      levels = c("Day 1", "Steepest 3-day\n↓ trend", "Steepest 3-day\n↑ trend")
    ),
    type = as.integer(type_name))

g <- ggplot(
  df_final, 
  aes(
    x = combo_name,                                           
    y = reorder(labtest, subgroup),              
    fill = result_norm_scale            
    )
  ) +
  facet_wrap(vars(factor(type, 
                         levels = c(1,2,3), 
                         labels = c("Day 1", "Steepest 3-day\n↓ trend", "Steepest 3-day\n↑ trend"))), nrow =  1) + 
  geom_tile(
    colour = "white",                    
    size = 0.5                           
  ) +
  labs(
    x = "Cycles",                         
    y = "Laboratory values",              
    title = ""                            
  ) +
  scale_fill_distiller(
    palette = "RdBu", 
    name = "Normalized\nvalues"           
  ) +
  scale_x_discrete(                                                                   
    breaks = unique(df$combo_name),                              
    labels = unique(df$combo_name),                           
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = "right",            
    legend.direction = "vertical",        
    legend.title = element_text(          
      colour = "black", 
      hjust = 0.5,
      size = 6,
      face = "bold"
    ), 
    legend.margin = margin(grid::unit(0, "cm")), 
    legend.text = element_text(          
      colour = "black", 
      size = 6
    ),
    legend.key.height = grid::unit(0.5, "cm"), 
    legend.key.width = grid::unit(0.5, "cm"),  
    plot.title.position = "panel",        
    plot.title = element_text(            
      colour = "black", 
      hjust = 0.4, 
      vjust = -1, 
      size = 10, 
      face = "bold"
    ),
    strip.text = element_text(size = 5), 
    strip.background = element_rect(fill = "grey90", colour = "black"),
    axis.title.y = element_text(
      size = 7,
      face = "bold",
      margin = margin(r = 8)
    ),
    axis.title.x = element_text(
      size = 7,
      face = "bold",
      margin = margin(t = 8)
    ),
    axis.text.x = element_text(           
      size = 5,                                                                        
      colour = "black",
      angle = 45,
      hjust = 1,
      margin = margin(t = 1)
    ),
    axis.text.y = element_text(           
      size = 5, 
      colour = "black"
    ),
    axis.ticks = element_line(linewidth = 0.5),
    plot.margin = margin(                 
      0.7, 0.4, 0.1, 0.2, "cm"
    ),
    panel.border = element_blank(),       
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    panel.background = element_rect(color = "black", fill = "white", linewidth = 1),
    legend.box.margin = margin(-5, -5, -5, -5), 
  )

# ---------
# Figure 1D
# ---------

cox_fun <- sapply(covariates,
                  function(x) as.formula(paste('Surv(df$follow_up_period, df$infection)~', x, sep="")))

model_list <- lapply(cox_fun, function(x){coxph(x, data = df)})

univ_results <- lapply(model_list,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$coefficients[5], digits=5)
                         beta<-signif(x$coef[1], digits=5);
                         HR <-signif(x$coef[2], digits=5);
                         HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 5)
                         HR.confint.upper <- signif(x$conf.int[1,"upper .95"], 5)
                         HR1 <- paste0(HR, " (",
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                         C <- signif(x$concordance[1], digits = 5)
                         res<-c(beta, HR, HR1, p.value, C)
                         names(res)<-c("beta", "HR", "HR (95% CI for HR)", "p.value", "C")
                         return(res)
                       })

univ_results2 <- t(as.data.frame(univ_results)) %>%
  as.data.frame %>%
  tibble::rownames_to_column("names") %>%
  mutate(
    beta = as.numeric(as.character(beta)),
    HR = as.numeric(as.character(HR)),
    p.value = as.numeric(as.character(p.value)),
    p.value_adj = p.adjust(p.value, method = "BH", n = nrow(.)),
    C = as.numeric(as.character(C))) 

df1 <- df %>% 
  dplyr::select(names, HR, p.value, p.value_adj, category, C) %>%                                       
  dplyr::mutate(                    
    HR_log10 = log10(HR),                                                                                                     
    p.value_adj_miinuslog10 = -log10(p.value_adj),                                                      
    p.value_scaled = ifelse(HR >= 1, p.value_adj_miinuslog10, -p.value_adj_miinuslog10),                
    order1 = 1:nrow(.),                                                                                
    names = gsub("_", " ", names),                                                                      
    sig = ifelse(p.value_adj >= 0.05, "ns",                                                             
                 ifelse(HR >= 1, "HR≥1", "HR<1")),
    category2 = factor(ifelse(sig != "ns", df$category, "ns"),                                          
                       levels = c("LABS", "GEN", "CYTO", "OTHER VARIABLES", "ns"), 
                       labels = c("LABS", "GEN", "CYTO", "OTHER VARIABLES", "ns")),
    point_size = abs(HR_log10)+3                                                                        
  ) 

df1 = df1 %>%
  dplyr::mutate(p.value_scaled_2 = pmin(p.value_scaled, 20, na.rm = TRUE),                             
                p.value_scaled_2 = pmax(p.value_scaled_2, -20, na.rm = TRUE),                           
                category3 = factor(ifelse(p.value_scaled_2 > 0 & p.value_adj<0.05, 2,                    
                                          ifelse(p.value_scaled_2 < 0 & p.value_adj<0.05, 1, 0))),      
                point_size1 = pmin(point_size, ceiling(sort(unique(point_size), decreasing = TRUE)[2]), na.rm = TRUE), 
                point_size1 = ifelse(p.value_adj > 0.05, 0.5, point_size1))

limits_df <- df1 %>%
  group_by(category) %>%                                
  summarize(ymax = max(abs(p.value_scaled))) %>%        
  mutate(ymin = -ymax) %>%                              
  ungroup()

df1 <- df1 %>%                                          
  left_join(limits_df)

g = ggplot(df1, aes(order1, p.value_scaled)) +                                                            
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +              
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +               
  geom_point(aes(fill = category3), size = df1$point_size1, shape = 21) +                                 
  scale_fill_manual(values = c("black", "#377EB8", "#E41A1C"),                                            
                    labels = c("ns", "HR<1", "HR≥1")) +                                                   
  guides(fill = guide_legend("Scaled log10 HR",                                                            
                             nrow = 1,
                             title.hjust = 0.5,
                             override.aes = list(size = 5)),
         size = "none") +
  labs(x = "Rank", y = "Scaled -log10 AdjP") +                                                            
  theme_bw() +                                                                                            
  theme(                                                                                                  
    plot.title = element_text(hjust = 0.5, size = 23, face = "bold", margin=margin(10,0,20,0)),               
    axis.title.y = element_text(face = "bold", size = 17, color = 'black', margin=margin(0,20,0,0)),      
    axis.title.x = element_text(hjust = 0.5, face = "bold", size = 17, color = 'black', margin=margin(20,0,0,0)), 
    axis.text.x = element_text(size = 15, colour = "black"),                                              
    strip.text = element_text(size = 15),                                                                 
    strip.background = element_rect(fill = "grey90", colour = "black"),                                   
    legend.title = element_text(size = 17, face = "bold", colour = "black"),                              
    legend.text = element_text(size = 17, colour = "black"),
    legend.direction = 'horizontal',                                                                      
    legend.key = element_rect(linewidth = 3),                                                             
    legend.key.size = unit(1.5, 'lines'),
    legend.position = "bottom",                                                                           
    panel.grid.major = element_blank(),                                    
    panel.grid.minor = element_blank(),   
  ) +
  coord_cartesian(clip = "off") +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  facet_wrap(vars(factor(category,                                                                        
                         levels = c("LABS", "CYTO", "GEN", "OTHER VARIABLES"), 
                         labels = c("LABORATORY\nVARIABLES", "CYTOMORPHOLOGICAL\nCHANGES", "GENETIC\nALTERATIONS", "OTHERS"))),
             scales = "free_y", 
             nrow =  1) +
  geom_label_repel(data = df1 %>%                                                                         
                     dplyr::group_by(category) %>%                                                        
                     dplyr::arrange(desc(df1$C)) %>%                                                      
                     dplyr::filter(sig != "ns") %>%                                                       
                     dplyr::slice(1:15) %>%                                                               
                     dplyr::ungroup(),
                   size = 4,                                                                              
                   aes(label = names_2),                                                                  
                   box.padding = 0.7,                                                                     
                   max.overlaps = Inf, 
                   force = 20, 
                   force_pull = 20) 

# ---------
# Figure 1E
# ---------

dfi = dfi %>% 
  dplyr::select(id, infection_binary2, countdown1, countdown2, follow_up_period, cycle_IND, time_to_first_neutropenia) %>%
  dplyr::mutate(
    risk_group = case_when(
      time_to_first_neutropenia == 1 ~ "1 day",
      time_to_first_neutropenia == 2 | time_to_first_neutropenia == 3 | time_to_first_neutropenia == 4 | time_to_first_neutropenia == 5 | time_to_first_neutropenia == 6 | time_to_first_neutropenia == 7 ~ "2–7 days",
      time_to_first_neutropenia >=7 ~ "≥8 days"
    ))

dfi$risk_group = factor(dfi$risk_group, levels = c("1 day", "2–7 days", "≥8 days"),  labels = c("1 d", "2–7 d", "≥8 d"))

fit1 = survfit(Surv(follow_up_period, infection_binary2) ~ risk_group, data = dfi)

pval1 <- surv_pvalue(fit1, data = dfi)$pval
p_label1 <- if (pval1 < 0.001) {
  "italic(p) < 0.001"
} else {
  sprintf("italic(p) == %s", formatC(pval1, format = "f", digits = 3))
}

g1 <- ggsurvplot(fit1,
                fun = "event",  
                data = dfi,
                palette = "Set1",   
                size = 1,  
                ggtheme = theme_classic() +
                  theme(legend.title = element_text(face = "bold"),   
                        legend.text  = element_text(face = "plain")),  
                font.main = c(15, "bold", "black"),  
                font.x = c(14, "bold", "black"), 
                font.y = c(14, "bold", "black"),  
                font.tickslab = c(12, "black"),  
                conf.int = FALSE,  
                pval = FALSE, 
                risk.table.fontsize = 4.4,  
                risk.table.pos = "out",  
                tables.y.text = FALSE,  
                tables.theme = theme_cleantable(),
                risk.table = c("absolute"), 
                risk.table.col = "strata", 
                risk.table.height = 0.25,  
                ylab = "Probability of FN",  
                xlab = "Days",  
                surv.scale = "percent",  
                break.x.by = 3,  
                ylim = c(0, 1),  
                xlim = c(0, 28),  
                censor = TRUE,  
                censor.shape = 108,  
                censor.size = 4,  
                font.legend= c(12, "black"),  
                legend.title = "Onset of neutropenia (IND)",                            
                legend.labs = levels(dfi$risk_group),  
                legend = "top")

g1$plot = g1$plot + 
  scale_x_continuous(
    breaks = seq(0, 28, by = 3),
    labels = function(x) x + 1 # as the first day of the follow-up period was labeled as "Day 1" instead of "Day 0"
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           size = 4,
           label = p_label1, parse = TRUE, hjust = 0) +
  guides(color = guide_legend(nrow = 1))

# ---------
# Figure 1F
# ---------

dfk$tertile = cut(
  dfk$time_to_first_neutropenia,
  breaks = quantile(dfk$time_to_first_neutropenia, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("low", "int", "high")
)

dfk = dfk %>% 
  dplyr::select(id, infection_binary2, countdown1, countdown2, follow_up_period, cycle_IND, time_to_first_neutropenia, tertile) %>%
  dplyr::mutate(
    risk_group = case_when(
      tertile == "low" ~ "≤11 days",
      tertile == "int" ~ "12–13 days",
      tertile == "high" ~ "≥14 days"))

dfk$risk_group = factor(dfk$risk_group, levels = c("≤11 days", "12–13 days", "≥14 days"),  labels = c("≤11 d", "12–13 d", "≥14 d"))

fit2 = survfit(Surv(follow_up_period, infection_binary2) ~ risk_group, data = dfk) 

pval2 <- surv_pvalue(fit2, data = dfk)$pval
p_label2 <- if (pval2 < 0.001) {
  "italic(p) < 0.001"
} else {
  sprintf("italic(p) == %s", formatC(pval2, format = "f", digits = 3))
}

g2 <- ggsurvplot(fit2,
                fun = "event",  
                data = dfk,
                palette = "Set1",  
                size = 1,   
                ggtheme = theme_classic() +
                  theme(legend.title = element_text(face = "bold"),   
                        legend.text  = element_text(face = "plain")),  
                font.main = c(15, "bold", "black"),  
                font.x = c(14, "bold", "black"),  
                font.y = c(14, "bold", "black"),  
                font.tickslab = c(12, "black"),  
                conf.int = FALSE,  
                pval = FALSE,  
                risk.table.fontsize = 4.4,   
                risk.table.pos = "out",  
                tables.y.text = FALSE,  
                tables.theme = theme_cleantable(),  
                risk.table = c("absolute"), 
                risk.table.col = "strata",  
                risk.table.height = 0.25,  
                ylab = "Probability of FN",  
                xlab = "Days", 
                surv.scale = "percent",  
                break.x.by = 3,  
                ylim = c(0, 1), 
                xlim = c(0, 28),  
                censor = TRUE,  
                censor.shape = 108,  
                censor.size = 4,  
                font.legend= c(12, "black"),  
                legend.title = "Onset of neutropenia (CONS)", 
                legend.labs = levels(dfk$risk_group),  
                legend = "top")

g2$plot = g2$plot + 
  scale_x_continuous(
    breaks = seq(0, 28, by = 3),
    labels = function(x) x + 1
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           size = 4,
           label = p_label2, parse = TRUE, hjust = 0) +
  guides(color = guide_legend(nrow = 1))

# ---------
# Figure 1G
# ---------

dfi$tertile = cut(
  dfi$neutropenia_days_050,
  breaks = quantile(dfi$neutropenia_days_050, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("low", "int", "high")
)

dfi = dfi %>% 
  dplyr::select(id, infection, follow_up_period, cycle, neutropenia_days_050, tertile) %>%
  dplyr::mutate( 
    risk_group2 = case_when(
      tertile == "low" ~ "T1:\n≤19 days",
      tertile == "int" ~ "T2:\n20-23 days",
      tertile == "high" ~ "T3:\n≥24 days"
    ))

dfi$risk_group2 = factor(dfi$risk_group2, levels = c("T3:\n≥24 days", "T2:\n20-23 days", "T1:\n≤19 days"),  labels = c("≥24 d", "20–23 d", "≤19 d"))

fit1 = survfit(Surv(follow_up_period, infection) ~ risk_group2, data = dfi) 

pval1 <- surv_pvalue(fit1, data = dfi)$pval
p_label1 <- if (pval1 < 0.001) {
  "italic(p) < 0.001"
} else {
  sprintf("italic(p) == %s", formatC(pval1, format = "f", digits = 3))
}

g1 <- ggsurvplot(fit1,
                fun = "event",  
                data = dfi,
                palette = c("#E41A1C", "#377EB8", "#4DAF4A"),  
                size = 1,   
                ggtheme = theme_classic() +
                  theme(legend.title = element_text(face = "bold"),   
                        legend.text  = element_text(face = "plain")),  
                font.main = c(15, "bold", "black"), 
                font.x = c(14, "bold", "black"),  
                font.y = c(14, "bold", "black"),  
                font.tickslab = c(12, "black"),  
                conf.int = FALSE,  
                pval = FALSE, 
                risk.table.fontsize = 4.4,  
                risk.table.pos = "out",  
                tables.y.text = FALSE,  
                tables.theme = theme_cleantable(),  
                risk.table = c("absolute"),  
                risk.table.col = "strata",  
                risk.table.height = 0.25,  
                ylab = "Probability of FN",  
                xlab = "Days",  
                surv.scale = "percent",  
                break.x.by = 3,  
                ylim = c(0, 1), 
                xlim = c(0, 28),  
                censor = TRUE,  
                censor.shape = 108,  
                censor.size = 4,  
                font.legend= c(12, "black"),  
                legend.title = "Duration of neutropenia (IND)",           
                legend.labs = levels(dfi$risk_group2),
                legend = "top")

g1$plot = g1$plot + 
  scale_x_continuous(
    breaks = seq(0, 28, by = 3),
    labels = function(x) x + 1
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           size = 4,
           label = p_label1, parse = TRUE, hjust = 0) +
  guides(color = guide_legend(nrow = 1))

# ---------
# Figure 1H
# ---------

dfk$tertile = cut(
  dfk$neutropenia_days_050,
  breaks = quantile(dfk$neutropenia_days_050, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("low", "int", "high")
)

dfk = dfk %>% 
  dplyr::select(id, infection, follow_up_period, cycle, neutropenia_days_050, tertile) %>%
  dplyr::mutate(
    risk_group2 = case_when(
      tertile == "low" ~ "T1:\n≤9 days",
      tertile == "int" ~ "T2:\n10-13 days",
      tertile == "high" ~ "T3:\n≥14 days"
    ))

dfk$risk_group2 = factor(dfk$risk_group2, levels = c("T3:\n≥14 days", "T2:\n10-13 days", "T1:\n≤9 days"), labels = c("≥14 d", "10–13 d", "≤9 d"))

fit2 = survfit(Surv(follow_up_period, infection) ~ risk_group2, data = dfk) 


pval2 <- surv_pvalue(fit2, data = dfk)$pval
p_label2 <- if (pval2 < 0.001) {
  "italic(p) < 0.001"
} else {
  sprintf("italic(p) == %s", formatC(pval2, format = "f", digits = 3))
}

g2 <- ggsurvplot(fit2,
                fun = "event",  
                data = dfk,
                palette = c("#E41A1C", "#377EB8", "#4DAF4A"),  
                size = 1,   
                ggtheme = theme_classic() +
                  theme(legend.title = element_text(face = "bold"),   
                        legend.text  = element_text(face = "plain")),  
                font.main = c(15, "bold", "black"),  
                font.x = c(14, "bold", "black"), 
                font.y = c(14, "bold", "black"),  
                font.tickslab = c(12, "black"),  
                conf.int = FALSE,  
                pval = FALSE,  
                risk.table.fontsize = 4.4,  
                risk.table.pos = "out", 
                tables.y.text = FALSE, 
                tables.theme = theme_cleantable(), 
                risk.table = c("absolute"),  
                risk.table.col = "strata",  
                risk.table.height = 0.25,  
                ylab = "Probability of FN",  
                xlab = "Days", 
                surv.scale = "percent",  
                break.x.by = 3,  
                ylim = c(0, 1),  
                xlim = c(0, 28),  
                censor = TRUE,  
                censor.shape = 108, 
                censor.size = 4,  
                font.legend = c(12, "black"),  
                legend.title = "Duration of neutropenia (CONS)", 
                legend.labs = levels(dfk$risk_group2),  
                legend = "top")

g2$plot = g2$plot + 
  scale_x_continuous(
    breaks = seq(0, 28, by = 3),
    labels = function(x) x + 1
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           size = 4,
           label = p_label2, parse = TRUE, hjust = 0) +
  guides(color = guide_legend(nrow = 1))

# ----------------------
# Supplementary Figure 1
# ----------------------

a_df = df %>% 
  dplyr::filter(overlaps == TRUE) %>%                                          
  dplyr::group_by(id, treatment_date) %>%        
  distinct(drug, .keep_all = TRUE) %>%           
  ungroup() %>%                                             
  dplyr::mutate(class = case_when(
    str_detect(atc_code, "J02AA") ~ "Polyenes", 
    str_detect(atc_code, "J02AC") ~ "Azoles",
    str_detect(atc_code, "J02AX") ~ "Echinocandins",
    str_detect(atc_code, "J01GB") ~ "Aminoglycosides",
    str_detect(atc_code, "J01C") ~ "Penicillins",
    str_detect(atc_code, c("J01DB|J01DC|J01DD|J01DE|J01DI")) ~ "Cephalosporins",
    str_detect(atc_code, "J01DF") ~ "Monobactams", 
    str_detect(atc_code, "J01DH") ~ "Carbapenems",
    str_detect(atc_code, "J01FA") ~ "Macrolides",
    str_detect(atc_code, "J01FF") ~ "Lincosamides", 
    str_detect(atc_code, "J01XA") ~ "Glycopeptides", 
    str_detect(atc_code, "J01XX08") ~ "Oxazolidinones", 
    str_detect(atc_code, "J01XX09") ~ "Lipopeptides", 
    str_detect(atc_code, "J01XB") ~ "Polymyxins", 
    str_detect(atc_code, "J01AA") ~ "Tetracyclines", 
    str_detect(atc_code, "J01MA") ~ "Fluoroquinolones",
    str_detect(atc_code, "J01EE") ~ "Sulfonamides+\nDiaminopyrimidines",
    str_detect(atc_code, "J01XE") ~ "Nitrofurans", 
    str_detect(atc_code, "J01XC") ~ "Fusidanes", 
    str_detect(atc_code, "J01XD") ~ "Nitroimidazoles", 
    TRUE ~ str_sub(drug, 1,6)) 
  )

ai = a_df %>%  dplyr::filter(cycle == "IND")
ai2 = as.data.frame(table(ai$class)) %>% dplyr::arrange(desc(Freq))
ai3 = ai2 %>% 
  dplyr::mutate(Var1 = if_else(Freq <11, "Others", Var1)) %>% 
  dplyr::group_by(Var1) %>% 
  summarise(Freq = sum(Freq), .groups = "drop") %>% 
  ungroup() %>% 
  dplyr::mutate(group = case_when(
    Var1 %in% c("Polyenes", "Azoles", "Echinocandins") ~ "af",
    TRUE ~ "ab"               
  )) %>% 
  dplyr::rename(Freq_IND = Freq) %>% 
  dplyr::arrange(desc(Freq_IND))

ak = a_both %>% dplyr::filter(sykli == "CONS")
ak2 = as.data.frame(table(ak$class)) %>% dplyr::arrange(desc(Freq))
ak3 = ak2 %>% 
  dplyr::mutate(Var1 = if_else(Freq <16, "Others", Var1)) %>% 
  dplyr::group_by(Var1) %>% 
  summarise(Freq = sum(Freq), .groups = "drop") %>% 
  ungroup() %>% 
  dplyr::mutate(group = case_when(
    Var1 %in% c("Polyenes", "Azoles", "Echinocandins") ~ "af",
    TRUE ~ "ab"               
  )) %>% 
  dplyr::rename(Freq_KONS = Freq) %>%
  dplyr::arrange(desc(Freq_KONS))

ai3 = ai3 %>% # to ensure the same groups are shown in the mirror plot
  add_row(Var1 = "Sulfonamides+\nDiaminopyrimidines", Freq_IND = 6, group = "ab") %>% 
  dplyr::mutate(Freq_IND = ifelse(Var1 == "Others", 6, Freq_IND), 
                Perc_IND = (Freq_IND/206) * 100)

ak3 = ak3 %>% 
  add_row(Var1 = "Glycopeptides", Freq_KONS = 9, group = "ab") %>% 
  dplyr::mutate(Freq_KONS = ifelse(Var1 == "Others", 21, Freq_KONS), 
                Perc_KONS = (Freq_KONS/440) * 100)

a = full_join(ai3, ak3, by = c("Var1", "group"))
a = a %>% 
  dplyr::mutate(Total = Freq_IND + Freq_KONS) %>% 
  dplyr::arrange(Total)

a_mirror = a %>%                                                         
  pivot_longer(cols = c(Perc_IND, Perc_KONS),                            
               names_to = "Phase", 
               values_to = "Value") %>% 
  mutate(Value = ifelse(Phase == "Perc_IND", -Value, Value),             
         group = ifelse(Var1 == "Others", "muu", group),
         Var1 = factor(Var1, levels = a$Var1))

a_mirror$Phase = factor(a_mirror$Phase, levels = c("Perc_IND", "Perc_KONS"), labels = c("IND", "CONS")) 

plot1 = 
  ggplot(a_mirror, aes(x = Var1, y = Value, fill = Phase)) +                       
  geom_col() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(-100,100),
                     breaks = seq(-100,100, by = 20),
                     labels = function(x) abs(x)
  ) +
  labs(x = "Antimicrobial Agents", 
       y = "Percentage of cycles (%)",
       fill = NULL) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 11, face = "bold", margin = margin(t = 12)),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank()
  )

# ----------------------
# Supplementary Figure 2
# ----------------------

g <- ggplot(df1, aes(x = variable, y = neutropenia_days_050, group = genetic_variable1, fill = genetic_variable1)) +
  geom_jitter(width = 0.2,
              size = 0.2,
              alpha = 0.5,
              color = "black") +
  geom_boxplot(varwidth = TRUE,
               alpha = 0.7,
               color = "black",
               outlier.shape = NA,
               notch = FALSE, 
               notchwidth = 0.5,
               outlier.color = "red",
               outlier.size = 2
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("NEG", "POS")),
                     p.adjust.method = "BH",
                     paired = FALSE,
                     label = "p.signif"
  ) +
  labs(
    y = "Neutropenia days",
    x = "Genetic Variable1"
  ) +
  scale_x_discrete(labels = n) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black", face = "bold"),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 14, face = "bold", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )

# ----------------------
# Supplementary Figure 3
# ----------------------

g <- ggplot(df1, aes(x = genetic_variable1, y = time_to_first_neutropenia, group = genetic_variable1, fill = genetic_variable1)) +
  geom_jitter(width = 0.2,
              size = 0.2,
              alpha = 0.5,
              color = "black") +
  geom_boxplot(varwidth = TRUE,
               alpha = 0.7,
               color = "black",
               outlier.shape = NA,
               notch = FALSE, 
               notchwidth = 0.5,
               outlier.color = "red",
               outlier.size = 2
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("NEG", "POS")),
                     p.adjust.method = "BH",
                     paired = FALSE,
                     label = "p.signif") +
  labs(
    y = "Onset of neutropenia (days)",
    x = "Genetic Variable1"
  ) +
  scale_x_discrete(labels = n) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black", face = "bold"),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 14, face = "bold", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )

bloodcount_variable1 <- ggplot(df1, 
               aes(
                 x = days, 
                 y = bloodcount_variable1,
                 color = genetic_variable1                                                                      
               )) +
  geom_smooth(method = "loess", 
              se=TRUE,         
              alpha = 0.25
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  guides(color = guide_legend(title = "Genetic Variable1"))  + 
  labs(x = "Days", 
       y = "Genetic Variable1 10⁹/L") +                            
  scale_x_continuous(breaks = seq(1,28, by = 3)) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank()
  )

# ----------------------
# Supplementary Figure 4
# ----------------------

g <- ggplot(df1, aes(x = variable, y = neutropenia_days_050, group = variable, fill = variable)) +
  geom_jitter(width = 0.2,
              size = 0.2,
              alpha = 0.5,
              color = "black") +
  geom_boxplot(varwidth = TRUE,
               alpha = 0.7,
               color = "black",
               outlier.shape = NA,
               notch = FALSE, 
               notchwidth = 0.5,
               outlier.color = "red",
               outlier.size = 2
   ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("A", "B")),
                     p.adjust.method = "BH",
                     paired = FALSE,
                     label = "p.signif") +
  labs(
    y = "Neutropenia days",
    x = "variable"
  ) +
  scale_x_discrete(labels = n) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black", face = "bold"),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 14, face = "bold", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )

fit1 = survfit(Surv(follow_up_period, infection) ~ etiology2, data = df) 

pval1 <- surv_pvalue(fit1, data = df)$pval
p_label1 <- if (pval1 < 0.001) {
  "italic(p) < 0.001"
} else {
  sprintf("italic(p) == %s", formatC(pval1, format = "f", digits = 3))
}

g <- ggsurvplot(fit1,
                fun = "event",  
                data = df,
                palette = c("#377EB8", "#E41A1C"),   
                size = 1,   
                ggtheme = theme_classic() +
                  theme(legend.title = element_text(face = "bold"),  
                        legend.text  = element_text(face = "plain")),  
                font.main = c(15, "bold", "black"),  
                font.x = c(14, "bold", "black"),  
                font.y = c(14, "bold", "black"),  
                font.tickslab = c(12, "black"),  
                conf.int = FALSE,  
                pval = FALSE, 
                risk.table.fontsize = 4.4,  
                risk.table.pos = "out",  
                tables.y.text = FALSE,  
                tables.theme = theme_cleantable(),  
                risk.table = c("absolute"),  
                risk.table.col = "strata",  
                risk.table.height = 0.25,  
                ylab = "Probability of FN",  
                xlab = "Days",  
                surv.scale = "percent",  
                break.x.by = 3,  
                ylim = c(0, 1),  
                xlim = c(0, 28),  
                censor = TRUE,  
                censor.shape = 108,  
                censor.size = 4,  
                font.legend= c(12, "black"),  
                legend.title = "Etiology", 
                legend.labs = levels(df$etiology2),  
                legend = "top")

g$plot = g$plot + 
  scale_x_continuous(
    breaks = seq(0, 28, by = 3),
    labels = function(x) x + 1
  ) +
  annotate("text",
           x = 0.1, y = 0.1,
           size = 4,
           label = p_label1, parse = TRUE, hjust = 0) +
  guides(color = guide_legend(nrow = 1)); g

