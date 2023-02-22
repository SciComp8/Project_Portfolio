dataset = dat.expr.full2390.traintop9
var.x = "BMI"
var.y = "ENSG00000164691.12"


dataset |> ggplot(mapping = aes(y = get(var.y), x = get(var.x))) + 
  geom_boxplot(mapping = aes(color = get(var.x))) +
  geom_jitter(mapping = aes(color = get(var.x))) + 
  stat_compare_means(method = "wilcox.test", 
                     label.x.npc = "right",
                     label.y = 2,
                     mapping = aes(label = ifelse(p < 1.e-3, "p < 0.001", sprintf("p = %4.3f", as.numeric(..p.format..))))) +
  facet_wrap(~factor(Data, levels = c("Train", "Test"))) + 
  labs(y = var.y, 
       x = "") +
  theme_bw(base_size = 14) + 
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF")) + 
  scale_x_discrete(limits = c("low", "high")) + 
  guides(color = "none")
