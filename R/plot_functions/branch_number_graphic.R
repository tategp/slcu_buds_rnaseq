# Made up data frame
plasticity <- rep(c("lowplas","highplas"),each = 100)
nitrate <- rep(c("Low Nitrate","High Nitrate","Low Nitrate","High Nitrate"), each = 50) 
count <- as.numeric(c(rnorm(50,2.9,1.5),rnorm(50,3.1,1.5),rnorm(50,1,1.5),rnorm(50,5,1.5)))
df <- data.frame(plasticity,nitrate,count)

nitratelvl <- c("Low Nitrate","High Nitrate")

# Plot the counts
ggplot(df, aes(x = nitrate, y = count, colour = plasticity)) +
  geom_pointrange(stat = "summary", fun.data = "mean_cl_boot", aes(group = plasticity), size = 1) +
  geom_line(stat = "summary", fun.y = "mean", aes(group = plasticity), size = 2) +
  scale_x_discrete(limits = nitratelvl) +
  labs(y = "Branch number") +
  scale_colour_manual(values=c("deepskyblue4", "red"), 
                      name  ="Plasticity", 
                      labels=c("Non-Plastic", "Plastic"))+
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = c(0.2,0.8)) 
