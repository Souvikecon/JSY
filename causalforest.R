library(haven)
library(ggplot2)
library(grf)
library(tidyverse)  # ggplot(), %>%, mutate(), and friends
library(aod)
library(estimatr)

set.seed(4596)

keep = c("birthweight","JSY","BPL","SC","ST","OBC","Hindu",
         "Muslim","age_mother","education_mother","birth_order",
         "wealth_index","rural","mother_height","female_child",
         "deadchild","pregnancy_complications","tetanus_injections",
         "iron_tablets", "antenatal_care","antenatal_visits",
         "prenatal_care_person","registered_card","female_head",
         "age_household_head","LPS")

df=read_dta("D:/JSY/Data/JSY_filtered.dta")%>% select(all_of(keep)) %>% na.omit()
covariate_names= c("BPL","SC","ST","OBC","Hindu",
                   "Muslim","age_mother","education_mother","birth_order",
                   "wealth_index","rural","mother_height","female_child",
                   "deadchild","pregnancy_complications","tetanus_injections",
                   "iron_tablets", "antenatal_care","antenatal_visits",
                   "prenatal_care_person","registered_card","female_head",
                   "age_household_head","LPS")

X=df[,3:26] %>% as.matrix()
Y=log(df$birthweight) %>% as.matrix()
W=df$JSY %>% as.matrix()


cf = causal_forest(X, Y , W,
                   tune.parameters = 'all')


hist(e.hat <- cf$W.hat) # Distribution of propensity scores

test.fit=test_calibration(cf) # Run best linear predictor analysis
test.fit

tau.hat = predict(cf)$predictions
hist(tau.hat) # Distribution of CATE

tau=predict(cf,estimate.variance = T)
tau= tau %>% mutate(sd=sqrt(variance.estimates)) %>%
  mutate(lower_ci=predictions-1.96*sd,
         upper_ci=predictions+1.96*sd) #dataframe of tauhats with CIs

df$predictions=tau.hat

ATE = average_treatment_effect(cf,target.sample="all")
ATT= average_treatment_effect(cf,target.sample="treated")

var_imp <- variable_importance(cf)
sorted_var_imp = sort(var_imp, decreasing=TRUE)
ranked.vars <- order(varimp, decreasing = TRUE)
# Top 5 variables according to this measure
colnames(X)[ranked.vars]


varimp_df= data.frame(varimp)
varimp_df$variables= names(df[,3:26])

p <- ggplot(varimp_df,aes(x = reorder(variables,+varimp), y = varimp)) + 
  geom_col(width = 0.7) + 
  xlab("Variables")+
  ylab("Variable importance")+
  theme_classic()
p+coord_flip()
ggsave(paste0(plotsdir,"/varimp.png"),dpi="print",width=10,height=7)


boxplot= function (data,X,Y){
  ggplot(data,aes(y ={{Y}}, x = factor({{X}}))) + 
    geom_boxplot()+
    ylab("CATE")+
    theme_classic()
}

boxplot(df,wealth_index,predictions)+xlab("Wealth index")
ggsave(paste0(plotsdir,"/Wealth index.png"),dpi="print",width=10,height=7)

boxplot(df,education_mother,predictions)+xlab("Mother's education")
ggsave(paste0(plotsdir,"/Mother's education.png"),dpi="print",width=15,height=7)

boxplot(df[df$antenatal_visits<=10,],antenatal_visits,predictions)+xlab("Antenatal_visits")
ggsave(paste0(plotsdir,"/Antenatal visits.png"),dpi="print",width=10,height=7)

df %>% group_by(mother_height) %>% summarise(mean_cate=mean(predictions)) %>%
  ggplot(aes(y =mean_cate, x = mother_height)) + 
  geom_point()+
  ylab("Avg CATE")+
  theme_classic()+
  xlab("Mother's height")
ggsave(paste0(plotsdir,"/Mother's height.png"),dpi="print",width=15,height=7)

boxplot(df,tetanus_injections,predictions)+xlab("Tetanus injections")
ggsave(paste0(plotsdir,"/Tetanus injections.png"),dpi="print",width=10,height=7)

boxplot(df,age_household_head,predictions)+xlab("Age of head")
ggsave(paste0(plotsdir,"/Age of head.png"),dpi="print",width=20,height=7)

boxplot(df,age_mother,predictions)+xlab("Mother's age")
ggsave(paste0(plotsdir,"/Mother's age.png"),dpi="print",width=20,height=7)

boxplot(df,prenatal_care_person,predictions)+xlab("Prenatal care person")
ggsave(paste0(plotsdir,"/Prenatal care person.png"),dpi="print",width=10,height=7)

boxplot(df,rural,predictions)+xlab("Rural")
ggsave(paste0(plotsdir,"/rural.png"),dpi="print",width=10,height=7)

boxplot(df,female_child,predictions)+xlab("female child")
ggsave(paste0(plotsdir,"/femalechild.png"),dpi="print",width=10,height=7)

boxplot(df,antenatal_care,predictions)+xlab("antenatal_care")
ggsave(paste0(plotsdir,"/antenatal care.png"),dpi="print",width=10,height=7)

boxplot(df,BPL,predictions)+xlab("BPL")
ggsave(paste0(plotsdir,"/BPL.png"),dpi="print",width=10,height=7)

boxplot(df,female_head,predictions)+xlab("female head")
ggsave(paste0(plotsdir,"/female head.png"),dpi="print",width=10,height=7)


heatmap= function(data,x,y,fill){
  ggplot(data,aes(x = {{x}} , y = {{y}}, fill = {{fill}})) +
    scale_fill_gradient(high = "red", low = "white")+
    geom_tile()
}

jpeg("heat_WI_edu.jpg",width=1200,height=800,res=200)
df %>% 
  ggplot(aes(x = wealth_index , y = education_mother, fill = predictions)) +
  scale_fill_gradient(high = "red", low = "white")+
  geom_tile() ## heatmap WI & education of mother
dev.off()

jpeg("heat_agemother_edu.jpg",width=1200,height=800,res=200)
df %>% 
  ggplot(aes(x = age_mother, y = education_mother, fill = predictions)) +
  scale_fill_gradient(high = "red", low = "white")+
  geom_tile() ## heatmap WI & education of mother
dev.off()

jpeg("heat_WI_antenatalvisits.jpg",width=1200,height=800,res=200)
df %>% filter(antenatal_visits<=10) %>%
  ggplot(aes(x = wealth_index , y = antenatal_visits, fill = predictions)) +
  scale_fill_gradient(high = "red", low = "white")+
  geom_tile() ## heatmap WI & antenatal visits
dev.off()

jpeg("heat_WI_tetanus.jpg",width=1200,height=800,res=200)
df %>% 
  ggplot(aes(x = wealth_index , y = tetanus_injections, fill = predictions)) +
  scale_fill_gradient(high = "red", low = "white")+
  geom_tile() ## heatmap WI & tetanus injections
dev.off()