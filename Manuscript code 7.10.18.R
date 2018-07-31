##Code for Joao

##Import data

##Note for some of the descriptives I used SPSS just because I already had the data in that format,
#so I don't have the code in here. 

rm(list = ls())

library(lattice)
library(ggplot2)
library(stats)

#change WD

setwd("~/Desktop/TrustRSA/")

#data for main analyses

data <- read.table(paste("Trust_n48_6.19.18.txt",sep=""), header=FALSE)
colnames(data) <-c("Participant","Age","Sex", "Amygdala", "AmygdalaCentromedial", "AmygdalaLaterobasal", "AmygdalaSuperficial", "NeutralProportion")
str(data)

#create new df to center age
datarev <- as.data.frame(cbind(data$Age))


center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}
data2 <- center_colmeans(datarev)
str(data2)

colnames(data2) <-c("Age")

#Main dataset used in analyses is data3
data3 <- as.data.frame(cbind(data$Participant, data2$Age, data$Sex, data$Amygdala, data$AmygdalaCentromedial, data$AmygdalaLaterobasal, data$AmygdalaSuperficial, data$NeutralProportion))

#change participant ID to numerical order for visualization in charts
data3[data3=="101"] <- 1
data3[data3=="102"] <- 2
data3[data3=="103"] <- 3
data3[data3=="105"] <- 4
data3[data3=="106"] <- 5
data3[data3=="107"] <- 6
data3[data3=="108"] <- 7
data3[data3=="109"] <- 8
data3[data3=="110"] <- 9
data3[data3=="111"] <- 10
data3[data3=="113"] <- 11
data3[data3=="114"] <- 12
data3[data3=="115"] <- 13
data3[data3=="116"] <- 14
data3[data3=="117"] <- 15
data3[data3=="118"] <- 16
data3[data3=="119"] <- 17
data3[data3=="120"] <- 18
data3[data3=="121"] <- 19
data3[data3=="122"] <- 20
data3[data3=="123"] <- 21
data3[data3=="125"] <- 22
data3[data3=="126"] <- 23
data3[data3=="127"] <- 24
data3[data3=="128"] <- 25
data3[data3=="129"] <- 26
data3[data3=="130"] <- 27
data3[data3=="131"] <- 28
data3[data3=="132"] <- 29
data3[data3=="134"] <- 30
data3[data3=="135"] <- 31
data3[data3=="137"] <- 32
data3[data3=="138"] <- 33
data3[data3=="140"] <- 34
data3[data3=="141"] <- 35
data3[data3=="142"] <- 36
data3[data3=="143"] <- 37
data3[data3=="144"] <- 38
data3[data3=="145"] <- 39
data3[data3=="146"] <- 40
data3[data3=="147"] <- 41
data3[data3=="148"] <- 42
data3[data3=="149"] <- 43
data3[data3=="150"] <- 44
data3[data3=="151"] <- 45
data3[data3=="152"] <- 46
data3[data3=="153"] <- 47
data3[data3=="156"] <- 48

colnames(data3) <-c("Participant","Age","Sex", "Amygdala", "AmygdalaCentromedial", "AmygdalaLaterobasal", "AmygdalaSuperficial", "NeutralProportion")
str(data3)

#summary stats
summary(data3$NeutralProportion)
summary(data$Age)
summary(data$Sex)

#regression models using poly for orthogonalization
library(stats)
library(caret)
library(broom)
library(nlme)

#CUBIC MODEL BEST FIT amy anatomical. Naming conventions => linear = linear model, model = cubic model
linear <- lm(data3$NeutralProportion ~data3$Amygdala, data=data3)
summary(linear)

model<-lm(data3$NeutralProportion ~poly(data3$Amygdala, 3))
summary(model)

anova(linear, model)

glance(linear)
glance(model)

#rmse for cubic model (rmse for k-fold is pushed out in the caret script)
mse<-mean(residuals(model)^2)
rmse<-sqrt(mse)
rmse

#model diagnostics
hist(residuals(model), col="darkgray", xlab="residuals(Amygdala Anatomical Cubic Model)", main="Histogram of Residuals",cex.lab=1.5, cex.main=1.5, cex.axis=1.2, ylim=range(0,20))
plot(fitted(model), 
     residuals(model), xlab="fitted", ylab="residuals", main="Amygdala Anatomical Cubic Model", cex.lab=1.5, cex.main=1.5,cex.axis=1.2)


#LINEAR MODEL BEST FIT amygdala centromedial 
linear2 <- lm(data3$NeutralProportion ~data3$AmygdalaCentromedial, data=data3)
summary(linear2)
model2<-lm(data3$NeutralProportion ~poly(data3$AmygdalaCentromedial, 3))
summary(model2)

#trying quadratic since cubic was n.s.
model2.2<-lm(data3$NeutralProportion ~poly(data3$AmygdalaCentromedial, 2))
summary(model2.2)

anova(linear2,model2)

hist(residuals(linear2), col="darkgray", xlab="residuals(Amygdala Centromedial Linear Model)", main="Histogram of Residuals",cex.lab=1.5, cex.main=1.5, cex.axis=1.2, ylim=range(0,20))
plot(fitted(linear2), 
     residuals(linear2), xlab="fitted", ylab="residuals", main="Amygdala Centromedial Linear Model", cex.lab=1.5, cex.main=1.5,cex.axis=1.2)

mse2<-mean(residuals(linear2)^2)
rmse2<-sqrt(mse2)
rmse2

#CUBIC MODEL BEST FIT amy laterobasal
linear3 <- lm(data3$NeutralProportion ~data3$AmygdalaLaterobasal, data=data3)
summary(linear3)
model3<-lm(data3$NeutralProportion ~poly(data3$AmygdalaLaterobasal, 3))
summary(model3)

anova(linear3,model3)

hist(residuals(model3), col="darkgray", xlab="residuals(Amygdala Laterobasal Cubic Model)", main="Histogram of Residuals",cex.lab=1.5, cex.main=1.5, cex.axis=1.2, ylim=range(0,20))
plot(fitted(model3), 
     residuals(model3), xlab="fitted", ylab="residuals", main="Amygdala Laterobasal Cubic Model", cex.lab=1.5, cex.main=1.5,cex.axis=1.2)

mse3<-mean(residuals(model3)^2)
rmse3<-sqrt(mse3)
rmse3

#CUBIC MODEL BEST FIT amy superficial

linear4 <- lm(data3$NeutralProportion ~data3$AmygdalaSuperficial, data=data3)
summary(linear4)
model4<-lm(data3$NeutralProportion ~poly(data3$AmygdalaSuperficial, 3))
summary(model4)

anova(linear4,model4)

hist(residuals(model4), col="darkgray", xlab="residuals(Amygdala Superficial Cubic Model)", main="Histogram of Residuals",cex.lab=1.5, cex.main=1.5, cex.axis=1.2, ylim=range(0,20))
plot(fitted(model4), 
     residuals(model4), xlab="fitted", ylab="residuals", main="Amygdala Superficial Cubic Model", cex.lab=1.5, cex.main=1.5,cex.axis=1.2)

mse4<-mean(residuals(model4)^2)
rmse4<-sqrt(mse4)
rmse4

# K-fold cross-validation

#install.packages("caret")
library(caret)

# Define train control for k fold cross validation
set.seed(48)
train_control <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions = TRUE, verboseIter = TRUE)

modelk <- train(NeutralProportion~poly(Amygdala, 3), data=data3, trControl=train_control, method="lm")
modelk2 <- train(NeutralProportion~poly(AmygdalaCentromedial,1), data=data3, trControl=train_control, method="lm")
modelk3 <- train(NeutralProportion~poly(AmygdalaLaterobasal,3), data=data3, trControl=train_control, method="lm")
modelk4 <- train(NeutralProportion~poly(AmygdalaSuperficial,3), data=data3, trControl=train_control, method="lm")
# Summarise Results
#RMSE - root mean squared error (std of residuals. vertical direction is y variable direction) has units of y associated
#Rsquared (proportion, higher better, no units associated)
#MAE - mean absolute error
print(modelk)
print(modelk2)
print(modelk3)
print(modelk4)

#AGE AND SEX NS
linearas <- lm(NeutralProportion ~Amygdala +Age + Sex, data=data3)
summary(linearas)
linear2as <- lm(NeutralProportion ~AmygdalaCentromedial+Age + Sex, data=data3)
summary(linear2as)
linear3as <- lm(NeutralProportion ~AmygdalaLaterobasal+Age + Sex, data=data3)
summary(linear3as)
linear4as <- lm(NeutralProportion ~AmygdalaSuperficial+Age + Sex, data=data3)
summary(linear4as)

#QUAD NS
quad <- lm(NeutralProportion ~Amygdala + AmygdalaSq, data=data3)
summary(quad)
quad2 <- lm(NeutralProportion ~AmygdalaCentromedial + AmygdalaCentromedialSq, data=data3)
summary(quad2)
quad3 <- lm(NeutralProportion ~AmygdalaLaterobasal + AmygdalaLaterobasalSq,  data=data3)
summary(quad3)
quad4 <- lm(NeutralProportion ~AmygdalaSuperficial + AmygdalaSuperficialSq, data=data3)
summary(quad4)

#FITTED LOESS LINES

bold.20.text <- element_text(face = "bold", color = "black", size = 20)
text.14.text <- element_text(color = "black", size = 14)


formula=data3$NeutralProportion ~poly(data3$Amygdala, 3, raw=TRUE)
p<-ggplot(data=data3,aes(Amygdala,NeutralProportion))
p<-p+geom_point(colour="black")
p<-p + geom_smooth(method="loess",se=TRUE, colour= "black") + theme_classic() 
p <- p + geom_smooth(method="lm", se=TRUE, fill = NA, formula = data3$NeutralProportion~poly(data3$Amygdala,3,raw=TRUE)) 
p <- p+ labs(x = "Amygdala Anatomical Similarity", y="Neutral Discrimination Score")
p <- p + theme(axis.title = bold.20.text, axis.text=text.14.text)
p

formula2=data3$NeutralProportion ~poly(data3$AmygdalaCentromedial, 1, raw=TRUE)
p2<-ggplot(data=data3,aes(AmygdalaCentromedial,NeutralProportion))
p2<-p2+geom_point(colour="black")
p2<-p2 + geom_smooth(method="loess",se=TRUE, colour= "black") + theme_classic() 
p2 <- p2 + geom_smooth(method="lm", se=TRUE, fill = NA, formula = data3$NeutralProportion~poly(data3$AmygdalaCentromedial,1,raw=TRUE)) 
p2 <- p2+ labs(x = "Amygdala Centromedial Similarity", y="Neutral Discrimination Score")
p2 <- p2 + theme(axis.title = bold.20.text, axis.text=text.14.text)
p2

formula3=data3$NeutralProportion ~poly(data3$AmygdalaLaterobasal, 3, raw=TRUE)
p3<-ggplot(data=data3,aes(AmygdalaLaterobasal,NeutralProportion))
p3<-p3+geom_point(colour="black")
p3<-p3 + geom_smooth(method="loess",se=TRUE, colour= "black") + theme_classic() 
p3 <- p3 + geom_smooth(method="lm", se=TRUE, fill = NA, formula = data3$NeutralProportion~poly(data3$AmygdalaLaterobasal,3,raw=TRUE)) 
p3 <- p3+ labs(x = "Amygdala Laterobasal Similarity", y="Neutral Discrimination Score")
p3 <- p3 + theme(axis.title = bold.20.text, axis.text=text.14.text)
p3

formula4=data3$NeutralProportion ~poly(data3$AmygdalaSuperficial, 3, raw=TRUE)
p4<-ggplot(data=data3,aes(AmygdalaSuperficial,NeutralProportion))
p4<-p4+geom_point(colour="black")
p4<-p4 + geom_smooth(method="loess",se=TRUE, colour= "black") + theme_classic() 
p4 <- p4 + geom_smooth(method="lm", se=TRUE, fill = NA, formula = data3$NeutralProportion~poly(data3$AmygdalaSuperficial,3,raw=TRUE)) 
p4 <- p4+ labs(x = "Amygdala Superficial Similarity", y="Neutral Discrimination Score")
p4 <- p4 + theme(axis.title = bold.20.text, axis.text=text.14.text)
p4

###########PLOTS WITH PREDICTED MODELS#############

data6<-as.data.frame(cbind(data3$Participant, data3$Amygdala, data3$NeutralProportion))
colnames(data6) <-c("Participant","Similarity", "Discrimination")
mod1<-lm(data6$Discrimination ~poly(data6$Similarity, 3))

data7<- as.data.frame(cbind(data3$Participant,data3$AmygdalaCentromedial, data3$NeutralProportion))
colnames(data7) <-c("Participant","Similarity", "Discrimination")
mod2<-lm(data7$Discrimination ~poly(data7$Similarity, 1))

data8<-as.data.frame(cbind(data3$Participant,data3$AmygdalaLaterobasal, data3$NeutralProportion))
colnames(data8) <-c("Participant","Similarity", "Discrimination")
mod3<-lm(data8$Discrimination ~poly(data8$Similarity, 3))

data9<-as.data.frame(cbind(data3$Participant,data3$AmygdalaSuperficial, data3$NeutralProportion))
colnames(data9) <-c("Participant","Similarity", "Discrimination")
mod4<-lm(data9$Discrimination ~poly(data9$Similarity, 3))

p<- ggplot(data6, aes(y=Discrimination, x=Similarity)) + 
  geom_point(alpha = .5) + 
  stat_smooth(method = "lm", formula = y ~ poly(x,3), colour= "black") + theme_classic() 
p<- p+ labs(x = "Amygdala Anatomical Similarity", y="Neutral Discrimination Score") + theme(axis.title = bold.20.text, axis.text=text.14.text)
p

p2<- ggplot(data7, aes(y=Discrimination, x=Similarity)) + 
  geom_point(alpha = .5) + 
  stat_smooth(method = "lm", formula = y ~ poly(x,1), colour= "black") + theme_classic() 
p2<- p2+ labs(x = "Amygdala Centromedial Similarity", y="Neutral Discrimination Score") + theme(axis.title = bold.20.text, axis.text=text.14.text)
p2

p3<- ggplot(data8, aes(y=Discrimination, x=Similarity)) + 
  geom_point(alpha = .5) + 
  stat_smooth(method = "lm", formula = y ~ poly(x,3), colour= "black") + theme_classic() 
p3<- p3+ labs(x = "Amygdala Laterobasal Similarity", y="Neutral Discrimination Score") + theme(axis.title = bold.20.text, axis.text=text.14.text)
p3


p4<- ggplot(data9, aes(y=Discrimination, x=Similarity)) + 
  geom_point(alpha = .5) + 
  stat_smooth(method = "lm", formula = y ~ poly(x,3), colour= "black") + theme_classic() 
p4<- p4+ labs(x = "Amygdala Superficial Similarity", y="Neutral Discrimination Score") + theme(axis.title = bold.20.text, axis.text=text.14.text)
p4

##################below is code for average activation for faces with 12 score (untrustworthy) and 67 score (trustworthy)
##Note: Discrimination4 is discrimination for neutral faces whereas "Discrimination" is overall for the task and 
#Discriminationtrustworthy is faces =67 combined whereas Discriminatonuntrustworthy is 12 faces combined

rm(list = ls())

library(lattice)
library(ggplot2)
library(stats)

setwd("~/Desktop/TrustRSA/")

data3 <- read.table(paste("trust_RSA_average_data.txt",sep=""), header=FALSE)
colnames(data3) <-c("Participant","Discrimination1","Discrimination2", "Discrimination3",
                    "Discrimination4","Discrimination5","Discrimination6","Discrimination7",
                    "DiscriminationUntrustworthy","DiscriminationTrustworthy","Discrimination",
                    "Amygdala12", "AmygdalaCentromedial12", "AmygdalaLaterobasal12", "AmygdalaSuperficial12", 
                    "Amygdala67", "AmygdalaCentromedial67", "AmygdalaLaterobasal67", "AmygdalaSuperficial67")


library(stats)
summary(data3$Discrimination)

#########################12 faces (untrustworthy)

#n.s. amygdala anatomical
linear <- lm(Discrimination4 ~Amygdala12, data=data3)
summary(linear)
model<-lm(data3$Discrimination4 ~poly(data3$Amygdala12, 3))
summary(model)

#n.s. amygdala centromedial
linear2 <- lm(Discrimination4 ~AmygdalaCentromedial12, data=data3)
summary(linear2)
model2<-lm(data3$Discrimination4 ~poly(data3$AmygdalaCentromedial12, 3))
summary(model2)
model2.2<-lm(data3$Discrimination4 ~poly(data3$AmygdalaCentromedial12, 2))
summary(model2.2)

#n.s. amygdala laterobasal
linear3 <- lm(Discrimination4 ~AmygdalaLaterobasal12, data=data3)
summary(linear3)
model3<-lm(data3$Discrimination4 ~poly(data3$AmygdalaLaterobasal12, 3))
summary(model3)

#n.s. amygdala superficial

linear4 <- lm(Discrimination4 ~AmygdalaSuperficial12, data=data3)
summary(linear4)
model4<-lm(data3$Discrimination4 ~poly(data3$AmygdalaSuperficial12, 3))
summary(model4)

##########################67 faces (trustworthy)

#n.s. amygdala anatomical
linear <- lm(Discrimination4 ~Amygdala67, data=data3)
summary(linear)
model<-lm(data3$Discrimination4 ~poly(data3$Amygdala67, 3))
summary(model)

#n.s. amygdala centromedial
linear2 <- lm(Discrimination4 ~AmygdalaCentromedial67, data=data3)
summary(linear2)
model2<-lm(data3$Discrimination4 ~poly(data3$AmygdalaCentromedial67, 3))
summary(model2)
model2.2<-lm(data3$Discrimination4 ~poly(data3$AmygdalaCentromedial67, 2))
summary(model2.2)

#n.s. amygdala laterobasal
linear3 <- lm(Discrimination4 ~AmygdalaLaterobasal67, data=data3)
summary(linear3)
model3<-lm(data3$Discrimination4 ~poly(data3$AmygdalaLaterobasal67, 3))
summary(model3)

#n.s. amygdala superficial

linear4 <- lm(Discrimination4 ~AmygdalaSuperficial67, data=data3)
summary(linear4)
model4<-lm(data3$Discrimination4 ~poly(data3$AmygdalaSuperficial67, 3))
summary(model4)
