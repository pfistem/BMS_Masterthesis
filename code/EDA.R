
library(readxl)
data <- read_excel("data/data.xlsx")



data[data$STAGE == 1,]

data$PARAMCD <- as.factor(data$PARAMCD)

summary(data)

filtered_data <- data[data$STAGE == 1,]
filtered_data <- filtered_data[filtered_data$PARAMCD == "PFS",]
summary(filtered_data)


boxplot(sqrt(AVAL)~ARM, data= filtered_data)

table(filtered_data$ARM)
table(filtered_data$USUBJID)

anova(lm(sqrt(AVAL)~ARM, data=filtered_data))



# Stage2 ------------------------------------------------------------------


filtered_data2 <- data[data$STAGE == 2,]
filtered_data2 <- filtered_data2[filtered_data2$PARAMCD == "PFS",]
summary(filtered_data2)


boxplot(sqrt(AVAL)~ARM, data= filtered_data2)

table(filtered_data2$ARM)
table(filtered_data2$USUBJI2$PARAMCD == "PFS",]
summary(filtered_data2)


boxplot(sqrt(AVAL)~ARM, data= filtered_data2)

table(filtered_data2$ARM)
table(filtered_data2$USUBJID)


# Stage 1 & 2 -------------------------------------------------------------


filtered_data3 <- data[data$PARAMCD == "PFS",]
filtered_data3 <- filtered_data3[filtered_data3$ARM < 3,]
summary(filtered_data3)

table(filtered_data3$USUBJID, filtered_data3$STAGE)














