#read on data

#based on 
#https://github.com/r-glennie/occuR
#https://r-glennie.github.io/occuR/
#remotes::install_github("r-glennie/occuR")
library(occuR)
library(tidyverse)

#myfolder

visit_data <- readRDS("splines/JAGS/listlengthDF_MTB_Sdanae.rds")

#add coordinates
load("mtbqsDF.RData")
visit_data$x <- mtbqsDF$x_MTB[match(visit_data$MTB, mtbqsDF$Value)]
visit_data$y <- mtbqsDF$y_MTB[match(visit_data$MTB, mtbqsDF$Value)]

#make coordinates smaller
visit_data$x <- visit_data$x/10000
visit_data$y <- visit_data$y/1000000
#multi-occasion occupancy

#correlation of site occupancy over space and time is induced by allowing occupancy probability to be a smooth function of space and time

#need to make visit data with columns: site, occassion and obs
visit_data <- visit_data[,c("MTB","Year","visit","Species",
                            "singleList","yday","CoarseNaturraum","x","y")]
names(visit_data)[1:4] <- c("site","occasion","visit","obs")

#need to make vist be indexed from i to n within each site and occasion
visit_data <- visit_data %>%
              group_by(site, occasion) %>%
              mutate(visit = as.numeric(as.factor(visit)))%>%
              ungroup()

visit_data$occasion <- as.numeric(as.factor(visit_data$occasion))

#need to make site data with "site" and "occasion"
site_data <- unique(visit_data[,c("site","occasion","CoarseNaturraum","x","y")])

#basic model
m0 <- fit_occu(list(psi ~ 1, p ~ 1), as.data.table(visit_data), as.data.table(site_data))
m0

#fixed effect model
m_s <- fit_occu(list(psi ~ CoarseNaturraum, p ~ singleList),
                     as.data.table(visit_data), as.data.table(site_data))
m_s

#year effects
m_t <- fit_occu(list(psi ~ -1 + factor(occasion), p ~ occasion + singleList),
                     as.data.table(visit_data), as.data.table(site_data))
m_t

modelSummary <- data.frame(parameter = names(m_t$res$par.fixed),
                           estimate = plogis(as.numeric(m_t$res$par.fixed)),
                           index = 1: length(names(m_t$res$par.fixed)))

psiSummary <- subset(modelSummary, parameter = "beta_psi")

#extract and plot predictions
qplot(index, estimate, data = psiSummary, geom = "line")

#model with spline
#only splines with basis of “cs” or “ts” are well defined for this package

#one dimention
m_spline <- fit_occu(list(psi ~ s(x,bs = "cs"), p ~ 1),
                     as.data.table(visit_data), as.data.table(site_data))
m_spline

#two dimension
m_spline2d <- fit_occu(list(psi ~ t2(x,y,bs = "ts", k=10), p ~ 1),
                     as.data.table(visit_data), as.data.table(site_data))
m_spline2d
#k=5 worked pretty well
#k=10 more wiggly.
#k=15 strange lines

#predictions
siteInfo_NAs <- readRDS("splines/siteInfo_NAs.rds")
#make coordinates smaller
siteInfo_NAs$x <- siteInfo_NAs$x_MTB/10000
siteInfo_NAs$y <- siteInfo_NAs$y_MTB/1000000

pred_xy <- predict(m_spline2d,
                   as.data.table(visit_data),
                   data.table(occasion = 1, x = siteInfo_NAs$x, y = siteInfo_NAs$y),
                   nboot = 1000)

summary(pred_xy$psi)
siteInfo_NAs$preds <- pred_xy$psi[,1]

ggplot(siteInfo_NAs) +
  geom_point(aes(x = x, y = y, colour = preds)) +
  theme_bw() +
  scale_colour_viridis_c("Occupancy")


# spatio-temporal effect

m_spline3d <- fit_occu(list(psi ~ t2(x, y, occasion, bs = c("ts", "cs"), k=c(5,2)), p ~ 1),
                       as.data.table(visit_data), as.data.table(site_data))


xgr <- rep(gr[,1], nocc)
ygr <- rep(gr[,2], nocc)
tgr <- rep(1:nocc, each = nrow(gr))
pred_xyt <- predict(m_spline2d, visit_data, data.table(occasion = tgr, x = xgr, y = ygr, hab = "arable"), nboot = 1000)

ggplot(data.frame(x = xgr, y = ygr, t = tgr, psi = pred_xyt$psi)) +
  geom_tile(aes(x = x, y = y, group = t, fill = psi)) +
  theme_bw() +
  facet_wrap(~t) +
  scale_x_continuous("x") +
  scale_y_continuous("y") +
  scale_fill_viridis_c("Occupancy")

