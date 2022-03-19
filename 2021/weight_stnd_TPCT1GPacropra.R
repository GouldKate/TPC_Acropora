library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

mydata_M50Cd<-read.csv("C:/Github/TPC_Acropora/R_output_Data_Frames/mydata_M50Cd.csv")

treat1<-mydata_M50Cd%>%
  filter(treatment=="T1")
treat2<-mydata_M50Cd%>%
  filter(treatment=="T2")
treat3<-mydata_M50Cd%>%
  filter(treatment=="T3")
treatC<-mydata_M50Cd%>%
  filter(treatment=="C")

treat1GP<-treat1%>%
  filter(rate.type=="GP")
treat2GP<-treat2%>%
  filter(rate.type=="GP")
treat3GP<-treat3%>%
  filter(rate.type=="GP")
treatCGP<-treatC%>%
  filter(rate.type=="GP")

treat1GR<-treat1%>%
  filter(rate.type=="R")
treat2GR<-treat2%>%
  filter(rate.type=="R")
treat3GR<-treat3%>%
  filter(rate.type=="R")
treatCGR<-treatC%>%
  filter(rate.type=="R")

ggplot(treat1GP, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')


ggplot(treat2GP, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')


ggplot(treat3GP, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')


ggplot(treatCGP, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')

ggplot(treat1GR, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')


ggplot(treat2GR, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')


ggplot(treat3GR, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')

ggplot(treatCGR, aes(temp.Cat, log.rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '')

##################################
##################################
data("chlorella_tpc")
d <- filter(chlorella_tpc, curve_id == 1)

d_ave <- filter(chlorella_tpc, process == 'adaptation', growth_temp == 33, flux == 'photosynthesis') %>%
  group_by(temp, flux) %>%
  dplyr::summarise(., sd = sd(rate),
            ave_rate = mean(rate)) %>%
  ungroup()

ggplot() +
  geom_linerange(aes(x = temp, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2)



treat1GPSSt = summarySE(treat1GP, measurevar="log.rate",
                       groupvars=c("temp.Cat", "treatment"))
treat1GPSSt

ggplot() +
  geom_linerange(aes(x = temp.Cat, ymin = log.rate - sd, ymax = log.rate + sd), treat1GPSSt) +
  geom_point(aes(temp.Cat, log.rate), treat1GPSSt, size = 2, shape = 21, fill = 'green4') +
  theme_classic() +
  labs(x ='Temperature (ºC)',
       y = 'Gross Photosynthesis',
       title = 'Sprinkler Shading') +
  geom_hline(aes(yintercept = 0), linetype = 2)


treat1GPSStg = summarySE(treat1GP, measurevar="log.rate",
                        groupvars=c("temp.Cat", "genotype"))
treat1GPSStg


ggplot() +
  geom_linerange(aes(x = temp.Cat, ymin = log.rate - sd, ymax = log.rate + sd, group=genotype), treat1GPSStg) +
  geom_point(aes(temp.Cat, log.rate), treat1GPSStg, size = 2, shape = 21, fill = 'green4') +
  theme_classic() +
  labs(x ='Temperature (ºC)',
       y = 'Gross Photosynthesis',
       title = 'Sprinkler Shading') +
  geom_hline(aes(yintercept = 0), linetype = 2)

d_ave <- treat1GP%>%
  group_by(temp.Cat, treatment) %>%
  dplyr::summarise(., sd = sd(log.rate),
                   ave_rate = mean(log.rate)) %>%
  ungroup()

ggplot() +
  geom_linerange(aes(x = temp.Cat, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp.Cat, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x ='Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2)
# fit kamykowski model using rTPC with and without weights
d_fits <- nest(d_ave, data = c(temp.Cat, ave_rate, sd)) %>%
  mutate(standard = map(data, ~nls_multstart(ave_rate~kamykowski_1985(temp = temp.Cat, tmin, tmax, a,b,c),
                                             data = .x,
                                             iter = c(4,4,4,4,4),
                                             start_lower = get_start_vals(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985') - 10,
                                             start_upper = get_start_vals(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985') + 10,
                                             lower = get_lower_lims(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985'),
                                             upper = get_upper_lims(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         weighted = map(data, ~nls_multstart(ave_rate~kamykowski_1985(temp = temp.Cat, tmin, tmax, a,b,c),
                                             data = .x,
                                             iter = c(4,4,4,4,4),
                                             start_lower = get_start_vals(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985') - 10,
                                             start_upper = get_start_vals(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985') + 10,
                                             lower = get_lower_lims(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985'),
                                             upper = get_upper_lims(.x$temp.Cat, .x$ave_rate, model_name = 'kamykowski_1985'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE,
                                             # include weights here!
                                             modelweights = 1/sd)))
d_fits

calc_params(d_fits)
# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', standard:weighted)

# get predictions using augment
newdata <- tibble(temp = seq(min(d_ave$temp.Cat), max(d_ave$temp.Cat), length.out = 8))

d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# get topt for each model
d_topt <- mutate(d_stack, topt = map_dbl(fit, get_topt), rmax = map_dbl(fit, get_rmax),
                 topt_text = paste(topt, 'ºC', sep = ' '))

# plot
library(ggrepel)

ggplot(d_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  geom_linerange(aes(x = temp.Cat, ymin = ave_rate - sd, ymax = ave_rate + sd), d_ave) +
  geom_point(aes(temp.Cat, ave_rate), d_ave, size = 2, shape = 21, fill = 'green4') +
  geom_label_repel(aes(topt, rmax, label = topt_text, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_topt) +
  geom_point(aes(topt, rmax, col = model_name), size = 4, d_topt) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Photosynthesis rates across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2) +
  ylim(c(0, 1))

####################################################
tttest <- treat1GP%>%
  group_by(temp.Cat, treatment, genotype) %>%
  dplyr::summarise(., sd = sd(log.rate),
                   ave_rate = mean(log.rate)) %>%
  ungroup()

#####################################################
#####################################################
#####################################################


library(boot)
library(car)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(patchwork)
library(minpack.lm)

data("bacteria_tpc")

d <- filter(bacteria_tpc, phage == 'nophage')

d<-treat1GP%>%
  group_by(treatment, temp.Cat, genotype)%>%
dplyr::summarise(., sd = sd(log.rate),
                 ave_rate = mean(log.rate)) %>%
  ungroup()

 
ggplot(d, aes(temp.Cat, ave_rate)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# fit Sharpe-Schoolfield model
# 
mod = 'sharpschoolhigh_1981'
start_vals <- get_start_vals(d$temp.Cat, d$ave_rate, model_name = 'sharpeschoolhigh_1981')
# get limits
low_lims <- get_lower_lims(d$temp.Cat, d$ave_rate, model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(d$temp.Cat, d$ave_rate, model_name = 'sharpeschoolhigh_1981')

d_fit <- nest(d, data = c(temp.Cat, ave_rate, treatment, genotype, sd)) %>%
  mutate(sharpeschoolhigh = map(data, ~nls_multstart(ave_rate~sharpeschoolhigh_1981(temp = temp.Cat, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(3,3,3,3),
                                                     start_lower = get_start_vals(.x$temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp.Cat, .x$ave_rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         # create new temperature data
         new_data = purrr::map(data, ~tibble(temp = seq(min(.x$temp.Cat), max(.x$temp.Cat), length.out = 96))),
         # predict over that data,
         preds =  map2(sharpeschoolhigh, new_data, ~augment(.x, newdata = .y)))


                  # unnest predictions
          d_preds <- dplyr::select(d_fit, preds) %>%
            unnest_legacy(preds)


          # plot data and predictions
          ggplot() +
            geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
            geom_point(aes(temp.Cat, ave_rate), d, size = 2, alpha = 0.5) +
            theme_bw(base_size = 12) +
            labs(x = 'Temperature (ºC)',
                 y = 'Growth rate',
                 title = 'Growth rate across temperatures')
          
          
          # refit model using nlsLM
          fit_nlsLM <- minpack.lm::nlsLM(ave_rate~sharpeschoolhigh_1981(temp = temp.Cat, r_tref,e,eh,th, tref = 15),
                                         data = d,
                                         start = coef(d_fit$sharpeschoolhigh[[1]]),
                                         lower = get_lower_lims(d$temp.Cat, d$ave_rate, model_name = 'sharpeschoolhigh_1981'),
                                         upper = get_upper_lims(d$temp.Cat, d$ave_rate, model_name = 'sharpeschoolhigh_1981'),
                                         weights = rep(1, times = nrow(d)))
          
          
          
          # bootstrap using case resampling
          boot1 <- Boot(fit_nlsLM, method = 'case')
          
          head(boot1$t)

          hist(boot1, layout = c(2,2))
          
          
          # create predictions of each bootstrapped model
          boot1_preds <- boot1$t %>%
            as.data.frame() %>%
            drop_na() %>%
            mutate(iter = 1:n()) %>%
            group_by_all() %>%
            do(data.frame(temp = seq(min(d$temp.Cat), max(d$temp.Cat), length.out = 100))) %>%
            ungroup() %>%
            mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))
          
          
          # calculate bootstrapped confidence intervals
          boot1_conf_preds <- group_by(boot1_preds, temp) %>%
            dplyr::summarise(conf_lower = quantile(pred, 0.025),
                      conf_upper = quantile(pred, 0.975)) %>%
            ungroup()
          # plot bootstrapped CIs
          p1 <- ggplot() +
            #geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
            geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
            geom_point(aes(temp.Cat, ave_rate), d, size = 2, alpha = 0.5) +
            theme_bw(base_size = 12) +
            labs(x = 'Temperature (ºC)',
                 y = 'Growth rate',
                 title = 'Growth rate across temperatures')
          
          
            
    start_vals
low_lims
upper_lims



fit <- nls_multstart(log.rate~sharpeschoolhigh_1981(temp = temp.Cat, r_tref,e,eh,th, tref = 15),
                     data = treat1GP,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')

fit

calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)


# predict new data
new_data <- data.frame(temp = seq(min(treat1GP$temp.Cat), max(treat1GP$temp.Cat), length.out = 96))
preds <- broom::augment(fit, newdata = new_data)

# plot data and model fit
# # plot data and model fit
ggplot(treat1GP, aes(temp.Cat, log.rate)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate')


ggplot(d, aes(temp, rate)) +
  geom_point() +
  geom_line(aes(temp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# fit Sharpe-Schoolfield model
T2GP_fit <-  group_by(treat2GP,treatment) %>%
  nest()%>%
  mutate(sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp.Cat, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(3,3,3,3),
                                                     start_lower = get_start_vals(.x$temp.Cat, .x$rate, model_name = 'sharpeschoolhigh_1981') - 1,
                                                     start_upper = get_start_vals(.x$temp.Cat, .x$rate, model_name = 'sharpeschoolhigh_1981') + 1,
                                                     lower = get_lower_lims(.x$temp.Cat, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp.Cat, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp.Cat), max(.x$temp.Cat), length.out = 100))),
         # predict over that data,
         preds =  map2(sharpeschoolhigh, new_data, ~augment(.x, newdata = .y)))


# fit every model formulation in rTPC
T1GP_fits <- nest(treat1GP, data = c(temp.Cat, log.rate)) %>%
  mutate(beta = map(data, ~nls_multstart(log.rate~beta_2012(temp = temp.Cat, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'beta_2012') - 1,
                                         start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'beta_2012') + 1,
                                        # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'beta_2012'),
                                        # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         boatman = map(data, ~nls_multstart(log.rate~boatman_2017(temp = temp.Cat, rmax, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4,4),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'boatman_2017') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'boatman_2017') + 1,
                                           # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'boatman_2017'),
                                           # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'boatman_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         briere2 = map(data, ~nls_multstart(log.rate~briere2_1999(temp = temp.Cat, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'briere2_1999') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'briere2_1999') + 1,
                                           # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'briere2_1999'),
                                           # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'briere2_1999'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         delong = map(data, ~nls_multstart(log.rate~delong_2017(temp = temp.Cat, c, eb, ef, tm, ehc),
                                           data = .x,
                                           iter = c(4,4,4,4,4),
                                           start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'delong_2017') - 1,
                                           start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'delong_2017') + 1,
                                          # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'delong_2017'),
                                          # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'delong_2017'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         flinn = map(data, ~nls_multstart(log.rate~flinn_1991(temp = temp.Cat, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'flinn_1991') - 1,
                                          start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'flinn_1991') + 1,
                                         #lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'flinn_1991'),
                                         #upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(log.rate~gaussian_1987(temp = temp.Cat, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'gaussian_1987') - 1,
                                             start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'gaussian_1987') + 1,
                                            # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'gaussian_1987'),
                                            # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         hinshelwood = map(data, ~nls_multstart(log.rate~hinshelwood_1947(temp = temp.Cat, a, e, b, eh),
                                                data = .x,
                                                iter = c(5,5,5,5),
                                                start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'hinshelwood_1947') - 1,
                                                start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'hinshelwood_1947') + 1,
                                               # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'hinshelwood_1947'),
                                               # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'hinshelwood_1947'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
         joehnk = map(data, ~nls_multstart(log.rate~joehnk_2008(temp = temp.Cat, rmax, topt, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4,4, 4),
                                           start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'joehnk_2008') - 1,
                                           start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'joehnk_2008') + 1,
                                          # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'joehnk_2008'),
                                          # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'joehnk_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         johnson_lewin = map(data, ~suppressWarnings(nls_multstart(log.rate~ johnsonlewin_1946(temp = temp.Cat, r0, e, eh, topt),
                                                                   data = .x,
                                                                   iter = c(4,4,4,4),
                                                                   start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'johnsonlewin_1946') - 1,
                                                                   start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'johnsonlewin_1946') + 1,
                                                                  # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'johnsonlewin_1946'),
                                                                  # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'johnsonlewin_1946'),
                                                                   supp_errors = 'Y',
                                                                   convergence_count = FALSE))),
         kamykowski = map(data, ~nls_multstart(log.rate~kamykowski_1985(temp = temp.Cat, tmin, tmax, a, b, c),
                                               data = .x,
                                               iter = c(4,4,4,4,4),
                                               start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'kamykowski_1985') - 1,
                                               start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'kamykowski_1985') + 1,
                                              #lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'kamykowski_1985'),
                                              #upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin2 = map(data, ~nls_multstart(log.rate~lactin2_1995(temp = temp.Cat, a, b, tmax, delta_t),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'lactin2_1995') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'lactin2_1995') + 1,
                                           # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'lactin2_1995'),
                                           # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'lactin2_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         modifiedgaussian = map(data, ~nls_multstart(log.rate~modifiedgaussian_2006(temp = temp.Cat, rmax, topt, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'modifiedgaussian_2006') - 1,
                                                     start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'modifiedgaussian_2006') + 1,
                                                    # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'modifiedgaussian_2006'),
                                                    # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(log.rate~oneill_1972(temp = temp.Cat, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'oneill_1972') - 1,
                                           start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'oneill_1972') + 1,
                                          # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'oneill_1972'),
                                          # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         pawar = map(data, ~nls_multstart(log.rate~pawar_2018(temp = temp.Cat, r_tref, e, eh, topt, tref = 15),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'pawar_2018') - 1,
                                          start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'pawar_2018') + 1,
                                         # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'pawar_2018'),
                                         # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'pawar_2018'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(log.rate~quadratic_2008(temp = temp.Cat, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'quadratic_2008') + 0.5,
                                             # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'quadratic_2008'),
                                             # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = map(data, ~nls_multstart(log.rate~ratkowsky_1983(temp = temp.Cat, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'ratkowsky_1983') - 1,
                                              start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'ratkowsky_1983') + 1,
                                             # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'ratkowsky_1983'),
                                             # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         rezende = map(data, ~nls_multstart(log.rate~rezende_2019(temp = temp.Cat, q10, a,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'rezende_2019') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'rezende_2019') + 1,
                                          # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'rezende_2019'),
                                          # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'rezende_2019'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         sharpeschoolfull = map(data, ~nls_multstart(log.rate~sharpeschoolfull_1981(temp = temp.Cat, r_tref,e,el,tl,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolfull_1981') - 1,
                                                     start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolfull_1981') + 1,
                                                    # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolfull_1981'),
                                                    # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolfull_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(log.rate~sharpeschoolhigh_1981(temp = temp.Cat, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981') - 1,
                                                     start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981') + 1,
                                                    # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981'),
                                                    # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoollow = map(data, ~nls_multstart(log.rate~sharpeschoollow_1981(temp = temp.Cat, r_tref,e,el,tl, tref = 15),
                                                    data = .x,
                                                    iter = c(4,4,4,4),
                                                    start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoollow_1981') - 1,
                                                    start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoollow_1981') + 1,
                                                    #lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoollow_1981'),
                                                    #upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoollow_1981'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         spain = map(data, ~nls_multstart(log.rate~spain_1982(temp = temp.Cat, a,b,c,r0),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'spain_1982') - 1,
                                          start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'spain_1982') + 1,
                                         #lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'spain_1982'),
                                         #upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'spain_1982'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         thomas1 = map(data, ~nls_multstart(log.rate~thomas_2012(temp = temp.Cat, a,b,c,topt),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2012') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2012') + 2,
                                            #lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2012'),
                                            #upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2012'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         thomas2 = map(data, ~nls_multstart(log.rate~thomas_2017(temp = temp.Cat, a,b,c,d,e),
                                            data = .x,
                                            iter = c(3,3,3,3,3),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2017') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2017') + 1,
                                           # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2017'),
                                           # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'thomas_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(log.rate~weibull_1995(temp = temp.Cat, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'weibull_1995') - 1,
                                            start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'weibull_1995') + 1,
                                           # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'weibull_1995'),
                                           # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))






d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(beta = map(data, ~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         briere2 = map(data, ~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         delong = map(data, ~nls_multstart(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                           data = .x,
                                           iter = c(4,4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         flinn = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                          data = .x,
                                          iter = c(5,5,5),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         hinshelwood = map(data, ~nls_multstart(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
                                                data = .x,
                                                iter = c(5,5,5,5),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') - 1,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') + 1,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
         joehnk = map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4,4, 4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         johnson_lewin = map(data, ~suppressWarnings(nls_multstart(rate~ johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                                   data = .x,
                                                                   iter = c(4,4,4,4),
                                                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') - 1,
                                                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') + 1,
                                                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                                   upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                                   supp_errors = 'Y',
                                                                   convergence_count = FALSE))),
         kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                               data = .x,
                                               iter = c(4,4,4,4,4),
                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') - 10,
                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') + 10,
                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
         lactin2 = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         modifiedgaussian = map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                              data = .x,
                                              iter = c(4,4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') - 10,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') + 10,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         rezende = map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         sharpeschoolfull = map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         sharpeschoollow = map(data, ~nls_multstart(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                    data = .x,
                                                    iter = c(4,4,4,4),
                                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') - 10,
                                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') + 10,
                                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)),
         spain = map(data, ~nls_multstart(rate~spain_1982(temp = temp, a,b,c,r0),
                                          data = .x,
                                          iter = c(4,4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') - 1,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') + 1,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
         thomas1 = map(data, ~nls_multstart(rate~thomas_2012(temp = temp, a,b,c,topt),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') - 1,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') + 2,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         thomas2 = map(data, ~nls_multstart(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                            data = .x,
                                            iter = c(3,3,3,3,3),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = temp, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))



# fit Sharpe-Schoolfield model
T1GP_fit <- nest(treat1GP, data = c(temp.Cat, log.rate)) %>%
  mutate(sharpeschoolhigh = map(data, ~nls_multstart(log.rate~sharpeschoolhigh_1981(temp = temp.Cat, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(3,3,3,3),
                                                     start_lower = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981') - 1,
                                                     start_upper = get_start_vals(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981') + 1,
                                                    # lower = get_lower_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981'),
                                                    # upper = get_upper_lims(.x$temp.Cat, .x$log.rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp.Cat), max(.x$temp.Cat), length.out = 100))),
         # predict over that data,
         preds =  map2(sharpeschoolhigh, new_data, ~augment(.x, newdata = .y)))


# keep just a single curve
d <- filter(bacteria_tpc, phage == 'nophage')
       
d_fit <- nest(d, data = c(temp, rate)) %>%
  mutate(sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(3,3,3,3),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(sharpeschoolhigh, new_data, ~augment(.x, newdata = .y)))

ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# If we want confidence bands around this prediction, 
# we can get those by resampling the data a number of times.
# The R package car contains the function Boot() that provides
# a wrapper for the widely used function boot::boot() that is tailored to bootstrapping regression models.
# 
# 
#nls_multstart() is designed to fit models across a wide possible parameter space, but as it samples multiple start parameters for each model, using it with bootstrapping becomes computationally expensive. Instead, we refit the model using minpack.lm::nlsLM(), using the coefficients of nls_multstart() as the start values. The Boot() function then refits the model 999 times and stores the model coefficients.

# refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                               data = d,
                               start = coef(d_fit$sharpeschoolhigh[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'sharpeschoolhigh_1981'),
                               weights = rep(1, times = nrow(d)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')
head(boot1$t)


#hist.boot() to look at the distribution of each parameter.

hist(boot1, layout = c(2,2))


create predictions for each of these models and through this confidence intervals around the original fitted predictions. We can then plot (1) the bootstrapped fits and (2) the confidence regions around the model predictions.

# create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# plot bootstrapped predictions
p2 <- ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_line(aes(temp, pred, group = iter), boot1_preds, col = 'blue', alpha = 0.007) +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

p1 + p2