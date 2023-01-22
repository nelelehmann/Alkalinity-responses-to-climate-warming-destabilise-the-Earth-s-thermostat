#------------------------------------------------------------------------------#
# This script is supplementary to:                                             #
# Alkalinity responses to climate warming destabilise the Earth's thermostat   #
# by Nele Lehmann, Tobias Stacke, Sebastian Lehmann, Hugues Lantuit,           #
# John Gosse, Chantal Mears, Jens Hartmann and Helmuth Thomas                  #
#------------------------------------------------------------------------------#

#### Dependent libraries ####
library(data.table)
library(openxlsx)
library(dplyr)
library(mgcv) # For GAM estimation
library(boot) # For cross validation (cv)
library(ggplot2)

#### Data preparation ####
dt <- read.xlsx("data/Lehmann-etal_2022_v6_NL.xlsx", startRow = 225) %>% 
  as.data.table()

# Adapting variable names to shorter ones 
old_names = c("AT.norm.[mol/m**4/a]", 
              "Sum.[%].(carbonate.containing.rock.typ...)", 
              "MAT.[°C].(Based.on.WorldClim.2)",
              "Erosion.rate.[mm/ka].(Based.on.OCTOPUS.data.base)",
              "Soil.thick.[cm].(Based.on.Global.depth.to.bedr...)",
              "Catch.area.[km**2]",
              "MAP.[mm].(Based.on.WorldClim.2)",
              "Elevation.[m.a.s.l.].(Based.on.OCTOPUS.data.base)",
              "Dams.extent.[%].(GOODD.(global.dataset.of.more...)",
              "Metam.[%].(Based.on.Global.lithological....)",
              "Sediments.[%].(unconsolidated,.Based.on.Glob...)",
              "Siliciclastic.sed.rocks.[%].(Based.on.Global.lithological....)",
              "Evaporite.[%].(Based.on.Global.lithological....)",
              "Pyroclastics.[%].(Based.on.Global.lithological....)",
              "Volc.rocks,.acidic.[%].(Based.on.Global.lithological....)",
              "Volc.rocks,.basic.[%].(Based.on.Global.lithological....)",
              "Volc.rocks,.intermed.[%].(Based.on.Global.lithological....)",
              "Plutonic.rocks,.acidic.[%].(Based.on.Global.lithological....)",
              "Plutonic.rocks,.basic.[%].(Based.on.Global.lithological....)",
              "Plutonic.rocks,.intermed.[%].(Based.on.Global.lithological....)",
              "Perm.snow.+.ice.[%].(Based.on.GlobCover)",
              "Mean.slope.gradient.[m/km].(Based.on.OCTOPUS.data.base)"
              )

new_names = c("alk_norm",
              "carbo_sum",
              "temp_ann",
              "ero_rate",
              "soil_thick",
              "area_2D",
              "precip_ann",
              "elev",
              "dams",
              "mt",
              "su",
              "ss",
              "ev",
              "py",
              "va",
              "vb",
              "vi",
              "pa",
              "pb",
              "pi",
              "snow",
              "slope"
              )

for(i in 1:length(old_names)){
  names(dt)[names(dt) == old_names[i]] <- new_names[i]
}

# Change soil regolith thickness from centimeter to meter
dt$soil_thick <- dt$soil_thick / 100

#### Model estimation ####
model_list <- list()

##### Models for comparison #####
model_list$m1 <-
  gam(alk_norm ~ carbo_sum,
      family = gaussian(link = "log"),
      data = dt)
model_list$m2 <- update(model_list$m1, . ~ . + poly(temp_ann, 3), 
                        family = gaussian(link = "log"))
model_list$m3 <- update(model_list$m2, . ~ . + poly(log(ero_rate), 2),
                        family = gaussian(link = "log"))
model_list$m4 <- update(model_list$m3, . ~ . + soil_thick + log(area_2D),
                        family = gaussian(link = "log"))

#####  Reference model ##### 
model_list$gam_reference <- gam(
  alk_norm ~
    s(carbo_sum) +
    s(temp_ann) +
    s(log(ero_rate)) +
    s(soil_thick) +
    s(log(area_2D)),
  family = gaussian(link = "log"),
  data = dt
)


#####  Final model ##### 
model_list$m5 <-
  gam(
    alk_norm ~ carbo_sum +
      poly(temp_ann, 5) +
      poly(log(ero_rate), 2) +
      soil_thick +
      log(area_2D),
    family = gaussian(link = "log"),
    data = dt
  )

##### Models with MAP for comparison #####
# Model m5 with MAP instead of MAT
model_list$m6 <- 
  gam(
  alk_norm ~ carbo_sum +
    precip_ann +
    poly(log(ero_rate), 2) +
    soil_thick +
    log(area_2D),
  family = gaussian(link = "log"),
  data = dt
)
# Model m5 with both MAT and MAP
model_list$m7 <- 
  gam(
    alk_norm ~ carbo_sum +
      precip_ann +
      poly(temp_ann, 5) +
      poly(log(ero_rate), 2) +
      soil_thick +
      log(area_2D),
    family = gaussian(link = "log"),
    data = dt
  )

#### Model outputs ####
for(model in model_list) {
  model %>%
    summary %>%
    print
}

##### Raw polynomials ####
# Estimation of m5 with raw polynomials instead of orthogonalized polynomials 
# for better interpretation of the coefficients 
m5_raw_pol <- gam(
  alk_norm ~ carbo_sum +
    poly(temp_ann, 5,raw = T) +
    poly(log(ero_rate), 2,raw = T) +
    soil_thick +
    log(area_2D),
  family = gaussian(link = "log"),
  data = dt
) 

m5_raw_pol %>% summary


#### Model fit ####
# creating a list for AIC, BIC, RSS, adjRSq & 
# mean squared prediction cross validation (cv) error
dt_model_fit <- data.table(
  model = character(),
  AIC = numeric(),
  BIC = numeric(),
  RSS = numeric(),
  adj_R_sq = numeric(),
  mean_sq_prediction_cv_error = numeric()
)

set.seed(42)

for (i in 1:length(model_list)) {
  # Cross validation
  tmp_cv <- cv.glm(data = dt, glmfit = model_list[[i]])
  dt_model_fit <-
    rbind(dt_model_fit,
          list(
            model = names(model_list[i]),
            AIC = round(AIC(model_list[[i]])),
            BIC = round(BIC(model_list[[i]])),
            RSS = (dt$alk_norm - model_list[[i]]$fitted.values) ^ 2 %>% 
              sum  %>%
              round(2),
            adj_R_sq = round(summary(model_list[[i]])$r.sq, 3),
            mean_sq_prediction_cv_error = round(tmp_cv$delta[2], 2)
          ))
}

dt_model_fit

#### Cross validation error for temperature bands ####
# Define temperature bands
tmp_interval <- function(tmp){
  if(tmp<0){
    return("[-3,0[")
  } else if( tmp >= 0 & tmp < 5){
    return("[0,5[")
  } else if( tmp >= 5 & tmp < 10){
    return("[5,10[")
  } else if( tmp >= 10 & tmp < 15){
    return("[10,15[")
  } else if( tmp >= 15 & tmp < 20){
    return("[15,20[")
  } else if( tmp >= 20 & tmp < 27){
    return("[20,27[")
  }
}

# Function for calculating the cv prediction error
get_cv_prediction_error <- 
  function(formula, dt) {
    p_error <- c()
    for (i in 1:nrow(dt)) {
      tmp_m <- gam(
        formula = formula,
        family = gaussian(link = "log"),
        data = dt[-i]
      )
      pred <- predict(tmp_m, newdata = dt[i], type = "response")
      p_error <- rbind(p_error, dt[i]$alk_norm - pred)
    }
    p_error <- p_error %>% as.vector
    return(p_error)
  }

# Apply the cv function to all models
dt_cv <- data.table(
  temp_ann = dt$temp_ann,
  tmp_interval = sapply(dt$temp_ann, tmp_interval),
  cv_pred_error_m1 = get_cv_prediction_error(formula = formula(model_list$m1),
                                             dt = dt),
  cv_pred_error_m2 = get_cv_prediction_error(formula = formula(model_list$m2),
                                             dt = dt),
  cv_pred_error_m3 = get_cv_prediction_error(formula = formula(model_list$m3),
                                             dt = dt),
  cv_pred_error_m4 = get_cv_prediction_error(formula = formula(model_list$m4),
                                             dt = dt),
  cv_pred_error_m5 = get_cv_prediction_error(formula = formula(model_list$m5),
                                             dt = dt),
  cv_pred_error_gam = get_cv_prediction_error(formula = formula(model_list$gam_reference),
                                             dt = dt)
)

# Get the mean squared prediction error for each temperature band
dt_cv[,
      .(cv_msq_prd_error_m1 = round(mean(cv_pred_error_m1 ^ 2), 3),
        cv_msq_prd_error_m2 = round(mean(cv_pred_error_m2 ^ 2), 3),
        cv_msq_prd_error_m3 = round(mean(cv_pred_error_m3 ^ 2), 3),
        cv_msq_prd_error_m4 = round(mean(cv_pred_error_m4 ^ 2), 3),
        cv_msq_prd_error_m5 = round(mean(cv_pred_error_m5 ^ 2), 3),
        cv_msq_prd_error_gam = round(mean(cv_pred_error_gam ^ 2), 3)
      ),
      tmp_interval]

# Overall mean squared prediction error
dt_cv[, lapply(.SD, function(x) round(mean(x ^ 2), 3)), .SDcols = names(dt_cv)[3:8]]

#### Feature importance checks #####
##### 1. Derive importance of each dependent variable by AIC Stepwise Algorithm #####
# Start point: Intercept only
# Scope: Feature set of m5
# Direction: Forward

intercept_only <- glm(
  alk_norm ~ 1,
  data = dt,
  family = gaussian(link = "log")
  )

forward_step_AIC <- step(intercept_only, 
                direction = 'forward', 
                scope = formula(model_list$m5), 
                trace = 0)

forward_step_AIC$anova

##### 2. Declare importance by permutation feature importance test ######

# Number of repetitions
K <- 10000

# Reference score (adj_R_sq) of model m5
s <- summary(model_list$m5)["r.sq"] %>% 
  unlist()

# Feature set 
feature_list <- c("carbo_sum", 
            "temp_ann",
            "ero_rate",
            "soil_thick",
            "area_2D") 

set.seed(42)

get_s_k <- function(feature_corrupted, dt2corrupt, feature, formula) {
  dt2corrupt[[feature]] <- feature_corrupted
  s_k <-
    summary(gam(formula = formula, data = dt2corrupt,
                gaussian(link = "log")))["r.sq"] %>% unlist()
  return(s_k)
}

get_importance4feature <- function(K, dt, formula, feature, s){
  feature_corrupted <- replicate(K, sample(dt[[feature]]))
  s_feature <- apply(
    X = feature_corrupted,
    MARGIN = 2,
    FUN = get_s_k,
    feature = feature,
    dt2corrupt = dt,
    formula = formula
  )
  importance_feature <- s - mean(s_feature)
  names(importance_feature) <- feature
  return(importance_feature)
}

importance_list <- lapply(
  X = feature_list,
  get_importance4feature,
  dt = dt,
  K = K,
  formula = formula(model_list$m5),
  s = s
) %>% 
  unlist()

# Results
importance_list


#### Sensitivity Analysis ####
# Exclude all data points with dam extent bigger than or equal to 10%
dt_nDam <- dt[dams < 10]

m5_nDam <- gam(formula(model_list$m5), 
               data = dt_nDam, 
               family = gaussian(link = "log")
               )
m5_nDam

# How does adding dams to the model change the fit
m_dam <-
  gam(
    alk_norm ~ carbo_sum +
      poly(temp_ann, 5, raw = T) +
      poly(log(ero_rate), 2, raw = T) +
      soil_thick +
      log(area_2D) +
      dams,
    family = gaussian(link = "log"),
    data = dt
  )
m_dam %>% summary()

#### Miscellaneous ####
##### Check significance of other rock types when added to m5 #####
rock_types <- c("mt",
                "su",
                "ss",
                "ev",
                "py",
                "va",
                "vb",
                "vi",
                "pa",
                "pb",
                "pi")

pv_rock_type <- adj_R_sq_rock_type <- c()

for(rock_type in rock_types){
  tmp <- update(model_list$m5, paste0(". ~ . + ", rock_type),
                family = gaussian(link = "log")) %>% summary
  pv_rock_type <- c(pv_rock_type, tmp$p.pv[names(tmp$p.pv) == rock_type])
  names(tmp$r.sq) <- rock_type
  adj_R_sq_rock_type <- c(adj_R_sq_rock_type, tmp$r.sq)
}

pv_rock_type
adj_R_sq_rock_type


#### Plot/fit erosion rate vs area and slope
ggplot(dt, aes(log(area_2D), log(ero_rate))) +
  geom_point() +
  geom_smooth() +
  geom_smooth(method = "lm", col = "red")

lm(log(ero_rate) ~ log(area_2D),dt) %>% summary()

ggplot(dt, aes(slope, log(ero_rate))) +
  geom_point() +
  geom_smooth() +
  geom_smooth(method = "lm", col = "red")

lm(log(ero_rate) ~ slope, dt) %>% summary()

ggplot(dt, aes(`Mean.slope.gradient.[m/km].(Based.on.OCTOPUS.data.base)` * log(area_2D), log(ero_rate))) +
  geom_point() +
  geom_smooth() +
  geom_smooth(method = "lm", col = "red")

lm(log(ero_rate) ~ log(area_2D)*slope, dt) %>% summary()

##### Check if feature importance of alkalinity flux differs from normalized alkalinity #####

dt$alk_flux <- dt$`AT.[µmol/l]` * dt$`Runoff.[mm].(mean.annual,.Based.on.UNH/GRD...)`

alk_flux_gam <- gam(
  alk_flux ~
    s(carbo_sum) +
    s(temp_ann) +
    s(log(ero_rate)) +
    s(soil_thick) +
    s(log(area_2D)),
  family = gaussian(link = "log"),
  data = dt
)

alk_flux_gam %>% 
  summary


###### permutation test ######


# Number of repetitions WARNING computation takes ~20h
K <- 500

# Reference score (adj_R_sq) of model m5
s <- summary(alk_flux_gam)["r.sq"] %>% 
  unlist()

# Feature set 
feature_list <- c("carbo_sum", 
                  "temp_ann",
                  "ero_rate",
                  "soil_thick",
                  "area_2D") 

set.seed(42)

importance_list <- lapply(
  X = feature_list,
  get_importance4feature,
  dt = dt,
  K = K,
  formula = formula(alk_flux_gam),
  s = s
) %>% 
  unlist()

# Results
importance_list



