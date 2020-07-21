Analysis of Changes in European Bumblebee Populations
=====================================================
Data can be found [here](https://drive.google.com/file/d/1TI7AL9AuLfAAXhklD11tqRPZtmKUV219/view?usp=sharing), [here](https://drive.google.com/file/d/1iMTjrwg-pjHXMzHaRe5vhVxlieKuPfyV/view?usp=sharing), and [here](https://drive.google.com/file/d/103u8Prq8OZcyK2olCwcSm50XQCOcqO3v/view?usp=sharing).

    sink("output.txt")
    data = read.csv("DryadSubmission/data_bombclimate.csv")
    data$period = ''
    data[which(data$periodFrom == 1901),]$period = "historic"
    data[which(data$periodFrom == 1975),]$period = "early"
    data[which(data$periodFrom == 1987),]$period = "middle"
    data[which(data$periodFrom == 1999),]$period = "late"
    data$period = as.factor(data$period)
    data = na.omit(data) 

    historic_EU = data[which(data$period == "historic" & data$continent == "EUR"),]
    early_EU = data[which(data$period == "early" & data$continent == "EUR"),]
    middle_EU = data[which(data$period == "middle" & data$continent == "EUR"),]
    late_EU = data[which(data$period == "late" & data$continent == "EUR"),]

    elevs = c(mean(historic_EU$elevation), mean(early_EU$elevation),mean(middle_EU$elevation), mean(late_EU$elevation))
    lats = c(mean(historic_EU$kmNorthEquator), mean(early_EU$kmNorthEquator), mean(middle_EU$kmNorthEquator), mean(late_EU$kmNorthEquator))
    minTemps = c(mean(historic_EU$minPeriodAnnualMeanT), mean(early_EU$minPeriodAnnualMeanT), mean(middle_EU$minPeriodAnnualMeanT), mean(late_EU$minPeriodAnnualMeanT))
    maxTemps = c(mean(historic_EU$maxPeriodAnnualMeanT), mean(early_EU$maxPeriodAnnualMeanT), mean(middle_EU$maxPeriodAnnualMeanT), mean(late_EU$maxPeriodAnnualMeanT))


    ### Simulations 
    par(mfrow=c(4, 2))
    historic_dist_sim = rnorm(10000, mean(historic_EU$kmNorthEquator), sd(historic_EU$kmNorthEquator))
    plot(density(historic_EU$kmNorthEquator), main = "Historic EU Actual Distribution")
    plot(density(historic_dist_sim), main = "Historic EU Simulated Distribution")
    early_dist_sim = rnorm(10000, mean(early_EU$kmNorthEquator), sd(early_EU$kmNorthEquator))
    plot(density(early_EU$kmNorthEquator), main = "Early EU Actual Distribution")
    plot(density(early_dist_sim), main = "Early EU Simulated Distribution")
    middle_dist_sim = rnorm(10000, mean(middle_EU$kmNorthEquator), sd(middle_EU$kmNorthEquator))
    plot(density(middle_EU$kmNorthEquator), main = "Middle EU Actual Distribution")
    plot(density(middle_dist_sim), main = "Middle EU Simulated Distribution")
    late_dist_sim = rnorm(10000, mean(late_EU$kmNorthEquator), sd(late_EU$kmNorthEquator))
    plot(density(late_EU$kmNorthEquator), main = "Late EU Actual Distribution")
    plot(density(late_dist_sim), main = "Late EU Simulated distribution")

    cat("K-S tests: \n")
    print(ks.test(historic_EU$kmNorthEquator, early_EU$kmNorthEquator)$p.value)
    print(ks.test(historic_dist_sim, early_dist_sim)$p.value)
    print(ks.test(early_EU$kmNorthEquator, middle_EU$kmNorthEquator)$p.value)
    print(ks.test(middle_dist_sim, early_dist_sim)$p.value)
    print(ks.test(middle_EU$kmNorthEquator, late_EU$kmNorthEquator)$p.value)
    print(ks.test(late_dist_sim, middle_dist_sim)$p.value)

    par(mfrow = c(3, 1))
    plot(density((early_dist_sim - historic_dist_sim)), main = "Difference Between Early and Historic")
    lines(density(rnorm(10000, 0, sd(early_dist_sim - historic_dist_sim))), col = 'red')
    plot(density((middle_dist_sim - early_dist_sim)), main = "Difference Between Middle and Early")
    lines(density(rnorm(10000, 0, sd(middle_dist_sim - early_dist_sim))), col = 'red')
    plot(density((late_dist_sim - middle_dist_sim)), main = "Difference Between Late and Middle")
    lines(density(rnorm(10000, 0, sd(late_dist_sim - middle_dist_sim))), col = 'red')

    stat_sim = function(data, INDEX){
        mu = mean(data[INDEX])
        return(mu)
    }
    cat("Tests with simulations: \n")
    boot_sim_late = boot((late_dist_sim - middle_dist_sim), stat = stat_sim, R = 1000)
    boot_sim_middle = boot((middle_dist_sim - early_dist_sim), stat = stat_sim, R = 1000)
    boot_sim_early = boot((early_dist_sim - historic_dist_sim), stat = stat_sim, R = 1000)
    cat(t.test(boot_sim_late$t, boot_sim_middle$t, alternative = "greater")$p.value, '\n')
    cat(wilcox.test(boot_sim_late$t, boot_sim_middle$t, alternative = "greater")$p.value, '\n')
    cat(t.test(boot_sim_middle$t, boot_sim_early$t, alternative = "greater")$p.value, '\n')
    cat(wilcox.test(boot_sim_middle$t, boot_sim_early$t, alternative = "greater")$p.value, '\n')
    ### Analysis
    ### Compute Proportion of Differences
    ###################################################################################################
    elev_pvals_hist = c()
    elev_pvals_early = c()
    elev_pvals_mid = c()

    cat("Elevation: \n")
    for(species in unique(historic_EU$species)){
        elev_pvals_hist = c(elev_pvals_hist, ks.test(historic_EU[which(historic_EU$species == species),]$elevation, early_EU[which(early_EU$species == species),]$elevation)$p.value)
        elev_pvals_early = c(elev_pvals_early, ks.test(early_EU[which(early_EU$species == species),]$elevation, middle_EU[which(middle_EU$species == species),]$elevation)$p.value)
        elev_pvals_mid = c(elev_pvals_mid, ks.test(late_EU[which(late_EU$species == species),]$elevation, middle_EU[which(middle_EU$species == species),]$elevation)$p.value)
    }
    cat(sum(p.adjust(elev_pvals_hist, n = length(elev_pvals_hist)) < 0.05)/length(elev_pvals_hist), '\n')
    cat(sum(p.adjust(elev_pvals_early, n = length(elev_pvals_early)) < 0.05)/length(elev_pvals_early), '\n')
    cat(sum(p.adjust(elev_pvals_mid, n = length(elev_pvals_mid)) < 0.05)/length(elev_pvals_mid), '\n')

    cat("Distance: \n")
    dist_pvals_hist = c()
    dist_pvals_early = c()
    dist_pvals_mid = c()
    for(species in unique(historic_EU$species)){
        dist_pvals_hist = c(dist_pvals_hist, ks.test(historic_EU[which(historic_EU$species == species),]$kmNorthEquator, early_EU[which(early_EU$species == species),]$kmNorthEquator)$p.value)
        dist_pvals_early = c(dist_pvals_early, ks.test(early_EU[which(early_EU$species == species),]$kmNorthEquator, middle_EU[which(middle_EU$species == species),]$kmNorthEquator)$p.value)
        dist_pvals_mid = c(dist_pvals_mid, ks.test(late_EU[which(late_EU$species == species),]$kmNorthEquator, middle_EU[which(middle_EU$species == species),]$kmNorthEquator)$p.value)
    }
    cat(sum(p.adjust(dist_pvals_hist, n = length(dist_pvals_hist)) < 0.05)/length(dist_pvals_hist), '\n')
    cat(sum(p.adjust(dist_pvals_early, n = length(dist_pvals_early)) < 0.05)/length(dist_pvals_early), '\n')
    cat(sum(p.adjust(dist_pvals_mid, n = length(dist_pvals_mid)) < 0.05)/length(dist_pvals_mid), '\n')

    cat("Min: \n")
    min_pvals_hist = c()
    min_pvals_early = c()
    min_pvals_mid = c()
    for(species in unique(historic_EU$species)){
        min_pvals_hist = c(min_pvals_hist, ks.test(historic_EU[which(historic_EU$species == species),]$minPeriodAnnualMeanT, early_EU[which(early_EU$species == species),]$minPeriodAnnualMeanT)$p.value)
        min_pvals_early = c(min_pvals_early, ks.test(early_EU[which(early_EU$species == species),]$minPeriodAnnualMeanT, middle_EU[which(middle_EU$species == species),]$minPeriodAnnualMeanT)$p.value)
        min_pvals_mid = c(min_pvals_mid, ks.test(late_EU[which(late_EU$species == species),]$minPeriodAnnualMeanT, middle_EU[which(middle_EU$species == species),]$minPeriodAnnualMeanT)$p.value)
    }
    cat(sum(p.adjust(min_pvals_hist, n = length(min_pvals_hist)) < 0.05)/length(min_pvals_hist), '\n')
    cat(sum(p.adjust(min_pvals_early, n = length(min_pvals_early)) < 0.05)/length(min_pvals_early), '\n')
    cat(sum(p.adjust(min_pvals_mid, n = length(min_pvals_mid)) < 0.05)/length(min_pvals_mid), '\n')

    cat("Max: \n")
    max_pvals_hist = c()
    max_pvals_early = c()
    max_pvals_mid = c()
    for(species in unique(historic_EU$species)){
        max_pvals_hist = c(max_pvals_hist, ks.test(historic_EU[which(historic_EU$species == species),]$maxPeriodAnnualMeanT, early_EU[which(early_EU$species == species),]$maxPeriodAnnualMeanT)$p.value)
        max_pvals_early = c(max_pvals_early, ks.test(early_EU[which(early_EU$species == species),]$maxPeriodAnnualMeanT, middle_EU[which(middle_EU$species == species),]$maxPeriodAnnualMeanT)$p.value)
        max_pvals_mid = c(max_pvals_mid, ks.test(late_EU[which(late_EU$species == species),]$maxPeriodAnnualMeanT, middle_EU[which(middle_EU$species == species),]$maxPeriodAnnualMeanT)$p.value)
    }
    cat(sum(p.adjust(max_pvals_hist, n = length(max_pvals_hist)) < 0.05)/length(max_pvals_hist), '\n')
    cat(sum(p.adjust(max_pvals_early, n = length(max_pvals_early)) < 0.05)/length(max_pvals_early), '\n')
    cat(sum(p.adjust(max_pvals_mid, n = length(max_pvals_mid)) < 0.05)/length(max_pvals_mid), '\n')

    ### Compute Differences for Distance from Equator
    ###################################################################################################

    stat_dist = function(data, INDEX){
        mu = mean(data[INDEX,]$kmNorthEquator)
        return(mu)
    }

    boots_historic = data.frame(rep(1, 1000))
    boots_early = data.frame(rep(1, 1000))
    boots_middle = data.frame(rep(1, 1000))
    boots_late = data.frame(rep(1, 1000))

    for(species in unique(historic_EU$species)){
        bt_dist = boot(historic_EU[which(historic_EU$species == species), ], stat = stat_dist, R = 1000)
        boots_historic[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(early_EU$species)){
        bt_dist = boot(early_EU[which(early_EU$species == species), ], stat = stat_dist, R = 1000)
        boots_early[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(middle_EU$species)){
        bt_dist = boot(middle_EU[which(middle_EU$species == species), ], stat = stat_dist, R = 1000)
        boots_middle[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(late_EU$species)){
        bt_dist = boot(late_EU[which(late_EU$species == species), ], stat = stat_dist, R = 1000)
        boots_late[species] = bt_dist$t
    }

    cat("Distance from equator: \n")
    cat("Late - Middle vs Middle - Early\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_late[species] - boots_middle[species]), abs(boots_middle[species] - boots_early[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_late[,species] - boots_middle[,species]), abs(boots_middle[,species] - boots_early[,species]), alternative = "greater")$p.value)
    }
    cat(sum(p.adjust(pvals, n = length(pvals)) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p, n = length(wilcox_p)) < 0.05)/length(wilcox_p), '\n')

    cat("Middle - Early vs Early - Historic\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_middle[species] - boots_early[species]), abs(boots_early[species] - boots_historic[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_middle[,species] - boots_early[,species]), abs(boots_early[,species] - boots_historic[,species]), alternative = "greater")$p.value)
    }

    cat(sum(p.adjust(pvals, n = length(pvals)) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p, n = length(wilcox_p)) < 0.05)/length(wilcox_p), '\n')


    ### Compute Differences for Elevation
    ###################################################################################################
    stat_elev = function(data, INDEX){
        mu = mean(data[INDEX,]$elevation)
        return(mu)
    }

    boots_historic = data.frame(rep(1, 1000))
    boots_early = data.frame(rep(1, 1000))
    boots_middle = data.frame(rep(1, 1000))
    boots_late = data.frame(rep(1, 1000))

    for(species in unique(historic_EU$species)){
        bt_dist = boot(historic_EU[which(historic_EU$species == species), ], stat = stat_elev, R = 1000)
        boots_historic[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(early_EU$species)){
        bt_dist = boot(early_EU[which(early_EU$species == species), ], stat = stat_elev, R = 1000)
        boots_early[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(middle_EU$species)){
        bt_dist = boot(middle_EU[which(middle_EU$species == species), ], stat = stat_elev, R = 1000)
        boots_middle[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(late_EU$species)){
        bt_dist = boot(late_EU[which(late_EU$species == species), ], stat = stat_elev, R = 1000)
        boots_late[species] = bt_dist$t
    }

    cat("Elevation: \n")
    cat("Late - Middle vs Middle - Early\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_late[species] - boots_middle[species]), abs(boots_middle[species] - boots_early[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_late[,species] - boots_middle[,species]), abs(boots_middle[,species] - boots_early[,species]), alternative = "greater")$p.value)

    }
    cat(sum(p.adjust(pvals) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p) < 0.05)/length(wilcox_p), '\n')

    cat("Middle - Early vs Early - Historic\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_middle[species] - boots_early[species]), abs(boots_early[species] - boots_historic[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_middle[,species] - boots_early[,species]), abs(boots_early[,species] - boots_historic[,species]), alternative = "greater")$p.value)
    }
    cat(sum(p.adjust(pvals) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p) < 0.05)/length(wilcox_p), '\n')

    ### Compute Differences for Min Annual Mean Temp
    ##########################################################################################################################################################

    stat_min = function(data, INDEX){
        mu = mean(data[INDEX,]$minPeriodAnnualMeanT)
        return(mu)
    }

    boots_historic = data.frame(rep(1, 1000))
    boots_early = data.frame(rep(1, 1000))
    boots_middle = data.frame(rep(1, 1000))
    boots_late = data.frame(rep(1, 1000))

    for(species in unique(historic_EU$species)){
        bt_dist = boot(historic_EU[which(historic_EU$species == species), ], stat = stat_min, R = 1000)
        boots_historic[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(early_EU$species)){
        bt_dist = boot(early_EU[which(early_EU$species == species), ], stat = stat_min, R = 1000)
        boots_early[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(middle_EU$species)){
        bt_dist = boot(middle_EU[which(middle_EU$species == species), ], stat = stat_min, R = 1000)
        boots_middle[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(late_EU$species)){
        bt_dist = boot(late_EU[which(late_EU$species == species), ], stat = stat_min, R = 1000)
        boots_late[species] = bt_dist$t
    }

    cat("Min: \n")
    cat("Late - Middle vs Middle - Early\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_late[species] - boots_middle[species]), abs(boots_middle[species] - boots_early[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_late[,species] - boots_middle[,species]), abs(boots_middle[,species] - boots_early[,species]), alternative = "greater")$p.value)
    }
    cat(sum(p.adjust(pvals) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p) < 0.05)/length(wilcox_p), '\n')

    cat("Middle - Early vs Early - Historic\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_middle[species] - boots_early[species]), abs(boots_early[species] - boots_historic[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_middle[,species] - boots_early[,species]), abs(boots_early[,species] - boots_historic[,species]), alternative = "greater")$p.value)

    }
    cat(sum(p.adjust(pvals) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p) < 0.05)/length(wilcox_p), '\n')

    ### Compute Differences for Max Annual Mean Temp
    ################################################################################################################################################
    stat_max = function(data, INDEX){
        mu = mean(data[INDEX,]$maxPeriodAnnualMeanT)
        return(mu)
    }

    boots_historic = data.frame(rep(1, 1000))
    boots_early = data.frame(rep(1, 1000))
    boots_middle = data.frame(rep(1, 1000))
    boots_late = data.frame(rep(1, 1000))

    for(species in unique(historic_EU$species)){
        bt_dist = boot(historic_EU[which(historic_EU$species == species), ], stat = stat_max, R = 1000)
        boots_historic[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(early_EU$species)){
        bt_dist = boot(early_EU[which(early_EU$species == species), ], stat = stat_max, R = 1000)
        boots_early[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(middle_EU$species)){
        bt_dist = boot(middle_EU[which(middle_EU$species == species), ], stat = stat_max, R = 1000)
        boots_middle[species] = bt_dist$t
    }

    ####################################################################################################
    for(species in unique(late_EU$species)){
        bt_dist = boot(late_EU[which(late_EU$species == species), ], stat = stat_max, R = 1000)
        boots_late[species] = bt_dist$t
    }

    cat("Max: \n")
    cat("Late - Middle vs Middle - Early\n")

    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_late[species] - boots_middle[species]), abs(boots_middle[species] - boots_early[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_late[,species] - boots_middle[,species]), abs(boots_middle[,species] - boots_early[,species]), alternative = "greater")$p.value)
    }
    cat(sum(p.adjust(pvals) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p) < 0.05)/length(wilcox_p), '\n')

    cat("Middle - Early vs Early - Historic\n")
    pvals = c()
    wilcox_p = c()
    for(species in unique(historic_EU$species)){
        pvals = c(pvals, t.test(abs(boots_middle[species] - boots_early[species]), abs(boots_early[species] - boots_historic[species]), alternative = "greater")$p.value)
        wilcox_p = c(wilcox_p, wilcox.test(abs(boots_middle[,species] - boots_early[,species]), abs(boots_early[,species] - boots_historic[,species]), alternative = "greater")$p.value)
    }
    cat(sum(p.adjust(pvals) < 0.05)/length(pvals), '\n')
    cat(sum(p.adjust(wilcox_p) < 0.05)/length(wilcox_p), '\n')
    sink()
