import scipy.stats as stats


def set_IIV_vacc(eff_lower=0.75, eff_upper=0.9, eff_SD=0.025):
    cur_eff = 0.8
    new_eff = stats.truncnorm.rvs((eff_lower - cur_eff) / eff_SD, (eff_upper - cur_eff) / eff_SD, loc=cur_eff,
                                  scale=eff_SD)
    return new_eff


new_eff = []
for index, val in enumerate(range(1000)):
    new_eff.append(set_IIV_vacc(eff_lower=0.75, eff_upper=0.9, eff_SD=0.025))
