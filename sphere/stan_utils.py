#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import pickle
import numpy as np
from pkg_resources import resource_filename


def load_model(model_type):
    target_path = resource_filename(
        "sphere",
        "stan_models/{0}.pkl".format(model_type)
    )
    with open(target_path, "rb") as f:
        model = pickle.load(f)
    return model


def summarize_fit(fit, pars):
    summary = fit.summary(pars=pars)
    summary_df = pd.DataFrame(summary["summary"],
                              index=summary["summary_rownames"],
                              columns=summary["summary_colnames"])
    return summary_df


def summarize_ofit(ofit, pars):
    r = []
    if pars is None:
        pars = ofit.keys()
    for k in ofit.keys():
        if k in pars:
            a = ofit[k]
            if a.size == 1:
                h = {"": k, "mle": a}
                r.append(h)
            else:
                for i, v in enumerate(a):
                    h = {"": "{0}[{1}]".format(k, i), "mle": v}
                    r.append(h)
    summary_df = pd.DataFrame(r)
    summary_df = summary_df.set_index("")
    return summary_df


def save_fit(fit, fod: str):
    with open(fod, "wb") as f:
        pickle.dump(fit, f)


def save_log_lik(fit, lld: str):
    log_lik = fit.extract("log_lik")["log_lik"]
    np.savetxt(lld, log_lik, delimiter="\t")


def load_log_lik(llp: str) -> np.ndarray:
    log_lik = np.loadtxt(llp, delimiter="\t")
    return log_lik


def get_waic(log_lik: np.ndarray, t="bda3") -> float:
    if t == "bda3":
        lppd = np.sum(np.log(np.mean(np.exp(log_lik), axis=0)))
        pwaic = np.sum(np.var(log_lik, axis=0))
        waic = -2.0 * lppd + 2.0 * pwaic
    elif t == "original":
        T = - np.mean(np.log(np.mean(np.exp(log_lik), axis=0)))
        V_div_N = np.mean(np.mean(np.power(log_lik, 2), axis=0)
                          - np.power(np.mean(log_lik, axis=0), 2))
        waic = T + V_div_N
    return waic


def sampling(model, stan_data: dict, pars: list, si, sw, sc, st, ss):
    fit = model.sampling(data=stan_data,
                         pars=pars,
                         iter=si,
                         warmup=sw,
                         chains=sc,
                         thin=st,
                         seed=ss,
                         n_jobs=sc)
    return fit


def optimizing(model, stan_data: dict, ss):
    fit = model.optimizing(data=stan_data,
                           init_alpha=1e-6,
                           tol_obj=1e-4,
                           seed=ss)
    return fit


def main():
    pass


if __name__ == '__main__':
    main()
