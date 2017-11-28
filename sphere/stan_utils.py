#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import pickle
import numpy as np


def save_model(output_path, model):
    with open(output_path, "wb") as f:
        pickle.dump(model, f)


def load_model(pmp):
    with open(pmp, "rb") as f:
        model = pickle.load(f)
    return model


def summarize_fit(fit, pars):
    summary = fit.summary(pars=pars)
    summary_df = pd.DataFrame(summary["summary"],
                              index=summary["summary_rownames"],
                              columns=summary["summary_colnames"])
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


def main():
    pass


if __name__ == '__main__':
    main()
