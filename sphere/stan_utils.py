#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import pickle
import numpy as np
from pkg_resources import resource_filename
from scipy.special import logsumexp


def load_model(model_type: str):
    target_path = resource_filename(
        "sphere",
        "stan_models/{0}.pkl".format(model_type)
    )
    with open(target_path, "rb") as f:
        model = pickle.load(f)
    return model


def summarize_fit(fit, pars: list):
    summary = fit.summary(pars=pars)
    summary_df = pd.DataFrame(summary["summary"],
                              index=summary["summary_rownames"],
                              columns=summary["summary_colnames"])
    return summary_df


def summarize_ofit(ofit, pars: list):
    r = []
    if pars is None:
        pars = ofit.keys()
    for k in pars:
        v = ofit[k]
        if hasattr(v.tolist(), "__iter__") is False:
            h = {"": "{0}[0]".format(k), "mle": v}
            r.append(h)
        else:
            for i, vv in enumerate(v):
                if hasattr(vv.tolist(), "__iter__") is False:
                    h = {"": "{0}[{1}]".format(k, i), "mle": vv}
                    r.append(h)
                else:
                    for j, vvv in enumerate(vv):
                        h = {"": "{0}[{1},{2}]".format(k, i, j), "mle": vvv}
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


def get_waic(log_lik: np.ndarray, t: str="bda3") -> float:
    S, n = log_lik.shape
    if t == "bda3":
        # See (Gelman, et al., "BDA3", 2013) Page 174
        # Using logsumexp function to overcome the overflow of log_lik
        lppd = np.sum(-np.log(S) + logsumexp(log_lik, axis=0))
        pwaic = np.sum(np.var(log_lik, axis=0))
        waic = -2.0 * lppd + 2.0 * pwaic
    elif t == "original":
        # See (Sumio Watanabe, 2010, JMLR) formula (4), (5), (6)
        T = - np.mean(-np.log(S) + logsumexp(log_lik, axis=0))
        fV = np.mean(np.var(log_lik, axis=0))
        waic = T + fV
    return waic


def sampling(model, stan_data: dict, pars: list, si, sw, sc, st, ss, n_jobs):
    fit = model.sampling(data=stan_data,
                         pars=pars,
                         iter=si,
                         warmup=sw,
                         chains=sc,
                         thin=st,
                         seed=ss,
                         n_jobs=n_jobs)
    return fit


def optimizing(model, stan_data: dict, ss: int, om: str=None, sh: int=5):
    # init_alpha must be lower to estimate non-normalizing model correctly
    fit = model.optimizing(data=stan_data,
                           init_alpha=1e-10,
                           iter=1e4,
                           refresh=1,
                           algorithm=om,
                           seed=ss)

    return fit


def main():
    pass


if __name__ == '__main__':
    main()
