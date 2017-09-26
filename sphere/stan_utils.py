#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import pickle


def save_model(output_path, model):
    with open(output_path, "wb") as f:
        pickle.dump(model, f)


def load_model(pmp):
    with open(pmp, "rb") as f:
        model = pickle.load(f)
    return model


def summarize_fit(fit):
    summary = fit.summary()
    summary_df = pd.DataFrame(summary["summary"],
                              index=summary["summary_rownames"],
                              columns=summary["summary_colnames"])
    return summary_df


def save_fit(fit, fod):
    with open(fod, "wb") as f:
        pickle.dump(fit, f)


def sampling(model, v_c, si, sw, sc, st, ss):
    stan_data = {
        "I": v_c.size,
        "D": v_c
    }
    fit = model.sampling(data=stan_data,
                         iter=si,
                         warmup=sw,
                         chains=sc,
                         thin=st,
                         seed=ss)
    return fit


def main():
    pass


if __name__ == '__main__':
    main()
