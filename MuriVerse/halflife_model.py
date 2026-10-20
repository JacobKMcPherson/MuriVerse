"""The MuriVerse half-life regressor: an equal-weight ensemble of ridge, RBF SVR and XGBoost.

On a few hundred records from ~100 publications, a single boosted-tree model over-fits the
study it was trained on. Averaging a heavily regularised linear model, a smooth kernel model
and a shallow tree ensemble was the most stable combination under held-out family AND
publication cross-validation, and none of the three has a hyperparameter that is tuned on
the data it is scored on: all are fixed here.

`fit` and `predict` take a raw feature matrix. Missing values are median-imputed for the two
scikit-learn members (the training medians travel with the fitted object) and passed through
untouched to XGBoost, which routes them natively.

A prediction from this model is conditional on the STUDY DESIGN, not just the molecule.
Route of administration reaches it as an ordinary feature column (`route_model`, one-hot over
im/ip/iv/sc in the study_context block), so the number it returns answers "what half-life will
this molecule show by this route", and changing the route column changes the answer.

SPECIES IS A COLUMN. The shipped compact set carries the animal explicitly, as the one-hot
pair `species=mouse` / `species=rat` in the study_context block, so a mouse and a rat dosed
with the same molecule under the same design receive different feature vectors and different
predictions. It also reaches the model a second time through `strain` (BALB/c, C57BL/6 and
ICR/Swiss imply mouse; Sprague-Dawley, Wistar and Wistar Hannover imply rat), but a caller
that leaves strain unknown is no longer asking a species-blind question: the species columns
answer it on their own. Species, like route, is therefore a required input - a caller must
say which animal it means, and a prediction is defined only for an animal the cohort
contains.

Earlier versions of this model had no species column at all. `strain` was the only route by
which the animal reached the matrix, so every row with `strain=unknown` - a large minority of
the cohort, and every molecule the model was applied to from outside it - was species-blind
by construction, and mouse and rat predictions were identical to the last decimal. That is
fixed; see STUDY_CONDITIONS in build_halflife_training_set.py.

Route in particular is not interchangeable. After an intravenous dose the terminal half-life
is an elimination half-life; after an extravascular dose it is one only when absorption is
the faster step, and when it is not the terminal slope reports absorption instead
(flip-flop kinetics). IV and SC labels therefore do not share an absolute scale: a model
fitted on one and applied to the other is offset by roughly half a log unit in the direction
flip-flop predicts. The caller is responsible for supplying the route it actually means, and
for not reading a prediction made for one route as a prediction for another. See
`evaluate_route_flipflop.py` for the measurement and `report_halflife_by_route.py` for
per-route performance, which is very uneven - the pooled score is carried by IV.
"""

from __future__ import annotations

import numpy as np
import xgboost as xgb
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

REFERENCE_DIMS = 80          # ridge penalty is per column: alpha = RIDGE_ALPHA at this width
RIDGE_ALPHA = 10.0
XGB_PARAMS = dict(objective="reg:squarederror", eta=0.05, subsample=0.9,
                  colsample_bytree=0.5, min_child_weight=5, max_depth=2)
XGB_ROUNDS = 150
MEMBERS = ("ridge", "svr", "xgboost")


class HalfLifeEnsemble:
    def __init__(self, seed: int = 0, nthread: int = 1):
        self.seed = seed
        self.nthread = nthread

    def fit(self, X: np.ndarray, y: np.ndarray) -> "HalfLifeEnsemble":
        alpha = RIDGE_ALPHA * max(1.0, X.shape[1] / REFERENCE_DIMS)
        self.ridge_ = make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True),
                                    StandardScaler(), Ridge(alpha=alpha)).fit(X, y)
        self.svr_ = make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True),
                                  StandardScaler(), SVR(C=1.0, epsilon=0.2)).fit(X, y)
        params = dict(XGB_PARAMS, seed=self.seed, nthread=self.nthread)
        self.booster_ = xgb.train(params, xgb.DMatrix(X, label=y), XGB_ROUNDS)
        return self

    def member_predictions(self, X: np.ndarray) -> dict[str, np.ndarray]:
        return {"ridge": self.ridge_.predict(X), "svr": self.svr_.predict(X),
                "xgboost": self.booster_.predict(xgb.DMatrix(X))}

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.mean(list(self.member_predictions(X).values()), axis=0)
