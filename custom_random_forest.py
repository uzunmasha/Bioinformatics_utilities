from sklearn.base import BaseEstimator
import numpy as np
from multiprocessing import Pool
from sklearn.tree import DecisionTreeClassifier

class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=None
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit_tree(self, args):
        i, X, y = args
        np.random.seed(self.random_state)
        feat_ids = np.random.choice(range(X.shape[1]), size=self.max_features, replace=False)

        sample_indices = np.random.choice(range(len(X)), size=len(X), replace=True)
        X_bootstrap = X[sample_indices][:, feat_ids]
        y_bootstrap = y[sample_indices]

        tree = DecisionTreeClassifier(max_depth=self.max_depth, random_state=self.random_state)
        tree.fit(X_bootstrap, y_bootstrap)

        return tree, feat_ids

    def fit(self, X, y, n_jobs=None):
        self.classes_ = sorted(np.unique(y))

        if n_jobs is None:
            n_jobs = 1

        with Pool(n_jobs) as pool:
            results = pool.map(self.fit_tree, [(i, X, y) for i in range(self.n_estimators)])

        self.trees, self.feat_ids_by_tree = zip(*results)

        return self

    def predict_proba_single_tree(self, args):
        tree, feat_ids, X = args
        X_subset = X[:, feat_ids]
        tree_probas = tree.predict_proba(X_subset)
        return tree_probas

    def predict_proba(self, X, n_jobs=None):
        if n_jobs is None:
            n_jobs = 1

        with Pool(n_jobs) as pool:
            tree_proba_list = pool.map(self.predict_proba_single_tree, [(tree, feat_ids, X) for tree, feat_ids in zip(self.trees, self.feat_ids_by_tree)])

        probas = np.sum(tree_proba_list, axis=0) / self.n_estimators
        return probas

    def predict(self, X, n_jobs=None):
        if n_jobs is None:
            n_jobs = 1

        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions