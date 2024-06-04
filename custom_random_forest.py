from sklearn.base import BaseEstimator
import numpy as np
from multiprocessing import Pool
from sklearn.tree import DecisionTreeClassifier

class RandomForestClassifierCustom(BaseEstimator):
    """
    Custom implementation of RandomForestClassifier using multiprocessing.

    Args:
        n_estimators (int, optional): Number of trees in the forest. Defaults to 10.
        max_depth (int, optional): Maximum depth of the tree. Defaults to None.
        max_features (int, optional): Number of features to consider when looking for the best split.
            Defaults to None.
        random_state (int, optional): Controls both the randomness of the bootstrapping of the samples
            used when building trees and the sampling of the features to consider when looking for
            the best split at each node. Defaults to None.

    Attributes:
        classes_ (list): List of classes found during fitting.
        trees (list): List of fitted DecisionTreeClassifier models.
        feat_ids_by_tree (list): List of feature indices used by each tree.
    """
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=None
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def _fit_tree(self, args):
        """
        Fit a single decision tree on a bootstrap sample of the data.

        Args:
            args (tuple): Tuple containing tree_id, feature matrix X, and target vector y.

        Returns:
            tuple: Fitted DecisionTreeClassifier and array of feature indices used by the tree.
        """
        tree_id, X, y = args
        np.random.seed(self.random_state)
        feat_ids = np.random.choice(range(X.shape[1]), size=self.max_features, replace=False)

        sample_indices = np.random.choice(range(len(X)), size=len(X), replace=True)
        X_bootstrap = X[sample_indices][:, feat_ids]
        y_bootstrap = y[sample_indices]

        tree = DecisionTreeClassifier(max_depth=self.max_depth, random_state=self.random_state + tree_id)
        tree.fit(X_bootstrap, y_bootstrap)

        return tree, feat_ids

    def fit(self, X, y, n_jobs=1):
        """
        Fit the random forest classifier on the training data.

        Args:
            X (np.ndarray): Feature matrix of shape (n_samples, n_features).
            y (np.ndarray): Target vector of shape (n_samples,).
            n_jobs (int, optional): Number of jobs to run in parallel. Defaults to 1.

        Returns:
            RandomForestClassifierCustom: Fitted random forest classifier.
        """
        self.classes_ = sorted(np.unique(y))

        with Pool(n_jobs) as pool:
            results = pool.map(self._fit_tree, [(i, X, y) for i in range(self.n_estimators)])

        self.trees, self.feat_ids_by_tree = zip(*results)

        return self

    def _predict_proba_single_tree(self, args):
        """
        Predict class probabilities for a single tree.

        Args:
            args (tuple): Tuple containing a fitted tree, feature indices, and feature matrix X.

        Returns:
            np.ndarray: Predicted class probabilities for the input samples.
        """
        tree, feat_ids, X = args
        X_subset = X[:, feat_ids]
        tree_probas = tree.predict_proba(X_subset)
        return tree_probas

    def predict_proba(self, X, n_jobs=1):
        """
        Predict class probabilities for the input samples.

        Args:
            X (np.ndarray): Feature matrix of shape (n_samples, n_features).
            n_jobs (int, optional): Number of jobs to run in parallel. Defaults to 1.

        Returns:
            np.ndarray: Predicted class probabilities of shape (n_samples, n_classes).
        """
        with Pool(n_jobs) as pool:
            tree_proba_list = pool.map(self._predict_proba_single_tree, [(tree, feat_ids, X) for tree, feat_ids in zip(self.trees, self.feat_ids_by_tree)])

        probas = np.sum(tree_proba_list, axis=0) / self.n_estimators
        return probas

    def predict(self, X, n_jobs=1):
        """
        Predict class labels for the input samples.

        Args:
            X (np.ndarray): Feature matrix of shape (n_samples, n_features).
            n_jobs (int, optional): Number of jobs to run in parallel. Defaults to 1.

        Returns:
            np.ndarray: Predicted class labels of shape (n_samples,).
        """
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions
