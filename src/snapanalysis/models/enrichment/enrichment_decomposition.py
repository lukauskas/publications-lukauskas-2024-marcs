import numpy as np
from matplotlib import pyplot as plt
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA

class EnrichmentDecoposition(BaseEstimator, TransformerMixin):
    def __init__(self, random_state=None,
                 std_threshold=2.0):
        self.random_state = random_state
        self.std_threshold = std_threshold

        super(EnrichmentDecoposition, self).__init__()

    @classmethod
    def vec_to_angle(cls, vec):
        angle = np.arctan2(vec[1], vec[0])

        half_pi = np.pi / 2.0

        if angle > half_pi:
            angle -= np.pi
        if angle < -half_pi:
            angle += np.pi

        return angle

    def inliers(self, X):
        X = np.asarray(X)

        assert X.shape[1] == 2

        sign = np.sign(X)
        inliers = sign[:, 0] != sign[:, 1]

        return inliers

    def fit(self, X, y=None):
        X = np.asarray(X)

        assert X.shape[1] == 2

        inliers = self.inliers(X)
        subdata = X[inliers]

        # First PCA
        pca = PCA(n_components=1, random_state=self.random_state)
        pca.fit(subdata)
        vec = pca.components_[0]
        self.angle_ = -self.vec_to_angle(vec)
        theta = self.angle_

        # Find outliers on residual dimension
        rotated = self.rotate_data(subdata, theta=theta)
        median = np.median(rotated[:, 1])
        std = np.std(rotated[:, 1], ddof=1)
        k = self.std_threshold
        outliers = rotated[:, 1] > median + k * std
        outliers |= rotated[:, 1] < median - k * std

        # Now that we have the new outliers, do another PCA
        reiter_subdata = subdata[~outliers]
        pca = PCA(n_components=1, random_state=self.random_state)
        pca.fit(reiter_subdata)
        vec = pca.components_[0]
        theta = -self.vec_to_angle(vec)

        self.angle_ = theta

        return self

    @classmethod
    def plot_angled_line(cls, angle, length=30, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()

        angle = -angle

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        x = np.cos(angle) * length
        y = np.sin(angle) * length

        ax.plot([-x, x], [-y, y], **kwargs)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    def plot_axis(self, **kwargs):
        self.plot_angled_line(self.angle_, **kwargs)

    @classmethod
    def rotate_data(cls, data, theta):
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        rotation = np.matrix([[cos_theta, -sin_theta],
                              [sin_theta, cos_theta]])

        rotated = np.matmul(rotation, data.T).T

        return np.asarray(rotated)

    def transform(self, X):
        theta = self.angle_
        return self.rotate_data(X, theta)
