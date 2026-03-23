import sys
import os
import numpy as np
import numpy.random as npr
import scipy
from scipy.stats import beta, norm
import pandas as pd
import math
import re
import timeit

# -----------------------------
# True density functions
# -----------------------------
def trig_density(testdata):
    testdata = np.asarray(testdata)
    n_test = testdata.shape[0]
    trueD = np.zeros(n_test)
    for i in range(n_test):
        trueD[i] = 1 + math.sin(2 * math.pi * testdata[i] - math.pi / 2)
    return trueD

def beta_density(testdata):
    testdata = np.asarray(testdata)
    n_test = testdata.shape[0]
    trueD = np.zeros(n_test)
    for i in range(n_test):
        trueD[i] = 0.3 * beta.pdf(testdata[i], 300, 10) + 0.7 * beta.pdf(testdata[i], 3, 12)
    return trueD

def circle_density(testdata):
    testdata = np.asarray(testdata)
    n_test = testdata.shape[0]
    trueD = np.zeros(n_test)
    for i in range(n_test):
        if testdata[i] >0.5 and testdata[i] <1:
            trueD[i] = (8.0/ math.pi) * np.sqrt(1 - (4* testdata[i] - 3) * (4*testdata[i] -3))
    return trueD

# The example is not used
def normAR1_density(testdata, mu, rho, sigma):
    testdata = np.asarray(testdata)

    if testdata.ndim == 2:
        n_test = testdata.shape[0]
        dim = testdata.shape[1]
    else:
        n_test = 1
        dim = testdata.shape[0]
        testdata = np.reshape(testdata, (-1, dim))

    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = sigma * sigma * pow(rho, abs(i - j))

    det = np.linalg.det(A)
    Theta = np.linalg.inv(A)

    trueD = np.zeros(n_test)
    for i in range(n_test):
        diff = testdata[i, :] - mu
        temp = np.matmul(diff, Theta)
        B = np.matmul(temp, diff)
        trueD[i] = math.exp(-0.5 * B) / math.sqrt(pow((2 * math.pi), dim) * det)

    return trueD

def mix5Dnorm_density(testdata):
    testdata = np.asarray(testdata)

    if testdata.ndim == 2:
        n_test = testdata.shape[0]
        dim = testdata.shape[1]
    else:
        n_test = 1
        dim = testdata.shape[0]
        testdata = np.reshape(testdata, (-1, dim))

    trueD = np.zeros(n_test)
    mu = np.array([0.25, 0.25])

    for i in range(n_test):
        y = testdata[i, 0:2]
        d1 = normAR1_density(y, mu, 0.36, 0.05)
        d1 *= norm.pdf(testdata[i, 2], 0.25, 0.05)

        d2 = 1.0
        for j in range(3):
            d2 *= norm.pdf(testdata[i, j], 0.75, 0.05)

        trueD[i] = 0.5 * d1 + 0.5 * d2

        j = 3
        while j < 5:
            trueD[i] *= norm.pdf(testdata[i, j], 0.5, 0.1)
            j += 1

        while j < dim:
            trueD[i] *= (0.5 * norm.pdf(testdata[i, j], 0.35, 0.1) +
                         0.5 * norm.pdf(testdata[i, j], 0.6, 0.05))
            j += 1

    return trueD
    

def test_trueD(X_test, dist):
    """
    dist can be 'trig', 'beta', 'circle', 'mix5Dnorm', or 'normAR1_Y' where Y is a digit.
    Example: 'normAR1_6' -> rho=0.6
    """
    match = re.search(r"normAR1_(\d)", dist)

    if dist == "trig":
        true_density = trig_density(X_test)
    elif dist == "beta":
        true_density = beta_density(X_test)
    elif dist == "circle":
        true_density = circle_density(X_test)
    elif dist == "mix5Dnorm":
        true_density = mix5Dnorm_density(X_test)
    elif match:
        rho = float(match.group(1)) / 10.0
        print(rho)
        # TODO: define mu and sigma appropriately
        # mu =
        # sigma =
        true_density = normAR1_density(X_test, mu, rho, sigma)  # noqa: F821
    else:
        print("Please check the distribution", file=sys.stderr)
        sys.exit(1)

    print(true_density[0:10])
    return true_density


def bd_selection(X_train, dim, K=2, para_min=0.1, para_max=5, para_len=49):
    para_list = np.linspace(para_min, para_max, para_len)
    test_lik = np.zeros(para_len)

    n = X_train.shape[0]
    npr.seed(124)
    rindex = npr.choice(n, size=n, replace=False)

    for i in range(para_len):
        para = para_list[i]
        temp_lik = np.zeros(K)

        for k in range(K):
            train_index = np.zeros(n, dtype=bool)
            train_loc = rindex[int(k * (n / K)): int((k + 1) * (n / K))]
            train_index[train_loc] = True
            test_index = ~train_index  # boolean complement

            if dim == 1:
                temp_train = X_train[train_index]
                temp_test = X_train[test_index]
                kernel_est = scipy.stats.gaussian_kde(temp_train, bw_method=n ** (-para / (dim + 4)))
                est_density = kernel_est(temp_test)
            else:
                temp_train = X_train[train_index, :]
                temp_test = X_train[test_index, :]
                kernel_est = scipy.stats.gaussian_kde(temp_train.T, bw_method=n ** (-para / (dim + 4)))
                est_density = kernel_est(temp_test.T)

            count = 0
            for value in est_density:
                if value > 0:
                    temp_lik[k] += -math.log(value)
                    count += 1

            if count > test_index.sum() * 0.95:
                temp_lik[k] /= count
            else:
                temp_lik[k] = 1e5

        test_lik[i] = temp_lik.mean()

    print(test_lik)
    min_index = test_lik.argmin()
    print("The selected bandwidth is:")
    print(para_list[min_index] / (4 + dim))
    return para_list[min_index]
    
def _scott_factor(n, d):
    return n ** (-1.0 / (d + 4.0))

def fit_gaussian_kde_explicit_fast(X_train, bw_factor=None, cov=None, chunk_size=1000, jitter=1e-10):
    """
    Fast explicit Gaussian KDE with full bandwidth covariance matrix.

    KDE:
        f(y) = (1/n) * sum_i N(y | x_i, H)
    where H = bw_factor^2 * Cov(X_train) (or user-provided cov)

    Parameters
    ----------
    X_train : array-like, shape (n, d) or (n,)
    bw_factor : float or None
        If None, Scott's factor n^{-1/(d+4)}.
    cov : array-like or None
        Optional covariance matrix (d,d). If None, uses sample covariance of X_train.
    chunk_size : int
        Number of test points to process per chunk (memory/time tradeoff).
    jitter : float
        Diagonal jitter added to H for numerical stability.

    Returns
    -------
    kde : callable
        kde(X_eval) -> densities (m,)
    """
    X = np.asarray(X_train, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    n, d = X.shape

    # Covariance of training data (unbiased ddof=1)
    if cov is None:
        C = np.cov(X, rowvar=False, ddof=1)
    else:
        C = np.asarray(cov, dtype=float)
    if d == 1:
        C = np.atleast_2d(C)

    if bw_factor is None:
        bw_factor = _scott_factor(n, d)

    H = (bw_factor ** 2) * C

    # Stabilize and Cholesky
    H = H + jitter * np.eye(d)
    try:
        L = np.linalg.cholesky(H)  # H = L L^T
    except np.linalg.LinAlgError:
        # Increase jitter if needed
        H = H + (100 * jitter) * np.eye(d)
        L = np.linalg.cholesky(H)

    # Whiten: solve L * z = x  => z = L^{-1} x
    # Do it for all rows efficiently by solving for transposed matrices.
    Xw = np.linalg.solve(L, X.T).T  # (n, d)

    # Precompute norms for distance expansion
    Xw_norm2 = np.sum(Xw * Xw, axis=1)  # (n,)

    detH = float(np.linalg.det(H))
    norm_const = 1.0 / (n * math.sqrt(((2.0 * math.pi) ** d) * detH))

    def kde(X_eval):
        Y = np.asarray(X_eval, dtype=float)
        if Y.ndim == 1:
            Y = Y[:, None]
        m, d2 = Y.shape
        if d2 != d:
            raise ValueError(f"Dimension mismatch: train dim={d}, eval dim={d2}")

        # Whiten test points
        Yw = np.linalg.solve(L, Y.T).T  # (m, d)
        Yw_norm2 = np.sum(Yw * Yw, axis=1)  # (m,)

        out = np.empty(m, dtype=float)

        # Chunk over Y to control memory: each chunk builds (b, n) distances
        for start in range(0, m, chunk_size):
            end = min(start + chunk_size, m)
            Yc = Yw[start:end]                 # (b, d)
            Yc_norm2 = Yw_norm2[start:end]     # (b,)

            # dist2 = ||y||^2 + ||x||^2 - 2 y·x
            # shapes: (b,1) + (1,n) - 2(b,d)(d,n) -> (b,n)
            cross = Yc @ Xw.T                  # (b, n)
            dist2 = Yc_norm2[:, None] + Xw_norm2[None, :] - 2.0 * cross

            # Sum kernels
            Ksum = np.exp(-0.5 * dist2).sum(axis=1)  # (b,)
            out[start:end] = norm_const * Ksum

        return out

    return kde


def kdeComparison(dist, trainDir, testDir, outDir, rep):
    test_file = os.path.join(testDir, "testdata.txt")

    data_test = pd.read_csv(test_file, header=None, sep="\t")
    dim = data_test.shape[1] - 1
    n_test = data_test.shape[0]

    print("Test data dimension is:\n")
    print(dim, n_test)

    X_test = np.array(data_test.values[:, 0:dim])
    if dim == 1:
        X_test = X_test.reshape(-1)

    print("The first few lines of test data are:\n")
    print(data_test.head())

    true_density = test_trueD(X_test, dist)

    os.makedirs(outDir, exist_ok=True)
    outfile = open(os.path.join(outDir, "kde_KL.txt"), "w")
    timefile = open(os.path.join(outDir, "kde_time.txt"), "w")

    KL_kernel = np.zeros(rep)

    for iter in range(rep):
        train_file = os.path.join(trainDir, f"rep_{iter}", f"data.txt")
        data_train = pd.read_csv(train_file, header=None, sep="\t")

        n = data_train.shape[0]
        print("Training sample size is:")
        print(n)

        X_train = np.array(data_train.values[:, 0:dim])
        if dim == 1:
            X_train = X_train.reshape(-1)
            start = timeit.default_timer()
            bd_para = bd_selection(X_train, dim, K=2, para_min=0.1, para_max=5, para_len=49)
            kernel_est = scipy.stats.gaussian_kde(X_train, bw_method=(n) ** (-bd_para / (dim + 4)))
            stop = timeit.default_timer()
            timefile.write(str(stop - start) + "\n")
            est_density = kernel_est(X_test)
        else:
            start = timeit.default_timer()
            # bd_para = bd_selection(X_train, dim)
            kernel_est = fit_gaussian_kde_explicit_fast(X_train, bw_factor=None, chunk_size=1000)
            #kernel_est = scipy.stats.gaussian_kde(X_train.T)  # , bw_method=n**(-bd_para/(dim+4)))
            stop = timeit.default_timer()
            timefile.write(str(stop - start) + "\n")
            est_density = kernel_est(X_test)
            #est_density = kernel_est(X_test.T)

        KL = 0.0
        k_count = 0

        for i in range(n_test):
            if (true_density[i] > 0 and est_density[i] > 0 and
                    not np.isnan(math.log(true_density[i] / est_density[i]))):
                KL += math.log(true_density[i] / est_density[i])
                k_count += 1
            elif est_density[i] == 0 and true_density[i] > 0:
                k_count += 1
            else:
                print(est_density[i], true_density[i])

            if i < 10:
                print(est_density[i], true_density[i])

        print(k_count)
        KL_kernel[iter] = KL / k_count
        outfile.write(str(KL_kernel[iter]) + "\n")

    outfile.close()
    timefile.close()
    return KL_kernel


def main(argv):
    if len(argv) < 6:
        print("Usage: python comparison.py <dist> <train_dir> <test_dir> <out_dir> <n_rep>",
              file=sys.stderr)
        sys.exit(1)

    dist = argv[1]
    trainDir = argv[2]
    testDir = argv[3]
    outDir = argv[4]
    rep = int(argv[5])

    kdeComparison(dist, trainDir, testDir, outDir, rep)


if __name__ == "__main__":
    main(sys.argv)
