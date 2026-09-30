import numpy as np
from scipy import stats


# ============================================================
# SETTINGS
# ============================================================

DATA_FILE = "/home/andrew/Documents/FS_MAG_Pipeline/MAGsMatched_results.txt"

ALPHA = 0.05
N_BOOTSTRAP = 100_000
N_PERMUTATIONS = 100_000

# Reproducible random-number generator
RNG = np.random.default_rng(42)


# ============================================================
# READ DATA
# ============================================================

def create_data_dict(data_file):
    """
    Read the accumulated results file.

    Returns:
        {
            "run1": [modern_count, deamSim_count],
            "run2": [modern_count, deamSim_count],
            ...
        }

    modern_count = false positives in the original data
    deamSim_count = false positives after deamination simulation
    """

    data_dict = {}

    with open(data_file, "r") as f_in:

        for line in f_in:

            line = line.strip()

            if not line:
                continue

            if not line.startswith("run"):
                continue

            filename, value = line.split(":")
            value = int(value)

            run_name = filename.split("_")[0]

            if run_name not in data_dict:
                data_dict[run_name] = [None, None]

            if "deamSim" in filename:
                data_dict[run_name][1] = value
            else:
                data_dict[run_name][0] = value

    return data_dict


# ============================================================
# EXTRACT PAIRED COUNTS
# ============================================================

def extract_paired_counts(data_dict):

    modern = []
    deamSim = []
    runs = []

    for run, values in sorted(data_dict.items()):

        modern_count = values[0]
        deam_count = values[1]

        if modern_count is None or deam_count is None:
            print(f"WARNING: incomplete pair for {run}")
            continue

        runs.append(run)
        modern.append(modern_count)
        deamSim.append(deam_count)

    return (
        np.asarray(modern, dtype=float),
        np.asarray(deamSim, dtype=float),
        runs
    )


# ============================================================
# DESCRIPTIVE STATISTICS
# ============================================================

def describe_data(modern, deamSim):

    differences = deamSim - modern

    print("\n" + "=" * 70)
    print("DESCRIPTIVE STATISTICS")
    print("=" * 70)

    print(f"Number of paired runs: {len(modern)}")

    print("\nModern:")
    print(f"  Mean:   {np.mean(modern):.3f}")
    print(f"  Median: {np.median(modern):.3f}")
    print(f"  SD:     {np.std(modern, ddof=1):.3f}")
    print(f"  Min:    {np.min(modern):.0f}")
    print(f"  Max:    {np.max(modern):.0f}")

    print("\nDeamSim:")
    print(f"  Mean:   {np.mean(deamSim):.3f}")
    print(f"  Median: {np.median(deamSim):.3f}")
    print(f"  SD:     {np.std(deamSim, ddof=1):.3f}")
    print(f"  Min:    {np.min(deamSim):.0f}")
    print(f"  Max:    {np.max(deamSim):.0f}")

    print("\nPaired difference (DeamSim - Modern):")
    print(f"  Mean:   {np.mean(differences):.3f}")
    print(f"  Median: {np.median(differences):.3f}")
    print(f"  SD:     {np.std(differences, ddof=1):.3f}")
    print(f"  Min:    {np.min(differences):.0f}")
    print(f"  Max:    {np.max(differences):.0f}")

    print("\nDirection of paired differences:")
    print(f"  DeamSim > Modern: {np.sum(differences > 0)}")
    print(f"  DeamSim = Modern: {np.sum(differences == 0)}")
    print(f"  DeamSim < Modern: {np.sum(differences < 0)}")


# ============================================================
# WILCOXON SIGNED-RANK TEST
# ============================================================

def wilcoxon_test(modern, deamSim):

    result = stats.wilcoxon(
        modern,
        deamSim,
        alternative="two-sided",
        zero_method="wilcox",
        method="auto"
    )

    print("\n" + "=" * 70)
    print("WILCOXON SIGNED-RANK TEST")
    print("=" * 70)

    print(f"Statistic: {result.statistic}")
    print(f"P-value:  {result.pvalue:.6g}")

    if result.pvalue < ALPHA:
        print("Result: statistically significant at alpha =", ALPHA)
    else:
        print("Result: not statistically significant at alpha =", ALPHA)

    return result


# ============================================================
# PAIRED BOOTSTRAP
# ============================================================

def paired_bootstrap(modern, deamSim, n_bootstrap=100_000):

    """
    Paired bootstrap.

    IMPORTANT:
    We resample RUNS, not individual observations independently.

    This preserves the pairing between modern and deamSim.

    Returns bootstrap distributions for:
        mean difference
        median difference
    """

    n = len(modern)

    mean_differences = np.empty(n_bootstrap)
    median_differences = np.empty(n_bootstrap)

    for i in range(n_bootstrap):

        indices = RNG.integers(
            0,
            n,
            size=n
        )

        modern_sample = modern[indices]
        deam_sample = deamSim[indices]

        differences = deam_sample - modern_sample

        mean_differences[i] = np.mean(differences)
        median_differences[i] = np.median(differences)

    return mean_differences, median_differences


def bootstrap_ci(values, confidence=0.95):

    lower = (1 - confidence) / 2
    upper = 1 - lower

    return (
        np.quantile(values, lower),
        np.quantile(values, upper)
    )


def report_bootstrap(modern, deamSim):

    print("\n" + "=" * 70)
    print("PAIRED BOOTSTRAP")
    print("=" * 70)

    mean_boot, median_boot = paired_bootstrap(
        modern,
        deamSim,
        N_BOOTSTRAP
    )

    observed_mean = np.mean(deamSim - modern)
    observed_median = np.median(deamSim - modern)

    mean_ci = bootstrap_ci(mean_boot)
    median_ci = bootstrap_ci(median_boot)

    print("\nMean paired difference:")
    print(f"  Observed: {observed_mean:.3f}")
    print(
        f"  95% bootstrap CI: "
        f"[{mean_ci[0]:.3f}, {mean_ci[1]:.3f}]"
    )

    print("\nMedian paired difference:")
    print(f"  Observed: {observed_median:.3f}")
    print(
        f"  95% bootstrap CI: "
        f"[{median_ci[0]:.3f}, {median_ci[1]:.3f}]"
    )

    return {
        "mean_bootstrap": mean_boot,
        "median_bootstrap": median_boot,
        "mean_ci": mean_ci,
        "median_ci": median_ci
    }


# ============================================================
# PAIRED PERMUTATION TEST
# ============================================================

def paired_permutation_test(modern, deamSim, n_permutations=100_000):

    """
    Paired randomization test.

    Under the null hypothesis that the two conditions are
    exchangeable within each run, randomly swap modern/deamSim
    within each pair.

    The test statistic is the mean paired difference.
    """

    observed_difference = np.mean(deamSim - modern)

    differences = deamSim - modern

    n = len(differences)

    # Generate random +/- signs for each paired difference.
    signs = RNG.choice(
        [-1, 1],
        size=(n_permutations, n)
    )

    permuted_means = np.mean(
        signs * differences,
        axis=1
    )

    # Two-sided p-value
    p_value = (
        np.sum(
            np.abs(permuted_means)
            >= abs(observed_difference)
        ) + 1
    ) / (n_permutations + 1)

    print("\n" + "=" * 70)
    print("PAIRED PERMUTATION TEST")
    print("=" * 70)

    print(f"Observed mean difference: {observed_difference:.3f}")
    print(f"P-value: {p_value:.6g}")

    if p_value < ALPHA:
        print("Result: statistically significant at alpha =", ALPHA)
    else:
        print("Result: not statistically significant at alpha =", ALPHA)

    return observed_difference, p_value


# ============================================================
# SIGN TEST
# ============================================================

def sign_test(modern, deamSim):

    """
    Paired sign test.

    Ignores the magnitude of the differences and considers only
    whether deamSim is greater or smaller than modern.

    Ties are excluded.
    """

    differences = deamSim - modern

    positive = np.sum(differences > 0)
    negative = np.sum(differences < 0)
    zero = np.sum(differences == 0)

    n_nonzero = positive + negative

    if n_nonzero == 0:
        p_value = 1.0
    else:
        result = stats.binomtest(
            positive,
            n_nonzero,
            p=0.5,
            alternative="two-sided"
        )

        p_value = result.pvalue

    print("\n" + "=" * 70)
    print("PAIRED SIGN TEST")
    print("=" * 70)

    print(f"Positive differences: {positive}")
    print(f"Negative differences: {negative}")
    print(f"Zero differences:     {zero}")

    print(f"P-value: {p_value:.6g}")

    if p_value < ALPHA:
        print("Result: statistically significant at alpha =", ALPHA)
    else:
        print("Result: not statistically significant at alpha =", ALPHA)

    return p_value


# ============================================================
# SENSITIVITY ANALYSIS
# ============================================================

def sensitivity_analysis(modern, deamSim, wilcoxon_result,
                         permutation_result, sign_result):

    print("\n" + "=" * 70)
    print("SENSITIVITY ANALYSIS")
    print("=" * 70)

    print()
    print("Method                         P-value")
    print("-" * 50)

    print(
        f"Wilcoxon signed-rank       "
        f"{wilcoxon_result.pvalue:.6g}"
    )

    print(
        f"Paired permutation         "
        f"{permutation_result[1]:.6g}"
    )

    print(
        f"Paired sign test           "
        f"{sign_result:.6g}"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    # --------------------------------------------------------
    # Read data
    # --------------------------------------------------------

    data_dict = create_data_dict(DATA_FILE)

    modern, deamSim, runs = extract_paired_counts(data_dict)

    # --------------------------------------------------------
    # Display paired data
    # --------------------------------------------------------

    print("=" * 70)
    print("PAIRED DATA")
    print("=" * 70)

    for run, m, d in zip(runs, modern, deamSim):

        print(
            f"{run:>6}   "
            f"Modern = {int(m):>3}   "
            f"DeamSim = {int(d):>3}   "
            f"Difference = {int(d - m):>4}"
        )

    # --------------------------------------------------------
    # Descriptive statistics
    # --------------------------------------------------------

    describe_data(modern, deamSim)

    # --------------------------------------------------------
    # Wilcoxon
    # --------------------------------------------------------

    wilcoxon_result = wilcoxon_test(
        modern,
        deamSim
    )

    # --------------------------------------------------------
    # Bootstrap
    # --------------------------------------------------------

    bootstrap_results = report_bootstrap(
        modern,
        deamSim
    )

    # --------------------------------------------------------
    # Paired permutation test
    # --------------------------------------------------------

    permutation_result = paired_permutation_test(
        modern,
        deamSim,
        N_PERMUTATIONS
    )

    # --------------------------------------------------------
    # Sign test
    # --------------------------------------------------------

    sign_result = sign_test(
        modern,
        deamSim
    )

    # --------------------------------------------------------
    # Sensitivity analysis
    # --------------------------------------------------------

    sensitivity_analysis(
        modern,
        deamSim,
        wilcoxon_result,
        permutation_result,
        sign_result
    )

    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    observed_difference = np.mean(deamSim - modern)

    print(
        f"Mean difference (DeamSim - Modern): "
        f"{observed_difference:.3f}"
    )

    print(
        f"Median difference (DeamSim - Modern): "
        f"{np.median(deamSim - modern):.3f}"
    )

    print(
        f"Wilcoxon p-value: "
        f"{wilcoxon_result.pvalue:.6g}"
    )

    print(
        f"Permutation p-value: "
        f"{permutation_result[1]:.6g}"
    )

    print(
        f"Sign-test p-value: "
        f"{sign_result:.6g}"
    )

    print(
        f"Bootstrap 95% CI for mean difference: "
        f"[{bootstrap_results['mean_ci'][0]:.3f}, "
        f"{bootstrap_results['mean_ci'][1]:.3f}]"
    )

    print(
        f"Bootstrap 95% CI for median difference: "
        f"[{bootstrap_results['median_ci'][0]:.3f}, "
        f"{bootstrap_results['median_ci'][1]:.3f}]"
    )


if __name__ == "__main__":
    main()

