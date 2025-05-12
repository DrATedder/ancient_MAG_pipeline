from scipy import stats

def create_data_dict(data_file):
    data_dict = {}
    with open(data_file, "r") as f_in:
        for lines in f_in:
            if lines.startswith("run") and "deamSim" not in lines:
                modern = int(lines.split(":")[1].strip())
            elif "deamSim" in lines:
                ancient = int(lines.split(":")[1].strip())
                data_dict[lines.split("_")[0]] = [modern, ancient]
    return data_dict

def paired_ttest(list_1, list_2):
    t_statistic, p_value = stats.ttest_rel(list_1, list_2)
    return t_statistic, p_value

alpha = 0.05
horrible_data = "/home/andrew/Documents/FS_MAG_Pipeline/MAGsMatched_results.txt"

ancient = []
modern = []

for x, y in create_data_dict(horrible_data).items():
    ancient.append(y[1])
    modern.append(y[0])


t_statistic, p_value = paired_ttest(modern, ancient)

if p_value < alpha:
    print("P_value is {0}, the groups are statistically different!".format(p_value))
else:
    print("P_value is {0}, the groups are not statistically different.".format(p_value))
