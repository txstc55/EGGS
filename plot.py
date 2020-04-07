import json
import matplotlib.pyplot as plt
import numpy as np
import sys
from math import ceil, log

json_file_name = "result_datas/all_result.json"

results = []
if len(sys.argv) == 2:
    with open(sys.argv[1], 'r') as j:
        results = json.load(j)
        json_file_name = sys.argv[1]
else:
    with open('result_datas/all_result.json', 'r') as j:
        results = json.load(j)

assemble = [float(x["ASSEMBLE"])/100.0/1000.0 for x in results]
compute = [float(x["COMPUTE"])/100.0/1000.0 for x in results]
solve = [float(x["SOLVE"])/100.0/1000.0 for x in results]

sum = []
for i in range(len(assemble)):
    sum.append(assemble[i]+compute[i]+solve[i])

# total_time =[assemble[i]+compute[i]+solve[i] for i in range(5)]

bars = np.add(solve, compute).tolist()

N = len(results)

ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, solve, width, color="#34495e")
p2 = plt.bar(ind, compute, width, bottom=solve, color="#2ecc71")
p3 = plt.bar(ind, assemble, width, bottom=bars, color="#3498db")


best_log_10 = ceil(log(max(sum), 10))
div = max(sum)/(10**(best_log_10-1))+1
max_y = div*(10**(best_log_10-1))
gap = 10**(best_log_10-1)

plt.ylabel('Time in Milliseconds')
# plt.title('Time Result For Each Method')
# plt.xlabel('Method')
plt.xticks(ind, [x["name"] for x in results])
plt.yticks(np.arange(0, max_y, gap))
plt.legend((p1[0], p2[0], p3[0]), ('Solve', 'Compute', 'Assemble'))


for r1, r2, r3 in zip(p1, p2, p3):
    h1 = r1.get_height()
    h2 = r2.get_height()
    h3 = r3.get_height()
    if (h1!=0.0):
        plt.text(r1.get_x() + r1.get_width() / 2., h1 / 2., "%.3f" % h1,
                ha="center", va="center", color="white", fontsize=5, fontweight="bold")
    plt.text(r2.get_x() + r2.get_width() / 2., h1 + h2 / 2., "%.3f" % h2,
             ha="center", va="center", color="white", fontsize=5, fontweight="bold")
    if (h3!=0.0):
        plt.text(r3.get_x() + r3.get_width() / 2., h1 + h2 + h3 / 2., "%.3f" %
                h3, ha="center", va="center", color="white", fontsize=5, fontweight="bold")
save_plot_name = json_file_name.split(".")[0].split("/")[-1]


# plt.gca().set_axis_off()
plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
                    hspace=0, wspace=0)
# plt.margins(0, 0)
# plt.gca().xaxis.set_major_locator(plt.NullLocator())
# plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.savefig("test_result_graphs/"+save_plot_name+".pdf", bbox_inches='tight',
            pad_inches=0, dpi=200)


# plt.savefig("test_result_graphs/"+save_plot_name+".png", dpi=200)


print(results)
