import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
file_list = ["pipe", "sypr", "syrk"]
file_list_15 = [x+"_15" for x in file_list]
file_list += file_list_15
file_list_full = ["result_datas/"+x+".json" for x in file_list]

x_range = [2, 3, 4, 5, 6]
x_range = [10**x for x in x_range]
for k in range(len(file_list_full)):
    f = file_list_full[k]
    with open(f, 'r') as j:
        results = json.load(j)
        print(results)
        time_result = {}
        pre_comp_time = []
        for r in results:
            time = r['time']
            time = [float(x) for x in time]
            time_result[r['method_name']] = time
            if (r['method_name'] == "OURS MULTI THREAD"):
                pre_comp_time = r["pre_comp"]

        pre_comp_time = [float(x) for x in pre_comp_time]

        multi_thread_factor = []
        single_thread_factor = []
        eigen_single_thread_factor = []
        pre_comp_time = []
        for i in range(5):
            multi_thread_factor.append(
                time_result["MKL MULTI THREAD"][i]/time_result["OURS MULTI THREAD"][i])
        plt.plot(x_range, multi_thread_factor,
                 label="Four Threads Speedup", color="#1B9CFC")
        if ("MKL SINGLE THREAD" in time_result):
            for i in range(5):
                single_thread_factor.append(
                    time_result["MKL SINGLE THREAD"][i]/time_result["OURS SINGLE THREAD"][i])
                eigen_single_thread_factor.append(
                    time_result["EIGEN SINGLE THREAD"][i]/time_result["OURS SINGLE THREAD"][i])
            plt.plot(x_range, single_thread_factor,
                     label="Single Thread Speedup", color="#55E6C1")

        plt.legend()
        plt.yticks(list(range(0, 13)))
        # plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
        #                     hspace=0, wspace=0)
        ax = plt.axes()
        ax.set_xscale('log')
        # ax.xaxis.set_major_locator(ticker.FixedLocator(x_range))
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + ".pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

        # plot the eigen difference
        plt.plot(x_range, eigen_single_thread_factor,
                 label="Single Thread Speedup Over Eigen", color="#d63031")
        ax = plt.axes()
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + "_eigen.pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()

        plt.plot(x_range, pre_comp_time,
                 label="Pre Computation Time", color="#b8e994")
        ax = plt.axes()
        ax.set_xscale('log')
        ax.xaxis.set_minor_locator(ticker.FixedLocator([]))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        save_plot_name = "test_result_graphs/" + \
            (file_list[k].split('.')[0]) + "_pre_comp.pdf"
        print(save_plot_name)
        plt.savefig(save_plot_name, bbox_inches='tight',
                    pad_inches=0, dpi=200)
        plt.close()
