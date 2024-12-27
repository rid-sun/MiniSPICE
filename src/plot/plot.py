import matplotlib.pyplot as plt
import numpy as np
import random
import math

def visualize(results, timepoints, labels):
    #
    y_axis_data = np.array(results)
    x_axis_data = np.array(timepoints)

    # 信息输出
    with open("curves.txt", "w") as file:
        file.write("timepoints")
        for i in labels:
            file.write(f",x_{i}")
        file.write("\n")
        for i in range(len(x_axis_data)):
            combined = np.concatenate((x_axis_data[i:i+1], y_axis_data[i]))
            np.savetxt(file, [combined], delimiter=",", fmt="%f")

    #
    colors = [
        (random.random(), random.random(), random.random()) 
        for _ in range(len(y_axis_data[0]))
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    for idx, color in enumerate(colors):
        ax.plot(x_axis_data, y_axis_data[:, idx], label=f"x_{labels[idx]}", color=color)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.15),
        ncol=math.floor(math.sqrt(len(y_axis_data))),
        handleheight=1.5,
        columnspacing=1.0,
        fontsize='medium'
    )
    
    plt.savefig('ptran.png')
