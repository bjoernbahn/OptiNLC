import numpy as np
import matplotlib.pyplot as plt

def asImage(input_file, output_file):

    data = np.loadtxt(input_file)
    t = data[:, 0]

    plt.figure(figsize=(10, 6))
    plt.plot(t, data[:, 1], label='Value 1')
    plt.plot(t, data[:, 2], label='Value 2')
    plt.plot(t, data[:, 3], label='Value 3')
    plt.plot(t, data[:, 4] , label='Value 4')

    plt.xlabel("Time (0.09 s steps)")
    plt.ylabel("Values")
    plt.title("Time series")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    #plt.show()
    plt.savefig(output_file)
    plt.close()

if __name__ == "__main__":
    asImage("eigen_data_max.txt", "max.png")
    asImage("eigen_data_min.txt", "min.png")