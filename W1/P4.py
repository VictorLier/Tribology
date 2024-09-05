import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('W1\Disc_MGF_virgin_R40_pos2_111031_P1167152.asc')
time = data[:, 0]
height = data[:, 1]

# a - plot
plt.figure()
plt.title("Height over time")
plt.plot(time, height)
plt.xlabel("Time [s]")
plt.ylabel("Height [mu]")


# b - Normalize
a, b = np.polyfit(time, height, 1) # polyfit

if True:
    plt.plot(time, a*time + b, label="Linear fit")
    plt.legend()

# remove linear trend
height = height - (a*time + b)

# Mean
mean = np.mean(height)*np.ones(len(height))

# Ten point
sorted_height = np.sort(height)
five_lowest = sorted_height[:5]
five_highest = sorted_height[-5:]

ten_mean = np.mean(np.concatenate((five_lowest, five_highest)))
ten_point = ten_mean*np.ones(len(height))



plt.figure()
plt.title("Height over time - Normalized")
plt.plot(time, height, label="Raw")
plt.plot(time, mean, label="Mean")
plt.plot(time, ten_point, label="Ten point")
# plt.plot(time, least_square, label="Least square")
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Height [mu]")
# plt.show()



# C - Roughness
Ra = np.mean(height) # Ra
Rq = np.sqrt(np.mean(height**2)) # Rq
Rt = np.max(height) - np.min(height) # Rt

print(f"Ra: {Ra:.2f} mu")
print(f"Rq: {Rq:.2f} mu")
print(f"Rt: {Rt:.2f} mu")

# D - Does it match turning?


print("Stop")