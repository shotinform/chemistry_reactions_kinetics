import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def vector_field(x, t):
    return x * (1 - x)

def plot_vector_field_and_curves():
    x = np.linspace(-2, 2, 20)
    y = np.linspace(-2, 2, 20)
    X, Y = np.meshgrid(x, y)
    U = X * (1 - X)
    V = np.zeros_like(U)  # Bez promene u y pravcu

    fig, ax = plt.subplots()
    ax.quiver(X, Y, U, V, scale=20)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Vektorsko polje i integralne krive')

    for x_start in np.linspace(-2, 2, 5):
        t_values = np.linspace(0, 2, 20)
        for t_start in t_values:
            trajectory, info_dict = integrate_curve([x_start, 0], t_start, full_output=True)
            ax.plot(trajectory[:, 0], trajectory[:, 1], color='b', alpha=0.3)

    plt.show()

def integrate_curve(x_start, t_start, full_output=False):
    def dx_dt(x, t):
        return vector_field(x, t)

    t = np.linspace(t_start, t_start + 2, 100)
    trajectory, info_dict = odeint(dx_dt, x_start, t, full_output=full_output)
    if full_output:
        return trajectory, info_dict
    else:
        return trajectory

if __name__ == "__main__":
    plot_vector_field_and_curves()
