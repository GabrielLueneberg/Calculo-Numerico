import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Pendulo:
    def __init__(self, n=3, thetas=None, thetaDots=None, g=-9.8):
        self.n = n
        self.thetas = np.full(n, 0.5 * np.pi) if thetas is None else thetas
        self.thetaDots = np.zeros(n) if thetaDots is None else thetaDots
        self.g = -g

    def A(self, thetas):
        M = np.zeros((self.n, self.n))
        for i in range(self.n):
            for j in range(self.n):
                M[i, j] = (self.n - max(i, j)) * np.cos(thetas[i] - thetas[j])
        return M

    def b(self, thetas, thetaDots):
        v = np.zeros(self.n)
        for i in range(self.n):
            b_i = 0
            for j in range(self.n):
                b_i -= (self.n - max(i, j)) * np.sin(thetas[i] - thetas[j]) * thetaDots[j] ** 2
            b_i -= self.g * (self.n - i) * np.sin(thetas[i])
            v[i] = b_i
        return v

    def f(self, thetas, thetaDots):
        A = self.A(thetas)
        b = self.b(thetas, thetaDots)
        return [thetaDots, np.linalg.solve(A, b)]

    def RK4(self, dt, thetas, thetaDots):
        k1 = self.f(thetas, thetaDots)
        k2 = self.f(thetas + 0.5 * dt * k1[0], thetaDots + 0.5 * dt * k1[1])
        k3 = self.f(thetas + 0.5 * dt * k2[0], thetaDots + 0.5 * dt * k2[1])
        k4 = self.f(thetas + dt * k3[0], thetaDots + dt * k3[1])

        thetaDeltas = (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * dt / 6
        thetaDotDeltas = (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * dt / 6

        return [thetas + thetaDeltas, thetaDots + thetaDotDeltas]

    def tick(self, dt):
        newState = self.RK4(dt, self.thetas, self.thetaDots)
        self.thetas = newState[0]
        self.thetaDots = newState[1]

    @property
    def coordinates(self):
        x = 0
        y = 0
        coords = []
        for i in range(len(self.thetas)):
            theta = self.thetas[i]
            x += np.sin(theta)
            y += np.cos(theta)
            coords.append((x, y))
        return coords

pendulo = Pendulo()

fig, ax = plt.subplots()

def draw(i):
    ax.clear()
    coords = pendulo.coordinates
    x1, y1 = 0.5, 0.5

    for coord in coords:
        x2 = 0.5 + coord[0] * 0.5 / pendulo.n
        y2 = 0.5 - coord[1] * 0.5 / pendulo.n
        ax.plot([x1, x2], [y1, y2], color='black', linewidth=2)
        ax.scatter(x2, y2, color='red')

        x1, y1 = x2, y2

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

def animate(i):
    draw(i)
    pendulo.tick(1/30)

ani = FuncAnimation(fig, animate, frames=100, interval=33)
plt.show()
