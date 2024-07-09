import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Pendulo:
    def __init__(self, n=8, thetas=None, thetaDots=None, g=-9.8):
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
                thetaDot_squared = np.clip(thetaDots[j], -100, 100) ** 2
                b_i -= (self.n - max(i, j)) * np.sin(thetas[i] - thetas[j]) * thetaDot_squared
            b_i -= self.g * (self.n - i) * np.sin(thetas[i])
            v[i] = b_i
        return v


    def f(self, thetas, thetaDots):
        A = self.A(thetas)
        b = self.b(thetas, thetaDots)
        return [thetaDots, np.linalg.solve(A, b)]

    def euler_explicito(self, dt):
        thetas_next = np.zeros(self.n)
        thetaDots_next = np.zeros(self.n)
        for i in range(self.n):
            f_i = self.f(self.thetas, self.thetaDots)  # Corrigido aqui
            thetas_next[i] = self.thetas[i] + self.thetaDots[i] * dt
            thetaDots_next[i] = self.thetaDots[i] + f_i[1][i] * dt  # Corrigido aqui
        self.thetas = thetas_next
        self.thetaDots = thetaDots_next

    def tick(self, dt):
        self.euler_explicito(dt)

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


