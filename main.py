from BezierCurve import BezierCurve
from GUI import GUI
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import keyboard as kb


def mpl_main():

    number = 8

    t = np.linspace(0, 10, number).reshape((number, 1))
    points = np.random.random((number, 2))*10
    vel = (np.random.random((number, 2))-0.5)*0.1
    fig, ax = plt.subplots()
    pl1, = plt.plot([], [], )
    pl2, = plt.plot([], [], )
    sc, = plt.plot([], [], '.')
    running = True
    def toggle():
        nonlocal running
        running = not running
    kb.add_hotkey('space', toggle)

    def init():
        ax.set_xlim(-1, 11)
        ax.set_ylim(-1, 11)
        return pl1, pl2, sc

    def update(frame):
        nonlocal points
        nonlocal vel
        if running:
            points += vel
        # vel += (np.random.random((8, 1)) - 0.5) * 0.05
        # points[:, 0] = 10 * (points[:, 0] - points[:, 0].min()) / (points[:, 0].max() - points[:, 0].min())
        # points[:, 1] = 10 * (points[:, 1] - points[:, 1].min()) / (points[:, 1].max() - points[:, 1].min())
        vel[points >= 10] = -abs(vel[points >= 10])
        vel[points <= 0] = abs(vel[points <= 0])
        curve1 = BezierCurve.get_cubic_curve_multiple_points(points, 60)
        curve2 = BezierCurve.get_quadratic_curve_multiple_points(points, 60)
        pl1.set_data(curve1[:, 0], curve1[:, 1])
        pl2.set_data(curve2[:, 0], curve2[:, 1])
        sc.set_data(points[:, 0], points[:, 1])
        # ax.set_xlim(curve[:, 0].min() - 1, curve[:, 0].max() + 1)
        # ax.set_ylim(curve[:, 1].min() - 1, curve[:, 1].max() + 1)
        ax.figure.canvas.draw()
        return pl1, pl2, sc

    ani = FuncAnimation(fig, update, init_func=init, blit=True, interval=10)
    plt.show()


def pg_main():
    gui = GUI()
    gui.run()


if __name__ == '__main__':
    pg_main()
