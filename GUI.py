from math import cos, sin, atan2
import pygame
from pygame import gfxdraw
from BezierCurve import BezierCurve
from numpy import array


class GUI:
    def __init__(self, width=600, height=400):
        pygame.init()
        self.width = width
        self.height = height
        self.screen = pygame.display.set_mode((self.width, self.height + 40), pygame.RESIZABLE)
        pygame.display.set_caption('Beizer Curve Drawer')
        self.background = pygame.Surface(self.screen.get_size())
        # fill the background white
        self.background.fill((255, 255, 255))
        self.background = self.background.convert()  # faster blitting
        self.screen.set_colorkey((255, 255, 255))
        self.surface = pygame.Surface(self.screen.get_size())
        self.surface.fill((255, 255, 255))
        self.points = []
        self.running = False
        self.FPS = 60
        self.delta = 5  # max distance to point to delete/move it
        self.curve_function_idx = 0
        self.curve_functions = [BezierCurve.get_quadratic_curve_multiple_points,
                                BezierCurve.get_cubic_curve_multiple_points]

    def run(self):
        self.running = True
        clock = pygame.time.Clock()
        while self.running:
            milliseconds = clock.tick(self.FPS)
            self.update()



    def update(self):
        self.handle_events()
        self.render()
        self.blit()

    def handle_events(self):
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                self.running = False  # pygame window closed by user
            if event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:
                    if [event.pos] != self.points[-1:]:
                        self.points.append(event.pos)
                if event.button == 3:
                    closest_idx = -1
                    closest_dist_sq = -1
                    for i, point in enumerate(self.points):
                        dist_sq = ((array(event.pos) - point) ** 2).sum()
                        if dist_sq < self.delta ** 2:
                            if dist_sq < closest_dist_sq or closest_dist_sq == -1:
                                closest_dist_sq = dist_sq
                                closest_idx = i
                    if closest_idx != -1:
                        self.points.pop(closest_idx)

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.curve_function_idx = (self.curve_function_idx + 1) % 2

    def render(self):
        self.surface = pygame.Surface(self.screen.get_size())
        self.surface.fill((255,)*3)
        points = array(self.points)
        curve = self.curve_functions[self.curve_function_idx](points, 60)
        if curve is not None:
            pygame.draw.aalines(self.surface, (0, 0, 0), False, curve, 1)
        for x,  y in points:
            gfxdraw.filled_circle(self.surface, x, y, 2, (255, 0, 0))

    def blit(self):
        self.screen.blit(self.surface, (0, 0))
        pygame.display.flip()

    def draw_line(self, point1, point2, thickness):
        center = (point1 + point2) / 2
        length = ((point2 - point1) ** 2).sum() ** 0.5  # Total length of line
        angle = atan2(point1[1] - point2[1], point1[0] - point2[0])
        UL = (center[0] + (length / 2.) * cos(angle) - (thickness / 2.) * sin(angle),
              center[1] + (thickness / 2.) * cos(angle) + (length / 2.) * sin(angle))
        UR = (center[0] - (length / 2.) * cos(angle) - (thickness / 2.) * sin(angle),
              center[1] + (thickness / 2.) * cos(angle) - (length / 2.) * sin(angle))
        BL = (center[0] + (length / 2.) * cos(angle) + (thickness / 2.) * sin(angle),
              center[1] - (thickness / 2.) * cos(angle) + (length / 2.) * sin(angle))
        BR = (center[0] - (length / 2.) * cos(angle) + (thickness / 2.) * sin(angle),
              center[1] - (thickness / 2.) * cos(angle) - (length / 2.) * sin(angle))
        pygame.gfxdraw.aapolygon(self.surface, (UL, UR, BR, BL), (0, 0, 0))
        pygame.gfxdraw.filled_polygon(self.surface, (UL, UR, BR, BL), (0, 0, 0))

