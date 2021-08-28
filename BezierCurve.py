from typing import Optional

import numpy as np


class BezierCurve:
    @staticmethod
    def get_quadratic_curve(points: np.ndarray, n: int, rate_limit: float = 1.0, right_part=False) -> np.ndarray:
        """
        Calculates points of a quadratic Bezier curve based on 3 points and returns it as a numpy array
        :param rate_limit:
        :param points: 3 points Bezier curve based on
        :param n: number of points in between p1 and p3
        :return: array of Bezier curve points
        """
        if points.shape != (3, 2):
            return None
        t = np.linspace(rate_limit, 1.0, n).reshape((1, n)) if right_part else \
            np.linspace(0.0, rate_limit, n).reshape((1, n))

        point1, point2, point3 = points.reshape((3, 2, 1))
        result = point1.dot((1 - t) ** 2) + 2 * point2.dot(t * (1 - t)) + point3.dot(t ** 2)
        return result.T

    @staticmethod
    def get_quadratic_curve_inc(points: np.ndarray, n: int, rate=None, limited=False, right_part=False) -> np.ndarray:
        """
        Calculates points of a quadratic Bezier curve based on 3 points, but curve includes second point and
        returns it as a numpy array
        :param right_part:
        :param limited: should curve be built only to selected/calculated rate
        :param points: 3 points Bezier curve based on
        :param n: number of points in between p1 and p3
        :return: array of Bezier curve points
        :param rate: at which point of a curve should point 2 be. If None - calculated with distances between points
        """
        if points.shape != (3, 2):
            return None
        point1, point2, point3 = points.reshape((3, 2, 1))
        if rate is None:
            distance12sq = ((point2 - point1) ** 2).sum()
            distance23sq = ((point3 - point2) ** 2).sum()
            rate = 1 / (1 + np.sqrt(distance23sq / distance12sq))  # same as dist12/(dist12+dist23)
        point2c = (point2 - point1 * (1 - rate) ** 2 - point3 * rate ** 2) / (2 * rate * (1 - rate))
        return BezierCurve.get_quadratic_curve(np.array([point1.T[0], point2c.T[0], point3.T[0]]), n,
                                               rate if limited else 1.0, right_part)

    @staticmethod
    def get_quadratic_curve_multiple_points(points: np.ndarray, n: int) -> np.ndarray:
        """
        Builds a Beizer curve based on many quadratic Beizer curves between every 3 neighbour points
        (second point included in a curve) and interpolates these curves
        :param points:
        :param n:
        :return:
        """
        if points.shape[0] < 2:
            return None
        result = []
        prev = None
        t = np.linspace(0.0, 1.0, n).reshape((n, 1))
        if points.shape[0] == 2:
            return np.linspace(points[0], points[1], n)
        for i in range(0, points.shape[0] - 2):
            added = BezierCurve.get_quadratic_curve_inc(points[i:i + 3], n, limited=True)
            if prev is None:
                prev = added
            added_corrected = (1 - t) * prev + t * added
            prev = BezierCurve.get_quadratic_curve_inc(points[i:i + 3], n, limited=True, right_part=True)
            result.append(added_corrected)
        result.append(BezierCurve.get_quadratic_curve_inc(points[-3:], n, limited=True, right_part=True))
        return np.concatenate(result, axis=0)

    @staticmethod
    def get_cubic_curve(points: np.ndarray, n: int, rate_limit_1: float = 1.0, rate_limit_2: float = 1.0,
                        part_to_get: Optional[int] = 0):
        """

        :param points:
        :param n:
        :param rate_limit_1:
        :param rate_limit_2:
        :param part_to_get:
        :return:
        """
        if points.shape != (4, 2):
            return None

        t = (np.linspace(0.0, 1.0, n) if part_to_get is None else
             np.linspace(0.0, rate_limit_1, n) if part_to_get == 0 else
             np.linspace(rate_limit_1, rate_limit_2, n) if part_to_get == 1 else
             np.linspace(rate_limit_2, 1.0, n)).reshape((1, n))

        point1, point2, point3, point4 = points.reshape((4, 2, 1))
        result = point1.dot((1 - t) ** 3) + 3 * point2.dot(t * (1 - t) ** 2) + 3 * point3.dot(t ** 2 * (1 - t)) + point4.dot(t ** 3)
        return result.T

    @staticmethod
    def get_cubic_curve_inc(points: np.ndarray, n: int, rate_1: float = None, rate_2: float = None,
                        part_to_get: int = None):
        """

        :param points:
        :param n:
        :param rate_1:
        :param rate_2:
        :param part_to_get:
        :return:
        """
        if points.shape != (4, 2):
            return None
        p1, p2, p3, p4 = points.reshape((4, 2, 1))
        if rate_1 is None or rate_2 is None:
            distance12sq = ((p2 - p1) ** 2).sum()
            distance23sq = ((p3 - p2) ** 2).sum()
            distance34sq = ((p4 - p3) ** 2).sum()
            rate_1 = 1 / (1 + np.sqrt(distance23sq / distance12sq) + np.sqrt(distance34sq / distance12sq))
            rate_2 = 1 / (1 + 1 / (np.sqrt(distance12sq / distance34sq) + np.sqrt(distance23sq / distance34sq)))

        t1, t2, s1, s2 = rate_1, rate_2, 1-rate_1, 1-rate_2

        point3c = ((p1*s1**3*s2**2*t2 - p1*s1**2*s2**3*t1 - p2*s2**2*t2 + p3*s1**2*t1 - p4*s1**2*t1*t2**3 + p4*s2**2*t1**3*t2)/(3*s1*s2*t1*t2*(s1*t2 - s2*t1)))
        point2c = ((-p1*s1**3*s2*t2**2 + p1*s1*s2**3*t1**2 + p2*s2*t2**2 - p3*s1*t1**2 + p4*s1*t1**2*t2**3 - p4*s2*t1**3*t2**2)/(3*s1*s2*t1*t2*(s1*t2 - s2*t1)))
        return BezierCurve.get_cubic_curve(np.array([p1.T[0], point2c.T[0], point3c.T[0], p4.T[0]]),
                                           n, rate_1, rate_2, part_to_get)

    @staticmethod
    def get_cubic_curve_multiple_points(points: np.ndarray, n: int) -> np.ndarray:
        """
        Builds a Beizer curve based on many quadratic Beizer curves between every 3 neighbour points
        (second point included in a curve) and interpolates these curves
        :param points:
        :param n:
        :return:
        """
        if points.shape[0] < 2:
            return None
        result = []
        prev = None
        preprev = None
        t = np.linspace(0.0, 1.0, n).reshape((n, 1))
        if points.shape[0] == 2:
            return np.linspace(points[0], points[1], n)
        if points.shape[0] == 3:
            return BezierCurve.get_quadratic_curve_inc(points, n)
        for i in range(0, points.shape[0] - 1):
            added = BezierCurve.get_cubic_curve_inc(points[i:i + 4], n, part_to_get=0)
            if added is None:
                if prev is None:
                    prev = preprev
                added = prev
            if prev is None:
                prev = added
            if preprev is None:
                preprev = prev
            added_corrected = 0.5 * ((1 - t) * preprev + t * added) + 0.5 * prev
            prev = BezierCurve.get_cubic_curve_inc(points[i:i + 4], n, part_to_get=1)
            preprev = BezierCurve.get_cubic_curve_inc(points[i-1:i + 3], n, part_to_get=2)
            result.append(added_corrected)
        return np.concatenate(result, axis=0)