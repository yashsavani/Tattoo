#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Tattoo Generator

Uses points drawn from a dirichlet distribution to generate a delaunay triangulation of a triangle
i.e. makes divides a triangle up into 11 smaller random triangles.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.special import huber
from tqdm import tqdm


def angle(vec1, vec2):
    """Angle between numpy vectors"""
    normalize = lambda vec: vec / np.linalg.norm(vec)
    return np.arccos(np.clip(np.dot(normalize(vec1), normalize(vec2)), -1.0, 1.0))

def volume(vec1, vec2):
    """Find volume of triangle using any two side vectors

    Uses half of the determinant of the matrix constructed with the side vectors as columns
    """
    return np.abs(0.5*np.linalg.det(np.array([vec1, vec2])))

def tri_distance(pnts, triangles, all_vols, all_angles, omega=0.2):
    """calculates the area of all the triangles in the Delaunay triangulation."""

    vols = []
    for triangle in triangles:
        points = [pnts[triangle[i]] for i in range(3)]

        vec1 = points[1] - points[0]
        vec2 = points[2] - points[0]
        vec3 = points[2] - points[1]
        angles = np.array([angle(vec1, vec2), angle(vec3, -vec1), angle(-vec2, -vec3)])
        angle_loss = np.mean(huber(1, (angles-1.047)/0.864))
        all_angles.append(angle_loss)

        # Use log because the distribution is right skewed
        normalized_vol = (np.log(volume(vec1, vec2)) + 3.40) / 0.62
        vols.append(omega*normalized_vol**2 + (1.0-omega)*angle_loss)
        all_vols.append(normalized_vol)

    return np.mean(vols)

def get_triangulations(iterations):
    """Get possible triangulations of a triangle sampled from a dirichlet distribution"""
    all_vols = []
    all_angles = []
    triangulations = []

    corners = np.array([[0, 0], [1, 0], [0.5, 0.75**0.5]])
    for _ in tqdm(range(iterations)):
        points = np.random.dirichlet([2]*3, 5).round(2)
        all_points = np.vstack([points, np.array([1, 0, 0]),
                                np.array([0, 1, 0]), np.array([0, 0, 1])])

        # Skip examples where points are drawn on top of each other
        euclidean_dists = 2*(np.sum(points**2, axis=1) - np.dot(points, points.T))
        if np.sum(euclidean_dists < 0.01) > points.shape[0]:
            continue

        # Skip if any points drawn outside the edges of the triangle (rounding error)
        if np.any(np.sum(points, axis=1) > 1):
            continue

        # Barycentric interpolation from points to planar coordinates
        planar_points = np.array([np.dot(p, corners) for p in all_points])

        # Create a Delaunay triangulation from the points
        tri = Delaunay(planar_points)
        omega = 0.2 # tested empirically
        dist = tri_distance(planar_points, tri.simplices.copy(), all_vols, all_angles, omega=omega)
        triangulations.append((planar_points, all_points, omega, dist))

    return (all_vols, all_angles, triangulations)

def plot_stats(var, title):
    """print statistics and plot histogram for a variable"""
    print title
    print 'mean = {}'.format(np.mean(var))
    print 'std = {}'.format(np.std(var))
    print 'median = {}'.format(np.median(var))
    plt.style.use('ggplot')
    plt.title(title)
    plt.hist(var, bins=100)
    plt.show()


def main():
    """Main Function"""


    all_vols, all_angles, triangulations = get_triangulations(10000)

    sort_tri = sorted(triangulations, key=lambda x: x[-1])

    # Uncomment if you want statistical information about the triangulation
    plot_stats(all_vols, 'Volume of triangles')
    plot_stats(all_angles, 'Angle loss for triangle')

    directory = 'figs'
    if not os.path.exists(directory):
        os.makedirs(directory)
    for iteration, (pnts, tri_pnts, omega, dist) in enumerate(sort_tri[:10]):
        print 'dist = {}'.format(dist)
        print tri_pnts
        print 'omega = {}'.format(omega)
        tri = Delaunay(pnts)
        plt.triplot(pnts[:, 0], pnts[:, 1], tri.simplices.copy(), color='black')
        # plt.plot(pnts[:, 0], pnts[:, 1], 'o')

        plt.axis('off')
        plt.axis('equal')
        plt.savefig('{}/{}.svg'.format(directory, iteration), format='svg')
        plt.show()

if __name__ == "__main__":
    main()
