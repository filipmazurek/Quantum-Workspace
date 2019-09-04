from forest_utils_ms import *
import tsp_qaoa_forest
import sys


def main():

    cities = [[0, 0], [0, 1]]

    distance_matrix = get_distance_matrix(cities)

    print()
    print("Distance Matrix: ")
    print(distance_matrix)

    solver = tsp_qaoa_forest.ForestTSPSolver(distance_matrix, steps=3, use_constraints=True)

    solution, distribution = solver.solve_tsp()

    print()
    print("Solution: ")
    print(solution)

    print()
    print("Distribution: ")
    print(distribution)


main()
