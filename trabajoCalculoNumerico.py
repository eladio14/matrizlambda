
import numpy as np

def lu_decomposition(matrix):
    n = len(matrix)
    lower = np.zeros((n, n))
    upper = np.zeros((n, n))

    for i in range(n):
        # Upper Triangular
        for k in range(i, n):
            sum = 0
            for j in range(i):
                sum += (lower[i][j] * upper[j][k])
            upper[i][k] = matrix[i][k] - sum

        # Lower Triangular
        for k in range(i, n):
            if i == k:
                lower[i][i] = 1
            else:
                sum = 0
                for j in range(i):
                    sum += (lower[k][j] * upper[j][i])
                lower[k][i] = (matrix[k][i] - sum) / upper[i][i]

    return lower, upper

def solve_system(lower, upper, b):
    n = len(lower)
    y = np.zeros(n)
    x = np.zeros(n)

    # Solve Ly = b
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += lower[i][j] * y[j]
        y[i] = b[i] - sum

    # Solve Ux = y
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i+1, n):
            sum += upper[i][j] * x[j]
        x[i] = (y[i] - sum) / upper[i][i]

    return x

# Get input from user
n = int(input("Enter the size of the square matrix: "))
matrix = np.zeros((n, n))
b = np.zeros(n)

print("Enter the coefficient matrix:")
for i in range(n):
    matrix[i] = list(map(float, input().split()))

print("Enter the vector of independent terms:")
b = list(map(float, input().split()))

# Perform LU decomposition
lower, upper = lu_decomposition(matrix)

# Solve the system
