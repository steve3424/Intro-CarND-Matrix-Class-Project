import math
from math import sqrt
import numbers
from copy import deepcopy

def zeroes(height, width):
        """
        Creates a matrix of zeroes.
        """
        g = [[0.0 for _ in range(width)] for __ in range(height)]
        return Matrix(g)

def identity(n):
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I

class Matrix(object):

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
 
    def determinant(self):
        """
        Calculates the determinant of any size matrix.
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate determinant of non-square matrix.")

        if self.h == 1:
            return self[0][0]
        elif self.h == 2:
            a = self[0][0]
            b = self[0][1]
            c = self[1][0]
            d = self[1][1]
            return a*d - b*c
        else:
            det = 0
            sign = 1
            for i,num in enumerate(self[0]):
                sub_matrix = self.get_sub_matrix(0, i)
                det = det + sign * num * sub_matrix.determinant()
                sign *= -1
            return det
    
    ## for determinant of matrix > 2x2
    def get_sub_matrix(self, row, col):
        M = deepcopy(self.g)
        del M[row]
        for i in range(len(M)):
            del M[i][col]
        return Matrix(M)

    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries).
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate the trace of a non-square matrix.")
        
        return sum([self[i][i] for i in range(self.h)])

    def inverse(self):
        """
        Calculates the inverse of any size matrix.
        """
        if not self.is_square():
            raise(ValueError, "Non-square Matrix does not have an inverse.")
            
        if self.h == 1:
            M = [[1/self[0][0]]]
            return Matrix(M)
        elif self.h == 2:
            tr_self = self.trace()
            I = identity(self.h)
            M = tr_self * I - self
            det = 1 / self.determinant()
            return det * M
        else:
            ## uses adjugate method for matrix > 2x2
            M = self.matrix_of_cofactors()
            adjugate = M.T()
            det = self.determinant()
            inverse = 1/det * adjugate
            return inverse
                
    ## for inverse of matrix > 2x2
    def matrix_of_cofactors(self):
        '''
        Creates matrix of minors and switches signs to return matrix of cofactors
        '''
        M = []
        sign = 1
        for i in range(self.h):
            new_row = []
            for j in range(self.w):
                sub_m = self.get_sub_matrix(i, j)
                new_row.append(sign * sub_m.determinant())
                sign *= -1
                
                ## for even widths
                if self.w % 2 == 0 and j == self.w - 1:
                    sign *= -1
            M.append(new_row)
        return Matrix(M)

    def T(self):
        """
        Returns a transposed copy of this Matrix.
        """
        M = []
        for j in range(self.w):
            new_row = [self[i][j] for i in range(self.h)]
            M.append(new_row)
        return Matrix(M)

    def is_square(self):
        return self.h == self.w

    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be added if the dimensions are the same")
        M = []
        for i in range(self.h):
            new_row = [x + y for x,y in zip(self[i], other[i])]
            M.append(new_row)
        return Matrix(M)

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        M = []
        for row in self:
            new_row = [-1 * col for col in row]
            M.append(new_row)
        return Matrix(M)

    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be subtracted if the dimensions are the same") 
            
        M = []
        for i in range(self.h):
            new_row = [x - y for x,y in zip(self[i], other[i])]
            M.append(new_row)
        return Matrix(M)

    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        if self.w != other.h:
            raise(ValueError, "Matrices must be of m*n and n*p dimensions to multiply!")
        
        M = []
        transpose = other.T()
        for i in range(self.h):
            new_row = []
            for j in range(transpose.h):
                dot_product = sum([x*y for x,y in zip(self[i], transpose[j])])
                new_row.append(dot_product)
            M.append(new_row)
        return Matrix(M)
            

    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.

        Example:

        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        
        if isinstance(other, numbers.Number):
            M = []
            for row in self:
                new_row = [other * row[i] for i in range(self.w)]
                M.append(new_row) 
            return Matrix(M)
            