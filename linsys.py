from decimal import Decimal, getcontext
# from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):
    '''Linear System'''

    PLANES_SAME_DIM_MSG = 'Planes in the system should live in the same dim.'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.PLANES_SAME_DIM_MSG)

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x
        except AssertionError:
            raise Exception(self.PLANES_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = []
        for i, p in enumerate(self.planes):
            temp.append('Equation {}: {}'.format(i+1, p))
        ret += '\n'.join(temp)
        return ret

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        # num_variables = self.dimension

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coords)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def swap_rows(self, row1, row2):
        '''swap row1 with row2'''
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        '''mulitply row with coefficient'''
        n = self[row].normal_vector
        k = self[row].constant_term
        new_normal_vector = n.times_scalar(coefficient)
        new_constant_term = k*coefficient
        self[row] = Plane(new_normal_vector, new_constant_term)

    def add_multiple_times_row_to_row(self, coefficient, to_add, be_added_to):
        '''add multiple times row to row'''
        n1 = self[to_add].normal_vector
        n2 = self[be_added_to].normal_vector
        k1 = self[to_add].constant_term
        k2 = self[be_added_to].constant_term

        new_normal_vector = n1.times_scalar(coefficient).plus(n2)
        new_constant_term = (k1*coefficient) + k2

        self[be_added_to] = Plane(new_normal_vector, new_constant_term)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

if __name__ == '__main__':
    # p0 = Plane(Vector(['1', '1', '1']), '1')
    # p1 = Plane(Vector(['0', '1', '0']), '2')
    # p2 = Plane(Vector(['1', '1', '-1']), '3')
    # p3 = Plane(Vector(['1', '0', '-2']), '2')

    # s = LinearSystem([p0, p1, p2, p3])

    # print(s.indices_of_first_nonzero_terms_in_each_row())
    # print('{},{},{},{}'.format(s[0], s[1], s[2], s[3]))
    # print(len(s))
    # print(s)

    # s[0] = p1
    # print(s)

    # print(MyDecimal('1e-9').is_near_zero())
    # print(MyDecimal('1e-11').is_near_zero())

    p0 = Plane(Vector(['1', '1', '1']), '1')
    p1 = Plane(Vector(['0', '1', '0']), '2')
    p2 = Plane(Vector(['1', '1', '-1']), '3')
    p3 = Plane(Vector(['1', '0', '-2']), '2')

    s = LinearSystem([p0, p1, p2, p3])

    s.swap_rows(0, 1)

    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 1 failed')

    s.swap_rows(1, 3)
    if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
        print('test case 2 failed')

    s.swap_rows(3, 1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 3 failed')

    s.multiply_coefficient_and_row(1, 0)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 4 failed')

    s.multiply_coefficient_and_row(-1, 2)
    if not (s[0] == p1 and
            s[1] == p0 and
            s[2] == Plane(Vector(['-1', '-1', '1']), '-3') and
            s[3] == p3):
        print('test case 5 failed')

    s.multiply_coefficient_and_row(10, 1)
    if not (s[0] == p1 and
            s[1] == Plane(Vector(['10', '10', '10']), '10') and
            s[2] == Plane(Vector(['-1', '-1', '1']), '-3') and
            s[3] == p3):
        print('test case 6 failed')

    s.add_multiple_times_row_to_row(0, 0, 1)
    if not (s[0] == p1 and
            s[1] == Plane(Vector(['10', '10', '10']), '10') and
            s[2] == Plane(Vector(['-1', '-1', '1']), '-3') and
            s[3] == p3):
        print('test case 7 failed')

    s.add_multiple_times_row_to_row(1, 0, 1)
    if not (s[0] == p1 and
            s[1] == Plane(Vector(['10', '11', '10']), '12') and
            s[2] == Plane(Vector(['-1', '-1', '1']), '-3') and
            s[3] == p3):
        print('test case 8 failed')

    s.add_multiple_times_row_to_row(-1, 1, 0)
    if not (s[0] == Plane(Vector(['-10', '-10', '-10']), '-10') and
            s[1] == Plane(Vector(['10', '11', '10']), '12') and
            s[2] == Plane(Vector(['-1', '-1', '1']), '-3') and
            s[3] == p3):
        print('test case 9 failed')
