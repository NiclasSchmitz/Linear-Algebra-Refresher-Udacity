from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

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
        '''calculate pivot variable indices'''
        num_equations = len(self)

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

    def add_multiple_times_row_to_row(self, coefficient, from_here, to_here):
        '''add multiple times from_here row to_here row'''
        n1 = self[from_here].normal_vector
        n2 = self[to_here].normal_vector
        k1 = self[from_here].constant_term
        k2 = self[to_here].constant_term

        new_normal_vector = n1.times_scalar(coefficient).plus(n2)
        new_constant_term = (k1*coefficient) + k2
        self[to_here] = Plane(new_normal_vector, new_constant_term)

    def compute_triangular_form(self):
        '''form linear equations into echelon form

        Assumptions for Test Cases:
        1. Swap with topmost row below current row
        2. Don't multiply rows by numbers
        3. Only add a multiple of a row to rows underneath
        '''
        system = deepcopy(self)
        num_equations = len(system)
        num_variables = system.dimension
        j = 0  # current variable
        for i in range(num_equations):
            while j < num_variables:
                c = MyDecimal(system[i][j])
                if c.is_near_zero():
                    swap_succeeded = system.swap_with_row_below(i, j)
                    if not swap_succeeded:
                        j += 1
                        continue
                system.clear_coefficients_below(i, j)
                j += 1
                break
        return system

    def swap_with_row_below(self, row, col):
        '''find row with a value of != 0 in column col and swap'''
        num_equations = len(self)

        for k in range(row+1, num_equations):
            coefficient = MyDecimal(self[k][col])
            if not coefficient.is_near_zero():
                self.swap_rows(row, k)
                return True
        return False

    def clear_coefficients_below(self, row, col):
        '''all coefficients in column col below this row are getting cleared'''
        num_equations = len(self)
        beta = MyDecimal(self[row][col])

        for k in range(row+1, num_equations):
            n = self[k].normal_vector
            gamma = n[col]
            alpha = -gamma/beta
            self.add_multiple_times_row_to_row(alpha, row, k)

    def compute_rref(self):
        '''reduced row echelon form'''
        tf = self.compute_triangular_form()

        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()

        for i in range(num_equations)[::-1]:
            j = pivot_indices[i]
            if j < 0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(i, j)
            tf.clear_coefficients_above(i, j)
        return tf

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0')/n[col]
        self.multiply_coefficient_and_row(beta, row)

    def clear_coefficients_above(self, row, col):
        '''all coefficients in column col above this row are getting cleared'''
        for k in range(row)[::-1]:
            n = self[k].normal_vector
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, k)


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

    print('###########################')
    print('Quiz: Coding Row Operations')

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

    print('############################')
    print('Quiz: Coding Triangular Form')

    p1 = Plane(Vector(['1', '1', '1']), '1')
    p2 = Plane(Vector(['0', '1', '1']), '2')
    s = LinearSystem([p1, p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2):
        print('test case 1 failed')

    p1 = Plane(Vector(['1', '1', '1']), '1')
    p2 = Plane(Vector(['1', '1', '1']), '2')
    s = LinearSystem([p1, p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == Plane(constant_term='1')):
        print('test case 2 failed')

    p1 = Plane(Vector(['1', '1', '1']), '1')
    p2 = Plane(Vector(['0', '1', '0']), '2')
    p3 = Plane(Vector(['1', '1', '-1']), '3')
    p4 = Plane(Vector(['1', '0', '-2']), '2')
    s = LinearSystem([p1, p2, p3, p4])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2 and
            t[2] == Plane(Vector(['0', '0', '-2']), '2') and
            t[3] == Plane()):
        print('test case 3 failed')

    p1 = Plane(Vector(['0', '1', '1']), '1')
    p2 = Plane(Vector(['1', '-1', '1']), '2')
    p3 = Plane(Vector(['1', '2', '-5']), '3')
    s = LinearSystem([p1, p2, p3])
    t = s.compute_triangular_form()
    if not (t[0] == Plane(Vector(['1', '-1', '1']), '2') and
            t[1] == Plane(Vector(['0', '1', '1']), '1') and
            t[2] == Plane(Vector(['0', '0', '-9']), '-2')):
        print('test case 4 failed')

    print('#################')
    print('Quiz: Coding RREF')

    p1 = Plane(Vector(['1', '1', '1']), '1')
    p2 = Plane(Vector(['0', '1', '1']), '2')
    s = LinearSystem([p1, p2])
    r = s.compute_rref()
    if not (r[0] == Plane(Vector(['1', '0', '0']), '-1') and
            r[1] == p2):
        print('test case 1 failed')

    p1 = Plane(Vector(['1', '1', '1']), '1')
    p2 = Plane(Vector(['1', '1', '1']), '2')
    s = LinearSystem([p1, p2])
    r = s.compute_rref()
    if not (r[0] == p1 and
            r[1] == Plane(constant_term='1')):
        print('test case 2 failed')

    p1 = Plane(Vector(['1', '1', '1']), '1')
    p2 = Plane(Vector(['0', '1', '0']), '2')
    p3 = Plane(Vector(['1', '1', '-1']), '3')
    p4 = Plane(Vector(['1', '0', '-2']), '2')
    s = LinearSystem([p1, p2, p3, p4])
    r = s.compute_rref()
    if not (r[0] == Plane(Vector(['1', '0', '0']), '0') and
            r[1] == p2 and
            r[2] == Plane(Vector(['0', '0', '-2']), '2') and
            r[3] == Plane()):
        print('test case 3 failed')

    p1 = Plane(Vector(['0', '1', '1']), '1')
    p2 = Plane(Vector(['1', '-1', '1']), '2')
    p3 = Plane(Vector(['1', '2', '-5']), '3')
    s = LinearSystem([p1, p2, p3])
    r = s.compute_rref()
    if not (r[0] == Plane(Vector(['1', '0', '0']),
            Decimal('23')/Decimal('9')) and
            r[1] == Plane(Vector(['0', '1', '0']),
            Decimal('7')/Decimal('9')) and
            r[2] == Plane(Vector(['0', '0', '1']),
            Decimal('2')/Decimal('9'))):
        print('test case 4 failed')
